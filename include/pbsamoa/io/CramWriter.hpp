#ifndef PBSAMOA_IO_CRAMWRITER_HPP
#define PBSAMOA_IO_CRAMWRITER_HPP

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/cram/CramStructs.hpp>
#include <pbsamoa/index/CraiIndex.hpp>

#include <filesystem>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief The reference of a reference-based CRAM file, given in windows: for example from a
/// packed form, so the text of the whole reference need not be in memory.
///
/// Several threads call Bases at the same time.
class CramReference
{
public:
    virtual ~CramReference() = default;

    /// \brief The number of reference sequences: one per @SQ line, in header order.
    [[nodiscard]] virtual std::int32_t NumSequences() const = 0;

    /// \brief The length of reference sequence refId.
    [[nodiscard]] virtual std::int64_t Length(std::int32_t refId) const = 0;

    /// \brief The bases [start, start + length) of reference sequence refId (inside it).
    ///
    /// The result can be a view of buffer, which the call can change. The case does not
    /// matter.
    [[nodiscard]] virtual std::string_view Bases(std::int32_t refId, std::int64_t start,
                                                 std::int64_t length,
                                                 std::string& buffer) const = 0;
};

struct CramWriterConfig
{
    /// \brief Block compression method to use for data blocks.
    CramBlockMethod BlockCompressionMethod{CramBlockMethod::GZIP};

    /// \brief Optional per-data-series block compression overrides.
    ///
    /// Keys are CRAM data-series columns (for example RN, QS, BA). Values
    /// override BlockCompressionMethod for those columns only.
    std::unordered_map<CramDataSeries, CramBlockMethod> DataSeriesCompressionMethods{};

    /// \brief Maximum number of records per slice.
    std::int32_t RecordsPerSlice{10000};

    /// \brief Maximum number of slices per container.
    ///
    /// A value of 1 preserves the original one-slice-per-container behavior.
    /// A value <= 0 means "unlimited" and only flushes at Close().
    std::int32_t SlicesPerContainer{1};

    /// \brief CRAM major version (must be 3).
    std::uint8_t MajorVersion{3};

    /// \brief CRAM minor version (0 or 1).
    std::uint8_t MinorVersion{0};

    /// \brief Write a samtools-compatible CRAI sidecar alongside the CRAM output.
    bool WriteCrai{false};

    /// \brief Optional override path for CRAI output. Defaults to
    /// "<output>.crai".
    std::optional<std::filesystem::path> CraiPath{};

    /// \brief Number of worker threads for parallel CRAM block compression.
    ///
    /// 0 disables parallel compression and uses synchronous encoding.
    std::size_t CompressionWorkers{0};

    /// \brief Optional compression level for method-specific codecs.
    ///
    /// For gzip this maps to libdeflate levels [0, 12].
    std::optional<int> CompressionLevel{};

    /// \brief Write through a temporary file and atomically rename on Close().
    bool UseTempFile{false};

    /// \brief The reference sequence of each @SQ line, in header order, for reference-based
    /// CRAM: a mapped record stores only its differences to the reference, the slices carry
    /// the MD5 of their reference bases, and @SQ lines without M5 get one. Readers then need
    /// the reference. Empty (the default): the records store all their bases.
    ///
    /// The case of the reference does not matter. The views must stay valid until Close().
    std::vector<std::string_view> ReferenceSequences{};

    /// \brief The reference in windows, in place of ReferenceSequences (set one of the two).
    std::shared_ptr<const CramReference> Reference{};
};

/// \brief The M5 value of a reference sequence: the MD5 of its bases in upper case, as 32
/// lowercase hex digits (SAM spec).
[[nodiscard]] std::string CramReferenceMd5(std::string_view sequence);

/// \brief The M5 value of reference sequence refId, from windows of 64 KiB.
[[nodiscard]] std::string CramReferenceMd5(const CramReference& reference, std::int32_t refId);

/// \brief A container that CramWriter::EncodeContainer encoded, for CramWriter::WriteContainer.
struct CramEncodedContainer
{
    /// \brief The compressed slices; their record counters start at 0 until the write.
    CramContainer Container;

    /// \brief The CRAI entries of each slice, without their file offsets (empty without
    /// WriteCrai).
    std::vector<std::vector<CraiEntry>> SliceIndex;
};

/// \brief Writes CRAM files from BamRecord, RawRecordView, and RawRecord inputs.
///
/// Follows the BamWriter pattern: pimpl, Write(...), Close().
/// Encodes records into CRAM containers with configurable compression.
class CramWriter
{
public:
    CramWriter(const std::filesystem::path& path, const SamHeader& header,
               const CramWriterConfig& config = {});
    ~CramWriter();

    CramWriter(const CramWriter&) = delete;
    CramWriter& operator=(const CramWriter&) = delete;
    CramWriter(CramWriter&&) noexcept;
    CramWriter& operator=(CramWriter&&) noexcept;

    /// \brief Write a single record.
    void Write(const BamRecord& record);

    /// \brief Write a single record (move overload).
    void Write(BamRecord&& record);

    /// \brief Write a single non-owning raw record view.
    void Write(const RawRecordView& record);

    /// \brief Write a single owning raw record.
    void Write(const RawRecord& record);

    /// \brief Write a batch of raw records.
    void WriteBatch(const RawRecordBatch& batch);

    /// \brief Encode records into one container (slices of RecordsPerSlice records) for
    /// WriteContainer.
    ///
    /// Several threads can encode at the same time, also while another thread writes: this
    /// uses the calling thread only, not the compression workers.
    [[nodiscard]] CramEncodedContainer EncodeContainer(std::span<const BamRecord> records) const;

    /// \brief Write a container of EncodeContainer after all records written before.
    void WriteContainer(CramEncodedContainer container);

    /// \brief Flush and close. Called automatically by destructor.
    void Close();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_IO_CRAMWRITER_HPP
