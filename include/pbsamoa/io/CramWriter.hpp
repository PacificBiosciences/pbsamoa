#ifndef PBSAMOA_IO_CRAMWRITER_HPP
#define PBSAMOA_IO_CRAMWRITER_HPP

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/cram/CramStructs.hpp>

#include <filesystem>
#include <memory>
#include <optional>
#include <unordered_map>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

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
};

/// \brief Writes CRAM files from BamRecord and RawRecord inputs.
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

    /// \brief Write a single raw record view.
    void Write(const RawRecord& record);

    /// \brief Write a batch of raw records.
    void WriteBatch(const RawRecordBatch& batch);

    /// \brief Flush and close. Called automatically by destructor.
    void Close();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_IO_CRAMWRITER_HPP
