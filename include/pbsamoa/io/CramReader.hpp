#ifndef PBSAMOA_IO_CRAMREADER_HPP
#define PBSAMOA_IO_CRAMREADER_HPP

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/cram/CramReference.hpp>
#include <pbsamoa/index/CraiIndex.hpp>

#include <filesystem>
#include <memory>
#include <optional>
#include <ranges>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

struct CramReaderConfig
{
    /// \brief Optional FASTA file for decoding reference-based CRAM slices.
    ///
    /// The reader maps the file as a FastaFile: it reads the bases it needs, using `<fasta>.fai`
    /// when it exists. Set ReferencePath or Reference, not both.
    std::filesystem::path ReferencePath{};

    /// \brief Optional reference for decoding reference-based CRAM slices, given in windows.
    ///
    /// One sequence per @SQ line, in header order. The reader keeps the bases of one window
    /// only. Several readers can share one reference.
    std::shared_ptr<const CramReference> Reference{};

    /// \brief Number of worker threads that decode the containers of the sequential read.
    ///
    /// The workers decode up to 2 containers each ahead of the caller, and the reader gives
    /// their records in file order. Query and ReadSlice use the workers for the blocks of a
    /// slice. 0 decodes all on the calling thread.
    std::size_t DecompressionWorkers{0};
};

/// \brief Reads CRAM files and yields BamRecord objects.
///
/// Follows the same pattern as SamReader: pimpl, ReadRecord(), and
/// range interface via Records(). Decodes CRAM containers, slices,
/// and blocks per the CRAMv3 specification.
class CramReader
{
public:
    explicit CramReader(const std::filesystem::path& path,
                        const CramReaderConfig& config = CramReaderConfig{});
    ~CramReader();

    CramReader(const CramReader&) = delete;
    CramReader& operator=(const CramReader&) = delete;
    CramReader(CramReader&&) noexcept;
    CramReader& operator=(CramReader&&) noexcept;

    /// \brief Access the parsed SAM header from the CRAM header container.
    const SamHeader& Header() const;

    /// \brief Read the next record. Returns nullopt at EOF.
    std::optional<BamRecord> ReadRecord();

    /// \brief Read the next record as an owning RawRecord. Returns nullopt at
    /// EOF.
    std::optional<RawRecord> ReadRawRecord();

    /// \brief Return records overlapping [beg, end) from CRAI-indexed CRAM
    /// slices.
    ///
    /// Coordinates follow BAM query semantics: beg is 0-based inclusive, end is
    /// 0-based exclusive. For unmapped queries, use refId = -1.
    std::vector<BamRecord> Query(const CraiIndex& index, std::int32_t refId, std::int32_t beg,
                                 std::int32_t end);

    /// \brief Return raw BAM-layout records overlapping [beg, end) from
    /// CRAI-indexed CRAM slices.
    std::vector<RawRecord> QueryRaw(const CraiIndex& index, std::int32_t refId, std::int32_t beg,
                                    std::int32_t end);

    /// \brief All records of the slice of a CRAI row (from CraiIndex::Query), in file order.
    ///
    /// Like Query, it reads on the stream of the reader, and the sequential iteration goes on
    /// unchanged. A caller that keeps decoded slices queries each slice once.
    std::vector<BamRecord> ReadSlice(const CraiEntry& entry);

    /// \brief ReadSlice with raw BAM-layout records.
    std::vector<RawRecord> ReadSliceRaw(const CraiEntry& entry);

    /// \brief Range interface for iteration.
    class RecordRange : public std::ranges::view_interface<RecordRange>
    {
    public:
        class Iterator
        {
        public:
            Iterator();
            explicit Iterator(CramReader* reader);

            const BamRecord& operator*() const;
            const BamRecord* operator->() const;
            Iterator& operator++();
            void operator++(int);
            bool operator==(const Iterator& other) const;

        private:
            CramReader* reader_{nullptr};
            std::optional<BamRecord> current_;
        };

        explicit RecordRange(CramReader* reader);
        Iterator begin();
        Iterator end();

    private:
        CramReader* reader_;
    };

    RecordRange Records();

    /// \brief Range interface for iteration over RawRecord output.
    class RawRecordRange : public std::ranges::view_interface<RawRecordRange>
    {
    public:
        class Iterator
        {
        public:
            Iterator();
            explicit Iterator(CramReader* reader);

            const RawRecord& operator*() const;
            const RawRecord* operator->() const;
            Iterator& operator++();
            void operator++(int);
            bool operator==(const Iterator& other) const;

        private:
            CramReader* reader_{nullptr};
            std::optional<RawRecord> current_;
        };

        explicit RawRecordRange(CramReader* reader);
        Iterator begin();
        Iterator end();

    private:
        CramReader* reader_;
    };

    RawRecordRange RawRecords();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_IO_CRAMREADER_HPP
