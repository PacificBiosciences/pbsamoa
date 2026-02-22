#ifndef PBSAMOA_IO_BAMRAWREADER_HPP
#define PBSAMOA_IO_BAMRAWREADER_HPP

#include <pbsamoa/core/Bgzf.hpp>
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/index/BaiIndex.hpp>
#include <pbsamoa/index/ZmwWhitelist.hpp>

#include <filesystem>
#include <memory>
#include <optional>
#include <ranges>
#include <string_view>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

class ByteLimit
{
public:
    constexpr ByteLimit();
    constexpr explicit ByteLimit(std::size_t bytes);
    constexpr std::size_t Value() const;

private:
    std::size_t bytes_;
};

struct ThrowPolicy
{
    [[noreturn]] void OnCorruptRecord(std::string_view message) const;
};

class SkipPolicy
{
public:
    void OnCorruptRecord(std::string_view message);
    std::size_t SkippedCount() const;

private:
    std::size_t skipped_{0};
};

// --- inline constexpr definitions ---

constexpr ByteLimit::ByteLimit() : bytes_{std::size_t{256} * 1024 * 1024} {}

constexpr ByteLimit::ByteLimit(std::size_t bytes) : bytes_{bytes} {}

constexpr std::size_t ByteLimit::Value() const { return bytes_; }

namespace Literals {

constexpr ByteLimit operator""_MiB(unsigned long long value)
{
    return ByteLimit{value * 1024 * 1024};
}

constexpr ByteLimit operator""_KiB(unsigned long long value) { return ByteLimit{value * 1024}; }

}  // namespace Literals

struct BamRawReaderConfig
{
    std::size_t BgzfWorkers{0};
    std::size_t RecordLimit{0};   // 0 = unlimited; >0 = stop after N records
    std::int32_t ChunkNum{0};     // 0 = no chunking; 1-based chunk number
    std::int32_t TotalChunks{0};  // 0 = no chunking; total number of chunks
    std::optional<ZmwWhitelist> Whitelist{};
};

class BamRawReader
{
public:
    /// \brief Construct a BAM raw reader.
    /// Default config uses synchronous decompression; set BgzfWorkers > 0 for parallel pipeline.
    explicit BamRawReader(const std::filesystem::path& path, BamRawReaderConfig config = {});

    ~BamRawReader();

    BamRawReader(const BamRawReader&) = delete;
    BamRawReader& operator=(const BamRawReader&) = delete;
    BamRawReader(BamRawReader&&) noexcept;
    BamRawReader& operator=(BamRawReader&&) noexcept;

    const SamHeader& Header() const;

    /// Returns an owning record. Safe to hold across ReadRecord/ReadBatch calls.
    std::optional<RawRecord> ReadRecord();

    /// Read batch bounded by memory budget. Returns nullopt at EOF.
    std::optional<RawRecordBatch> ReadBatch(ByteLimit limit = ByteLimit{});

    void Seek(VirtualOffset offset);
    VirtualOffset Tell() const;

    /// Range interface: for (const auto& view : reader.Records()) { ... }
    class RecordRange : public std::ranges::view_interface<RecordRange>
    {
    public:
        class Iterator
        {
        public:
            using difference_type = std::ptrdiff_t;
            using value_type = RawRecord;
            using iterator_concept = std::input_iterator_tag;

            Iterator();
            explicit Iterator(BamRawReader* reader);

            const RawRecord& operator*() const;
            const RawRecord* operator->() const;
            Iterator& operator++();
            void operator++(int);
            bool operator==(const Iterator& other) const;

        private:
            BamRawReader* reader_{nullptr};
            std::optional<RawRecord> current_;
        };

        explicit RecordRange(BamRawReader* reader);
        Iterator begin();
        Iterator end();

    private:
        BamRawReader* reader_;
    };

    RecordRange Records();

    /// Range over records matching a BAI region query.
    class QueryRange : public std::ranges::view_interface<QueryRange>
    {
    public:
        class Iterator
        {
        public:
            using difference_type = std::ptrdiff_t;
            using value_type = RawRecord;
            using iterator_concept = std::input_iterator_tag;

            Iterator();
            Iterator(BamRawReader* reader, std::vector<Chunk> chunks, std::int32_t refId,
                     std::int32_t beg, std::int32_t end);

            const RawRecord& operator*() const;
            const RawRecord* operator->() const;
            Iterator& operator++();
            void operator++(int);
            bool operator==(const Iterator& other) const;

        private:
            BamRawReader* reader_{nullptr};
            std::vector<Chunk> chunks_;
            std::size_t chunkIdx_{0};
            std::int32_t refId_{-1};
            std::int32_t beg_{0};
            std::int32_t end_{0};
            std::optional<RawRecord> current_;

            void Advance();
        };

        QueryRange(BamRawReader* reader, std::vector<Chunk> chunks, std::int32_t refId,
                   std::int32_t beg, std::int32_t end);
        Iterator begin();
        Iterator end();

    private:
        BamRawReader* reader_;
        std::vector<Chunk> chunks_;
        std::int32_t refId_;
        std::int32_t beg_;
        std::int32_t end_;
    };

    QueryRange Query(const BaiIndex& index, std::int32_t refId, std::int32_t beg, std::int32_t end);

    /// Range over records matching a ZMW whitelist (seek-per-record via index).
    /// Owns a dedicated sync reader (pipeline mode would restart per seek).
    class WhitelistRange : public std::ranges::view_interface<WhitelistRange>
    {
    public:
        class Iterator
        {
        public:
            using difference_type = std::ptrdiff_t;
            using value_type = RawRecord;
            using iterator_concept = std::input_iterator_tag;

            Iterator();
            Iterator(BamRawReader* reader, std::vector<std::int64_t> offsets);

            const RawRecord& operator*() const;
            const RawRecord* operator->() const;
            Iterator& operator++();
            void operator++(int);
            bool operator==(const Iterator& other) const;

        private:
            BamRawReader* reader_{nullptr};
            std::vector<std::int64_t> offsets_;
            std::size_t idx_{0};
            std::optional<RawRecord> current_;

            void Advance();
        };

        WhitelistRange(std::unique_ptr<BamRawReader> reader, std::vector<std::int64_t> offsets);
        Iterator begin();
        Iterator end();

    private:
        std::unique_ptr<BamRawReader> ownedReader_;
        std::vector<std::int64_t> offsets_;
    };

    WhitelistRange Whitelist(const ZmwWhitelist& whitelist);

    /// \brief Snapshot of BGZF pipeline metrics. Thread-safe, lock-free.
    /// In sync mode (BgzfWorkers == 0), only RecordsConsumed is populated.
    BgzfMetrics GetMetrics() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;

    bool RefillBuffer();
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_IO_BAMRAWREADER_HPP
