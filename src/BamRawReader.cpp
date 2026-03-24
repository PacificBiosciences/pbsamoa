#include <pbsamoa/io/BamRawReader.hpp>

#include <pbsamoa/core/Bgzf.hpp>
#include <pbsamoa/index/BaiIndex.hpp>
#include <pbsamoa/index/ZmwIndex.hpp>
#include <pbsamoa/index/ZmwWhitelist.hpp>

#include <algorithm>
#include <format>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {
namespace {

void AppendRecordToBatch(const RawRecord& record, std::vector<std::byte>& batchBuffer,
                         std::vector<RawRecordBatch::RecordExtent>& extents, std::size_t& batchSize)
{
    const auto raw{record.RawData()};
    const auto offset{static_cast<std::uint32_t>(std::size(batchBuffer))};
    batchBuffer.insert(std::ranges::end(batchBuffer), std::ranges::begin(raw),
                       std::ranges::end(raw));
    extents.emplace_back(offset, static_cast<std::uint32_t>(std::size(raw)));
    batchSize += std::size(raw);
}

template <typename ReadNext>
std::pair<std::optional<RawRecordBatch>, std::size_t> ReadRawBatch(ReadNext readNext,
                                                                   ByteLimit limit,
                                                                   std::size_t remaining)
{
    const std::size_t byteLimit{limit.Value()};
    std::vector<std::byte> batchBuffer;
    batchBuffer.reserve(byteLimit);
    std::vector<RawRecordBatch::RecordExtent> extents;
    std::size_t batchSize{0};

    while (std::size(extents) < remaining) {
        if ((!std::empty(extents)) && (batchSize >= byteLimit)) {
            break;
        }

        const auto record{readNext()};
        if (!record) {
            break;
        }
        AppendRecordToBatch(*record, batchBuffer, extents, batchSize);
    }

    if (std::empty(extents)) {
        return {std::nullopt, 0};
    }

    return {RawRecordBatch{std::move(batchBuffer), std::move(extents)}, std::size(extents)};
}

std::optional<std::size_t> ToRecordLimit(std::size_t recordLimit)
{
    if (recordLimit == 0) {
        return std::nullopt;
    }
    return recordLimit;
}

void AccumulatePoolMetrics(PoolMetrics& target, const PoolMetrics& source)
{
    target.QueueDepth += source.QueueDepth;
    target.PeakQueueDepth = std::max(target.PeakQueueDepth, source.PeakQueueDepth);
    target.ActiveWorkers += source.ActiveWorkers;
    target.PeakActiveWorkers = std::max(target.PeakActiveWorkers, source.PeakActiveWorkers);
    target.ResultQueueDepth += source.ResultQueueDepth;
    target.PeakResultQueueDepth =
        std::max(target.PeakResultQueueDepth, source.PeakResultQueueDepth);
}

void AccumulateBgzfMetrics(BgzfMetrics& target, const BgzfMetrics& source)
{
    target.BytesRead += source.BytesRead;
    target.BlocksRead += source.BlocksRead;
    target.BytesDecompressed += source.BytesDecompressed;
    AccumulatePoolMetrics(target.Pool, source.Pool);
    target.RecordsProduced += source.RecordsProduced;
    target.RecordsConsumed += source.RecordsConsumed;
    target.IoStalls += source.IoStalls;
    target.ConsumerStalls += source.ConsumerStalls;
    target.ReaderStalls += source.ReaderStalls;
    target.IoReadNs += source.IoReadNs;
    target.DecompressNs += source.DecompressNs;
    target.RecordParseNs += source.RecordParseNs;
}

bool RecordOverlapsQuery(const RawRecord& view, std::int32_t refId, std::int32_t beg,
                         std::int32_t end)
{
    return (view.RefId() == refId) && (view.Pos() < end) &&
           ((view.Pos() + view.ReferenceLength()) > beg);
}

bool RecordIsPastQuery(const RawRecord& view, std::int32_t refId, std::int32_t end)
{
    return (view.RefId() > refId) || ((view.RefId() == refId) && (view.Pos() >= end));
}

}  // namespace

struct BamRawReader::Impl
{
    struct CollectionState
    {
        BamCollection collection;
        BamRawReaderConfig config;
        std::size_t nextFileIndex{0};
        std::unique_ptr<BamRawReader> currentReader_;
        std::optional<std::size_t> recordLimit;
        std::size_t recordsRead{0};
        BgzfMetrics completedMetrics{};

        explicit CollectionState(BamCollection c, BamRawReaderConfig cfg)
            : collection{std::move(c)}, config{cfg}, recordLimit{ToRecordLimit(cfg.RecordLimit)}
        {
            const bool hasChunking{(cfg.ChunkNum > 0) && (cfg.TotalChunks > 0)};
            if (hasChunking || cfg.Whitelist) {
                throw std::invalid_argument{
                    "BamRawReaderConfig: chunking and whitelist are single-file only"};
            }
        }

        ~CollectionState() { FinishCurrentReader(); }

        bool AtRecordLimit() const { return recordLimit && (recordsRead >= *recordLimit); }

        std::size_t RemainingRecordBudget() const
        {
            if (!recordLimit) {
                return std::numeric_limits<std::size_t>::max();
            }
            return *recordLimit - recordsRead;
        }

        void FinishCurrentReader()
        {
            if (!currentReader_) {
                return;
            }
            AccumulateBgzfMetrics(completedMetrics, currentReader_->GetMetrics());
            currentReader_.reset();
        }

        bool OpenNextReader()
        {
            FinishCurrentReader();
            if (nextFileIndex >= collection.Size()) {
                return false;
            }

            const BamFile& file{collection.Files()[nextFileIndex]};
            currentReader_ = std::make_unique<BamRawReader>(
                file.Filename(), BamRawReaderConfig{.BgzfWorkers = config.BgzfWorkers});
            ++nextFileIndex;
            return true;
        }

        BgzfMetrics Metrics() const
        {
            BgzfMetrics metrics{completedMetrics};
            if (currentReader_) {
                AccumulateBgzfMetrics(metrics, currentReader_->GetMetrics());
            }
            return metrics;
        }
    };

    std::filesystem::path path;
    std::unique_ptr<BgzfReader> bgzf;

    std::optional<std::size_t> recordLimit;
    std::size_t recordsRead{0};
    std::int32_t numZmws_{-1};
    std::unique_ptr<CollectionState> collection;

    bool IsCollection() const { return static_cast<bool>(collection); }

    bool AtRecordLimit() const
    {
        if (collection) {
            return collection->AtRecordLimit();
        }
        return recordLimit && (recordsRead >= *recordLimit);
    }

    std::size_t RemainingRecordBudget() const
    {
        if (collection) {
            return collection->RemainingRecordBudget();
        }
        if (!recordLimit) {
            return std::numeric_limits<std::size_t>::max();
        }
        return *recordLimit - recordsRead;
    }

    explicit Impl(const std::filesystem::path& p, BamRawReaderConfig config)
        : path{p}
        , bgzf{std::make_unique<BgzfReader>(p, config.BgzfWorkers)}
        , recordLimit{ToRecordLimit(config.RecordLimit)}
    {
        const bool hasChunking{(config.ChunkNum > 0) && (config.TotalChunks > 0)};
        if (hasChunking && config.Whitelist) {
            throw std::invalid_argument{
                "BamRawReaderConfig: Whitelist and ChunkNum/TotalChunks are mutually "
                "exclusive"};
        }

        if ((config.ChunkNum > 0) && (config.TotalChunks > 0)) {
            ApplyChunkConfig(config.ChunkNum, config.TotalChunks);
        }
    }

    explicit Impl(BamCollection c, BamRawReaderConfig config)
        : collection{std::make_unique<CollectionState>(std::move(c), config)}
    {
    }

    void ApplyChunkConfig(std::int32_t chunkNum, std::int32_t totalChunks)
    {
        if (totalChunks < 1) {
            throw std::invalid_argument{
                std::format("TotalChunks must be >= 1, got {}", totalChunks)};
        }
        if ((chunkNum < 1) || (chunkNum > totalChunks)) {
            throw std::invalid_argument{
                std::format("ChunkNum must be in [1, {}], got {}", totalChunks, chunkNum)};
        }

        const ZmwIndex index{ZmwIndex::Open(path)};
        const std::vector<ZmwIdentity> unique{index.UniqueZmws()};
        const std::ptrdiff_t numZmws{std::ssize(unique)};

        if (numZmws == 0) {
            throw std::runtime_error{std::format("ZMW index is empty for: {}", path.string())};
        }

        const std::ptrdiff_t startIdx{(numZmws * static_cast<std::ptrdiff_t>(chunkNum - 1)) /
                                      totalChunks};
        const std::ptrdiff_t endIdx{(numZmws * static_cast<std::ptrdiff_t>(chunkNum)) /
                                    totalChunks};

        recordsRead = 0;

        if (startIdx >= endIdx) {
            recordLimit = std::size_t{0};
            return;
        }

        numZmws_ = static_cast<std::int32_t>(endIdx - startIdx);

        const std::span<const ZmwIdentity> chunkZmws{std::data(unique) + startIdx,
                                                     static_cast<std::size_t>(endIdx - startIdx)};
        const std::vector<std::int64_t> chunkOffsets{index.Find(chunkZmws)};

        recordLimit = std::optional{std::size(chunkOffsets)};

        const VirtualOffset startOffset(index.FirstOffset(unique[startIdx]));
        bgzf->Seek(startOffset);
    }

    ~Impl() = default;
    Impl(const Impl&) = delete;
    Impl& operator=(const Impl&) = delete;
};

BamRawReader::BamRawReader(const std::filesystem::path& path, BamRawReaderConfig config)
    : impl_{std::make_unique<Impl>(path, config)}
{
}

BamRawReader::BamRawReader(BamCollection collection, BamRawReaderConfig config)
    : impl_{std::make_unique<Impl>(std::move(collection), config)}
{
}

BamRawReader::BamRawReader(std::vector<std::filesystem::path> paths, BamRawReaderConfig config)
    : BamRawReader{BamCollection{std::move(paths)}, config}
{
}

BamRawReader::BamRawReader(std::vector<BamFile> files, BamRawReaderConfig config)
    : BamRawReader{BamCollection{std::move(files)}, config}
{
}

BamRawReader::~BamRawReader() = default;
BamRawReader::BamRawReader(BamRawReader&&) noexcept = default;
BamRawReader& BamRawReader::operator=(BamRawReader&&) noexcept = default;

const SamHeader& BamRawReader::Header() const
{
    if (impl_->collection) {
        return impl_->collection->collection.Header();
    }
    return impl_->bgzf->Header();
}

std::optional<RawRecord> BamRawReader::ReadRecord()
{
    if (impl_->AtRecordLimit()) {
        return std::nullopt;
    }

    if (!impl_->collection) {
        auto result{impl_->bgzf->ReadRecord()};
        if (result) {
            ++impl_->recordsRead;
        }
        return result;
    }

    auto& collection{*impl_->collection};
    while (true) {
        if (!collection.currentReader_ && !collection.OpenNextReader()) {
            return std::nullopt;
        }

        std::optional<RawRecord> result{collection.currentReader_->ReadRecord()};
        if (result) {
            ++collection.recordsRead;
            return result;
        }

        if (!collection.OpenNextReader()) {
            return std::nullopt;
        }
    }
}

std::optional<RawRecordBatch> BamRawReader::ReadBatch(ByteLimit limit)
{
    if (impl_->AtRecordLimit()) {
        return std::nullopt;
    }

    const std::size_t remaining{impl_->RemainingRecordBudget()};
    if (!impl_->collection) {
        auto [batch, recordCount] =
            ReadRawBatch([this]() { return impl_->bgzf->ReadRecord(); }, limit, remaining);
        impl_->recordsRead += recordCount;
        return batch;
    }

    return ReadRawBatch([this]() { return ReadRecord(); }, limit, remaining).first;
}

void BamRawReader::Seek(VirtualOffset offset)
{
    if (impl_->collection) {
        throw std::logic_error{"BamRawReader::Seek is not supported for multi-BAM input"};
    }
    impl_->bgzf->Seek(offset);
}

VirtualOffset BamRawReader::Tell() const
{
    if (impl_->collection) {
        throw std::logic_error{"BamRawReader::Tell is not supported for multi-BAM input"};
    }
    return impl_->bgzf->Tell();
}

// --- RecordRange ---

BamRawReader::RecordRange::RecordRange(BamRawReader* reader) : reader_{reader} {}

BamRawReader::RecordRange::Iterator BamRawReader::RecordRange::begin() { return Iterator{reader_}; }

BamRawReader::RecordRange::Iterator BamRawReader::RecordRange::end() { return Iterator{}; }

// --- Iterator ---

BamRawReader::RecordRange::Iterator::Iterator() = default;

BamRawReader::RecordRange::Iterator::Iterator(BamRawReader* reader) : reader_{reader}
{
    current_ = reader_->ReadRecord();
    if (!current_) {
        reader_ = nullptr;
    }
}

const RawRecord& BamRawReader::RecordRange::Iterator::operator*() const { return *current_; }

const RawRecord* BamRawReader::RecordRange::Iterator::operator->() const { return &*current_; }

BamRawReader::RecordRange::Iterator& BamRawReader::RecordRange::Iterator::operator++()
{
    current_ = reader_->ReadRecord();
    if (!current_) {
        reader_ = nullptr;
    }
    return *this;
}

void BamRawReader::RecordRange::Iterator::operator++(int) { ++(*this); }

bool BamRawReader::RecordRange::Iterator::operator==(const Iterator& other) const
{
    return reader_ == other.reader_;
}

BamRawReader::RecordRange BamRawReader::Records() { return RecordRange{this}; }

// --- QueryRange ---

BamRawReader::QueryRange::QueryRange(BamRawReader* reader, std::vector<Chunk> chunks,
                                     std::int32_t refId, std::int32_t beg, std::int32_t end)
    : reader_{reader}, chunks_{std::move(chunks)}, refId_{refId}, beg_{beg}, end_{end}
{
}

BamRawReader::QueryRange::Iterator BamRawReader::QueryRange::begin()
{
    return Iterator{reader_, chunks_, refId_, beg_, end_};
}

BamRawReader::QueryRange::Iterator BamRawReader::QueryRange::end() { return Iterator{}; }

// --- QueryRange::Iterator ---

BamRawReader::QueryRange::Iterator::Iterator() = default;

BamRawReader::QueryRange::Iterator::Iterator(BamRawReader* reader, std::vector<Chunk> chunks,
                                             std::int32_t refId, std::int32_t beg, std::int32_t end)
    : reader_{reader}, chunks_{std::move(chunks)}, chunkIdx_{0}, refId_{refId}, beg_{beg}, end_{end}
{
    if (std::empty(chunks_)) {
        reader_ = nullptr;
        return;
    }
    reader_->Seek(chunks_[0].Begin);
    Advance();
}

const RawRecord& BamRawReader::QueryRange::Iterator::operator*() const { return *current_; }

const RawRecord* BamRawReader::QueryRange::Iterator::operator->() const { return &*current_; }

BamRawReader::QueryRange::Iterator& BamRawReader::QueryRange::Iterator::operator++()
{
    Advance();
    return *this;
}

void BamRawReader::QueryRange::Iterator::operator++(int) { ++(*this); }

bool BamRawReader::QueryRange::Iterator::operator==(const Iterator& other) const
{
    return reader_ == other.reader_;
}

void BamRawReader::QueryRange::Iterator::Advance()
{
    const auto advanceChunk = [this]() -> bool {
        ++chunkIdx_;
        if (chunkIdx_ >= std::size(chunks_)) {
            return false;
        }
        reader_->Seek(chunks_[chunkIdx_].Begin);
        return true;
    };

    while (chunkIdx_ < std::size(chunks_)) {
        auto view{reader_->ReadRecord()};
        if (!view) {
            if (!advanceChunk()) {
                break;
            }
            continue;
        }

        const VirtualOffset currentPos{reader_->Tell()};
        if (currentPos > chunks_[chunkIdx_].End) {
            if (!advanceChunk()) {
                break;
            }
            continue;
        }

        if (RecordOverlapsQuery(*view, refId_, beg_, end_)) {
            current_ = std::move(view);
            return;
        }

        if (RecordIsPastQuery(*view, refId_, end_)) {
            break;
        }
    }

    reader_ = nullptr;
    current_.reset();
}

BamRawReader::QueryRange BamRawReader::Query(const BaiIndex& index, std::int32_t refId,
                                             std::int32_t beg, std::int32_t end)
{
    if (impl_->collection) {
        throw std::logic_error{
            "BamRawReader::Query(index, ...) is single-file "
            "only; use a collection-aware reader"};
    }
    std::vector<Chunk> chunks{index.Query(refId, beg, end)};
    return QueryRange{this, std::move(chunks), refId, beg, end};
}

// --- WhitelistRange ---

BamRawReader::WhitelistRange::Iterator::Iterator() = default;

BamRawReader::WhitelistRange::Iterator::Iterator(BamRawReader* reader,
                                                 std::vector<std::int64_t> offsets)
    : reader_{reader}, offsets_{std::move(offsets)}
{
    if (std::empty(offsets_)) {
        reader_ = nullptr;
        return;
    }
    Advance();
}

const RawRecord& BamRawReader::WhitelistRange::Iterator::operator*() const { return *current_; }

const RawRecord* BamRawReader::WhitelistRange::Iterator::operator->() const { return &*current_; }

BamRawReader::WhitelistRange::Iterator& BamRawReader::WhitelistRange::Iterator::operator++()
{
    ++idx_;
    Advance();
    return *this;
}

void BamRawReader::WhitelistRange::Iterator::operator++(int) { ++(*this); }

bool BamRawReader::WhitelistRange::Iterator::operator==(const Iterator& other) const
{
    return reader_ == other.reader_;
}

void BamRawReader::WhitelistRange::Iterator::Advance()
{
    if (idx_ >= std::size(offsets_)) {
        reader_ = nullptr;
        current_.reset();
        return;
    }

    reader_->Seek(VirtualOffset(offsets_[idx_]));
    auto record{reader_->ReadRecord()};
    if (!record) {
        reader_ = nullptr;
        current_.reset();
        return;
    }
    current_ = std::move(record);
}

BamRawReader::WhitelistRange::WhitelistRange(std::unique_ptr<BamRawReader> reader,
                                             std::vector<std::int64_t> offsets)
    : ownedReader_{std::move(reader)}, offsets_{std::move(offsets)}
{
}

BamRawReader::WhitelistRange::Iterator BamRawReader::WhitelistRange::begin()
{
    return Iterator{ownedReader_.get(), offsets_};
}

BamRawReader::WhitelistRange::Iterator BamRawReader::WhitelistRange::end() { return Iterator{}; }

BamRawReader::WhitelistRange BamRawReader::Whitelist(const ZmwWhitelist& whitelist)
{
    if (impl_->collection) {
        throw std::logic_error{"BamRawReader::Whitelist is single-file only"};
    }

    const ZmwIndex index{ZmwIndex::Open(impl_->path)};
    std::vector<std::int64_t> offsets{whitelist.Resolve(index)};
    auto syncReader{std::make_unique<BamRawReader>(impl_->path)};
    return WhitelistRange{std::move(syncReader), std::move(offsets)};
}

std::int32_t BamRawReader::NumZmws() const
{
    if (impl_->collection) {
        return -1;
    }
    if (impl_->numZmws_ >= 0) {
        return impl_->numZmws_;
    }
    // Fallback: open the ZMW index to compute the total count.
    // This only runs for non-chunked single-file reads where
    // ApplyChunkConfig did not already populate the count.
    try {
        impl_->numZmws_ = static_cast<std::int32_t>(ZmwIndex::Open(impl_->path).NumZmws());
    } catch (...) {
    }
    return impl_->numZmws_;
}

// --- GetMetrics ---

BgzfMetrics BamRawReader::GetMetrics() const
{
    if (impl_->collection) {
        return impl_->collection->Metrics();
    }
    return impl_->bgzf->GetMetrics();
}

// --- ThrowPolicy / SkipPolicy ---

void ThrowPolicy::OnCorruptRecord(std::string_view message) const
{
    throw std::runtime_error{std::string{message}};
}

void SkipPolicy::OnCorruptRecord(std::string_view /*message*/) { ++skipped_; }

std::size_t SkipPolicy::SkippedCount() const { return skipped_; }

}  // namespace Samoa
}  // namespace PacBio
