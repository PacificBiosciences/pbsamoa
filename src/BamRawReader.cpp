#include <pbsamoa/io/BamRawReader.hpp>

#include <pbsamoa/core/Bgzf.hpp>
#include <pbsamoa/index/BaiIndex.hpp>
#include <pbsamoa/index/ZmwIndex.hpp>
#include <pbsamoa/index/ZmwWhitelist.hpp>

#include <algorithm>
#include <format>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

struct BamRawReader::Impl
{
    std::filesystem::path path;
    BgzfReader bgzf;

    // Chunking support
    bool hasRecordLimit{false};
    std::size_t recordLimit{0};
    std::size_t recordsRead{0};

    explicit Impl(const std::filesystem::path& p, BamRawReaderConfig config)
        : path{p}
        , bgzf{p, config.BgzfWorkers}
        , hasRecordLimit{config.RecordLimit > 0}
        , recordLimit{config.RecordLimit}
    {
        // Validate mutual exclusion
        const bool hasChunking{(config.ChunkNum > 0) && (config.TotalChunks > 0)};
        if (hasChunking && config.Whitelist.has_value()) {
            throw std::invalid_argument{
                "BamRawReaderConfig: Whitelist and ChunkNum/TotalChunks are mutually exclusive"};
        }

        // Apply chunking if requested
        if ((config.ChunkNum > 0) && (config.TotalChunks > 0)) {
            ApplyChunkConfig(config.ChunkNum, config.TotalChunks);
        }
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

        hasRecordLimit = true;
        recordsRead = 0;

        if (startIdx >= endIdx) {
            // This chunk has no assigned ZMWs.
            recordLimit = 0;
            return;
        }

        const std::span<const ZmwIdentity> chunkZmws{std::data(unique) + startIdx,
                                                     static_cast<std::size_t>(endIdx - startIdx)};
        const std::vector<std::int64_t> chunkOffsets{index.Find(chunkZmws)};

        recordLimit = std::size(chunkOffsets);

        // Seek to the first record in this chunk
        const VirtualOffset startOffset(index.FirstOffset(unique[startIdx]));
        bgzf.Seek(startOffset);
    }

    ~Impl() = default;
    Impl(const Impl&) = delete;
    Impl& operator=(const Impl&) = delete;
};

BamRawReader::BamRawReader(const std::filesystem::path& path, BamRawReaderConfig config)
    : impl_{std::make_unique<Impl>(path, config)}
{
}

BamRawReader::~BamRawReader() = default;
BamRawReader::BamRawReader(BamRawReader&&) noexcept = default;
BamRawReader& BamRawReader::operator=(BamRawReader&&) noexcept = default;

const SamHeader& BamRawReader::Header() const { return impl_->bgzf.Header(); }

std::optional<RawRecord> BamRawReader::ReadRecord()
{
    if (impl_->hasRecordLimit && (impl_->recordsRead >= impl_->recordLimit)) {
        return std::nullopt;
    }

    auto result{impl_->bgzf.ReadRecord()};
    if (result.has_value()) {
        ++impl_->recordsRead;
    }
    return result;
}

std::optional<RawRecordBatch> BamRawReader::ReadBatch(ByteLimit limit)
{
    if (impl_->hasRecordLimit && (impl_->recordsRead >= impl_->recordLimit)) {
        return std::nullopt;
    }

    const std::size_t remaining{impl_->hasRecordLimit ? (impl_->recordLimit - impl_->recordsRead)
                                                      : std::numeric_limits<std::size_t>::max()};

    std::vector<std::byte> batchBuffer;
    batchBuffer.reserve(limit.Value());
    std::vector<RawRecordBatch::RecordExtent> extents;
    std::size_t batchSize{0};

    while (std::size(extents) < remaining) {
        if ((!std::empty(extents)) && (batchSize >= limit.Value())) {
            break;
        }

        const auto rec{impl_->bgzf.ReadRecord()};
        if (!rec.has_value()) {
            break;
        }

        const auto raw{rec->RawData()};
        const std::uint32_t offset = std::size(batchBuffer);
        batchBuffer.insert(std::ranges::end(batchBuffer), std::ranges::begin(raw),
                           std::ranges::end(raw));
        extents.push_back(
            RawRecordBatch::RecordExtent{offset, static_cast<std::uint32_t>(std::size(raw))});
        batchSize += std::size(raw);
    }

    if (std::empty(extents)) {
        return std::nullopt;
    }

    impl_->recordsRead += std::size(extents);
    return RawRecordBatch{std::move(batchBuffer), std::move(extents)};
}

void BamRawReader::Seek(VirtualOffset offset) { impl_->bgzf.Seek(offset); }

VirtualOffset BamRawReader::Tell() const { return impl_->bgzf.Tell(); }

// --- RecordRange ---

BamRawReader::RecordRange::RecordRange(BamRawReader* reader) : reader_{reader} {}

BamRawReader::RecordRange::Iterator BamRawReader::RecordRange::begin() { return Iterator{reader_}; }

BamRawReader::RecordRange::Iterator BamRawReader::RecordRange::end() { return Iterator{}; }

// --- Iterator ---

BamRawReader::RecordRange::Iterator::Iterator() = default;

BamRawReader::RecordRange::Iterator::Iterator(BamRawReader* reader) : reader_{reader}
{
    current_ = reader_->ReadRecord();
    if (!current_.has_value()) {
        reader_ = nullptr;
    }
}

const RawRecord& BamRawReader::RecordRange::Iterator::operator*() const { return *current_; }

const RawRecord* BamRawReader::RecordRange::Iterator::operator->() const { return &*current_; }

BamRawReader::RecordRange::Iterator& BamRawReader::RecordRange::Iterator::operator++()
{
    current_ = reader_->ReadRecord();
    if (!current_.has_value()) {
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
    // Seek to start of first chunk
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
    while (chunkIdx_ < std::size(chunks_)) {
        // Try to read the next record
        auto view{reader_->ReadRecord()};
        if (!view.has_value()) {
            // EOF or error - try next chunk
            ++chunkIdx_;
            if (chunkIdx_ < std::size(chunks_)) {
                reader_->Seek(chunks_[chunkIdx_].Begin);
            }
            continue;
        }

        // Check if we've passed the end of the current chunk
        const VirtualOffset currentPos{reader_->Tell()};
        if (currentPos > chunks_[chunkIdx_].End) {
            // Move to next chunk
            ++chunkIdx_;
            if (chunkIdx_ < std::size(chunks_)) {
                reader_->Seek(chunks_[chunkIdx_].Begin);
            }
            continue;
        }

        // Filter: must be on the correct reference and overlap [beg, end)
        if ((view->RefId() == refId_) && (view->Pos() < end_) &&
            (view->Pos() + view->ReferenceLength() > beg_)) {
            current_ = std::move(view);
            return;
        }

        // If we've passed the query region on this reference, we can stop
        if ((view->RefId() > refId_) || ((view->RefId() == refId_) && (view->Pos() >= end_))) {
            break;
        }
    }

    // No more matching records
    reader_ = nullptr;
    current_.reset();
}

BamRawReader::QueryRange BamRawReader::Query(const BaiIndex& index, std::int32_t refId,
                                             std::int32_t beg, std::int32_t end)
{
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
    auto rec{reader_->ReadRecord()};
    if (!rec.has_value()) {
        reader_ = nullptr;
        current_.reset();
        return;
    }
    current_ = std::move(rec);
}

BamRawReader::WhitelistRange::WhitelistRange(std::unique_ptr<BamRawReader> reader,
                                             std::vector<std::int64_t> offsets)
    : ownedReader_{std::move(reader)}, offsets_{std::move(offsets)}
{
}

BamRawReader::WhitelistRange::Iterator BamRawReader::WhitelistRange::begin()
{
    return Iterator{ownedReader_.get(), std::move(offsets_)};
}

BamRawReader::WhitelistRange::Iterator BamRawReader::WhitelistRange::end() { return Iterator{}; }

BamRawReader::WhitelistRange BamRawReader::Whitelist(const ZmwWhitelist& whitelist)
{
    const ZmwIndex index{ZmwIndex::Open(impl_->path)};
    std::vector<std::int64_t> offsets{whitelist.Resolve(index)};
    // Use a dedicated sync reader for seek-per-record iteration.
    // Pipeline mode would restart per seek, which is both slow and
    // incorrect for within-block offsets in the header block.
    auto syncReader{std::make_unique<BamRawReader>(impl_->path)};
    return WhitelistRange{std::move(syncReader), std::move(offsets)};
}

// --- GetMetrics ---

BgzfMetrics BamRawReader::GetMetrics() const { return impl_->bgzf.GetMetrics(); }

// --- ThrowPolicy / SkipPolicy ---

void ThrowPolicy::OnCorruptRecord(std::string_view message) const
{
    throw std::runtime_error{std::string{message}};
}

void SkipPolicy::OnCorruptRecord(std::string_view /*message*/) { ++skipped_; }

std::size_t SkipPolicy::SkippedCount() const { return skipped_; }

}  // namespace Samoa
}  // namespace PacBio
