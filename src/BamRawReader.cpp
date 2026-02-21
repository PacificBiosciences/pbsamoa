#include <pbsamoa/io/BamRawReader.hpp>

#include "BinaryUtils.hpp"

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

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>

namespace PacBio {
namespace Samoa {

namespace {

constexpr std::size_t RECORD_BLOCK_SIZE_FIELD{4};
constexpr std::size_t INITIAL_HEADER_BUF_SIZE{MAX_DECOMPRESSED_BLOCK_SIZE * 4};

}  // namespace

struct BamRawReader::Impl
{
    std::filesystem::path path;

    // Sync path (nullopt in pipeline mode)
    std::optional<BgzfReader> bgzf;
    SamHeader header{};

    // Record buffer: holds decompressed BGZF data for sync record iteration
    std::vector<std::byte> buffer{};
    std::size_t bufferSize{0};  // valid bytes in buffer
    std::size_t bufferPos{0};   // current read position in buffer
    bool eof{false};

    // Pipeline path (null in sync mode)
    std::unique_ptr<BgzfPipeline> pipeline;

    // Chunking support
    std::size_t recordLimit{0};  // 0 = unlimited
    std::size_t recordsRead{0};

    explicit Impl(const std::filesystem::path& p, BamRawReaderConfig config)
        : path{p}, recordLimit{config.RecordLimit}
    {
        if (config.BgzfWorkers > 0) {
            pipeline = std::make_unique<BgzfPipeline>(path, config.BgzfWorkers);
            pipeline->ParseHeader();
        } else {
            bgzf.emplace(path);
            ParseHeader();
        }

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

        const double chunkSize{1.0 * numZmws / totalChunks};
        const std::ptrdiff_t startIdx{(chunkNum == 1) ? 0
                                                      : std::lround(chunkSize * (chunkNum - 1))};
        const std::ptrdiff_t endIdx{(chunkNum == totalChunks) ? numZmws
                                                              : std::lround(chunkSize * chunkNum)};

        const std::vector<ZmwIdentity> chunkZmws{std::begin(unique) + startIdx,
                                                 std::begin(unique) + endIdx};
        const std::vector<std::int64_t> chunkOffsets{index.Find(chunkZmws)};

        recordLimit = std::size(chunkOffsets);
        recordsRead = 0;

        // Seek to the first record in this chunk
        const VirtualOffset startOffset(index.FirstOffset(unique[startIdx]));
        if (pipeline) {
            pipeline->Seek(startOffset);
        } else {
            bgzf->Seek(startOffset);
            bufferPos = 0;
            bufferSize = 0;
            eof = false;

            const std::uint16_t withinBlock{startOffset.WithinBlockOffset()};
            if (withinBlock > 0) {
                // Need to decompress the block and skip ahead
                const std::size_t needed{bufferSize + MAX_DECOMPRESSED_BLOCK_SIZE};
                if (std::size(buffer) < needed) {
                    buffer.resize(needed);
                }
                const std::optional<std::size_t> bytesRead{bgzf->ReadBlock(
                    std::span<std::byte>{buffer}.subspan(bufferSize, MAX_DECOMPRESSED_BLOCK_SIZE))};
                if (bytesRead.has_value() && (*bytesRead > 0)) {
                    bufferSize += *bytesRead;
                    if (bufferSize >= withinBlock) {
                        bufferPos = withinBlock;
                    }
                }
            }
        }
    }

    ~Impl() = default;
    Impl(const Impl&) = delete;
    Impl& operator=(const Impl&) = delete;

    void ParseHeader()
    {
        // Read BGZF blocks until we have the complete BAM header
        std::vector<std::byte> headerBuf(INITIAL_HEADER_BUF_SIZE);
        std::size_t headerLen{0};
        std::vector<std::byte> blockBuf(MAX_DECOMPRESSED_BLOCK_SIZE);

        while (true) {
            const std::optional<std::size_t> bytesRead{
                bgzf->ReadBlock(std::span<std::byte>{blockBuf})};

            if (!bytesRead.has_value()) {
                throw std::runtime_error{"Failed to read BGZF block while parsing BAM header"};
            }
            if (*bytesRead == 0) {
                break;  // EOF
            }

            // Grow header buffer if needed
            if (headerLen + *bytesRead > std::size(headerBuf)) {
                headerBuf.resize(headerLen + *bytesRead + MAX_DECOMPRESSED_BLOCK_SIZE);
            }

            std::ranges::copy_n(std::data(blockBuf), *bytesRead, std::data(headerBuf) + headerLen);
            headerLen += *bytesRead;

            // Check if we have enough data for the complete header
            const std::size_t hdrSize{ComputeHeaderSize(std::data(headerBuf), headerLen)};
            if (hdrSize > 0) {
                // Parse the header
                header = SamHeader::FromBamHeaderBlock(
                    std::span<const std::byte>{std::data(headerBuf), hdrSize});

                // Copy remaining bytes (after header) into the record buffer
                const std::size_t remaining{headerLen - hdrSize};
                if (remaining > 0) {
                    buffer.resize(std::ranges::max(remaining, MAX_DECOMPRESSED_BLOCK_SIZE));
                    std::ranges::copy_n(std::data(headerBuf) + hdrSize, remaining,
                                        std::data(buffer));
                    bufferSize = remaining;
                } else {
                    buffer.resize(MAX_DECOMPRESSED_BLOCK_SIZE);
                }
                return;
            }
        }

        // If we got here, we hit EOF before completing the header
        if (headerLen == 0) {
            throw std::runtime_error{"Empty BAM file: no BGZF blocks"};
        }

        // Try to parse what we have (may be a header-only file)
        const std::size_t hdrSize{ComputeHeaderSize(std::data(headerBuf), headerLen)};
        if (hdrSize > 0) {
            header = SamHeader::FromBamHeaderBlock(
                std::span<const std::byte>{std::data(headerBuf), hdrSize});
            buffer.resize(MAX_DECOMPRESSED_BLOCK_SIZE);
            eof = true;
        } else {
            throw std::runtime_error{"Incomplete BAM header"};
        }
    }
};

BamRawReader::BamRawReader(const std::filesystem::path& path, BamRawReaderConfig config)
    : impl_{std::make_unique<Impl>(path, config)}
{
}

BamRawReader::~BamRawReader() = default;
BamRawReader::BamRawReader(BamRawReader&&) noexcept = default;
BamRawReader& BamRawReader::operator=(BamRawReader&&) noexcept = default;

const SamHeader& BamRawReader::Header() const
{
    if (impl_->pipeline) {
        return impl_->pipeline->Header();
    }
    return impl_->header;
}

bool BamRawReader::RefillBuffer()
{
    // Compact unconsumed data to buffer front only when >50% consumed.
    // This halves memmove frequency vs. compacting on every call.
    const std::size_t remaining{impl_->bufferSize - impl_->bufferPos};
    if ((impl_->bufferPos > 0) &&
        ((remaining == 0) || (impl_->bufferPos > impl_->bufferSize / 2))) {
        if (remaining > 0) {
            std::memmove(std::data(impl_->buffer), std::data(impl_->buffer) + impl_->bufferPos,
                         remaining);
        }
        impl_->bufferPos = 0;
        impl_->bufferSize = remaining;
    }

    if (impl_->eof) {
        return remaining > 0;
    }

    // Ensure buffer can hold current data + one full decompressed block
    const std::size_t needed{impl_->bufferSize + MAX_DECOMPRESSED_BLOCK_SIZE};
    if (std::size(impl_->buffer) < needed) {
        impl_->buffer.resize(needed);
    }

    const std::optional<std::size_t> bytesRead{
        impl_->bgzf->ReadBlock(std::span<std::byte>{impl_->buffer}.subspan(
            impl_->bufferSize, MAX_DECOMPRESSED_BLOCK_SIZE))};

    if (!bytesRead.has_value()) {
        // Read error - treat as EOF
        impl_->eof = true;
        return remaining > 0;
    }
    if (*bytesRead == 0) {
        impl_->eof = true;
        return remaining > 0;
    }

    impl_->bufferSize += *bytesRead;
    return true;
}

std::optional<RawRecord> BamRawReader::ReadRecord()
{
    // Check record limit
    if ((impl_->recordLimit > 0) && (impl_->recordsRead >= impl_->recordLimit)) {
        return std::nullopt;
    }

    // Pipeline path: delegate directly
    if (impl_->pipeline) {
        auto result{impl_->pipeline->ReadRecord()};
        if (result.has_value()) {
            ++impl_->recordsRead;
        }
        return result;
    }

    // Sync path: parse from internal buffer, wrap in owning RawRecord
    // Ensure we have at least 4 bytes for block_size
    while ((impl_->bufferSize - impl_->bufferPos) < RECORD_BLOCK_SIZE_FIELD) {
        if (!RefillBuffer()) {
            return std::nullopt;  // EOF
        }
    }

    const std::byte* pos{std::data(impl_->buffer) + impl_->bufferPos};
    const std::uint32_t blockSize{ReadU32LE(pos)};

    if (blockSize == 0) {
        return std::nullopt;
    }

    const std::size_t totalRecordBytes{RECORD_BLOCK_SIZE_FIELD + blockSize};

    // Ensure we have all the record data
    while ((impl_->bufferSize - impl_->bufferPos) < totalRecordBytes) {
        if (!RefillBuffer()) {
            return std::nullopt;  // Truncated record at EOF
        }
    }

    // Create owning view from the record data (after block_size field)
    pos = std::data(impl_->buffer) + impl_->bufferPos + RECORD_BLOCK_SIZE_FIELD;
    impl_->bufferPos += totalRecordBytes;
    ++impl_->recordsRead;

    return RawRecord{std::span<const std::byte>{pos, blockSize}};
}

std::optional<RawRecordBatch> BamRawReader::ReadBatch(ByteLimit limit)
{
    // Check record limit
    if ((impl_->recordLimit > 0) && (impl_->recordsRead >= impl_->recordLimit)) {
        return std::nullopt;
    }

    const std::size_t remaining{(impl_->recordLimit > 0) ? (impl_->recordLimit - impl_->recordsRead)
                                                         : std::numeric_limits<std::size_t>::max()};

    // Pipeline path: accumulate records from ReadRecord()
    if (impl_->pipeline) {
        std::vector<std::byte> batchBuffer;
        batchBuffer.reserve(limit.Value());
        std::vector<RawRecordBatch::RecordExtent> extents;
        std::size_t batchSize{0};

        while (std::size(extents) < remaining) {
            if ((!std::empty(extents)) && (batchSize >= limit.Value())) {
                break;
            }

            const auto rec{impl_->pipeline->ReadRecord()};
            if (!rec.has_value()) {
                break;
            }

            const auto raw{rec->RawData()};
            const std::uint32_t offset{static_cast<std::uint32_t>(std::size(batchBuffer))};
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

    // Sync path: read directly from internal buffer (avoids extra copy)
    std::vector<std::byte> batchBuffer;
    batchBuffer.reserve(limit.Value());
    std::vector<RawRecordBatch::RecordExtent> extents;
    extents.reserve(std::max(std::size_t{1}, limit.Value() / 10000));
    std::size_t batchSize{0};

    while (std::size(extents) < remaining) {
        // Check if we've reached the budget (always include at least 1 record)
        if ((!std::empty(extents)) && (batchSize >= limit.Value())) {
            break;
        }

        // Ensure 4 bytes for block_size
        while ((impl_->bufferSize - impl_->bufferPos) < RECORD_BLOCK_SIZE_FIELD) {
            if (!RefillBuffer()) {
                // EOF - return what we have
                if (std::empty(extents)) {
                    return std::nullopt;
                }
                impl_->recordsRead += std::size(extents);
                return RawRecordBatch{std::move(batchBuffer), std::move(extents)};
            }
        }

        const std::byte* pos{std::data(impl_->buffer) + impl_->bufferPos};
        const std::uint32_t blockSize{ReadU32LE(pos)};

        if (blockSize == 0) {
            break;
        }

        const std::size_t totalRecordBytes{RECORD_BLOCK_SIZE_FIELD + blockSize};

        // Ensure we have all the record data
        while ((impl_->bufferSize - impl_->bufferPos) < totalRecordBytes) {
            if (!RefillBuffer()) {
                // Truncated record at EOF
                if (std::empty(extents)) {
                    return std::nullopt;
                }
                impl_->recordsRead += std::size(extents);
                return RawRecordBatch{std::move(batchBuffer), std::move(extents)};
            }
        }

        // Copy record data into batch buffer
        const std::byte* recordData{std::data(impl_->buffer) + impl_->bufferPos +
                                    RECORD_BLOCK_SIZE_FIELD};
        const std::uint32_t offset = std::size(batchBuffer);
        batchBuffer.insert(std::ranges::end(batchBuffer), recordData, recordData + blockSize);
        extents.push_back(RawRecordBatch::RecordExtent{offset, blockSize});

        impl_->bufferPos += totalRecordBytes;
        batchSize += blockSize;
    }

    if (std::empty(extents)) {
        return std::nullopt;
    }

    impl_->recordsRead += std::size(extents);
    return RawRecordBatch{std::move(batchBuffer), std::move(extents)};
}

void BamRawReader::Seek(VirtualOffset offset)
{
    if (impl_->pipeline) {
        impl_->pipeline->Seek(offset);
        impl_->bufferPos = 0;
        impl_->bufferSize = 0;
        impl_->eof = false;

        // Handle within-block offset for pipeline path too
        const std::uint16_t withinBlock{offset.WithinBlockOffset()};
        if (withinBlock > 0) {
            RefillBuffer();
            if (impl_->bufferSize >= withinBlock) {
                impl_->bufferPos = withinBlock;
            }
        }
        return;
    }

    impl_->bgzf->Seek(offset);
    impl_->bufferPos = 0;
    impl_->bufferSize = 0;
    impl_->eof = false;

    // Handle within-block offset: decompress the block and skip to the
    // requested position within it.
    const std::uint16_t withinBlock{offset.WithinBlockOffset()};
    if (withinBlock > 0) {
        RefillBuffer();
        if (impl_->bufferSize >= withinBlock) {
            impl_->bufferPos = withinBlock;
        }
    }
}

VirtualOffset BamRawReader::Tell() const
{
    if (impl_->pipeline) {
        return impl_->pipeline->Tell();
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

BgzfMetrics BamRawReader::GetMetrics() const
{
    if (impl_->pipeline) {
        BgzfMetrics m{impl_->pipeline->GetMetrics()};
        return m;
    }

    // Sync mode: only RecordsConsumed is meaningful
    BgzfMetrics m{};
    m.RecordsConsumed = impl_->recordsRead;
    return m;
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
