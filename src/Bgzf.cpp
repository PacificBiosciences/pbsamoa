#include <pbsamoa/core/Bgzf.hpp>

#include "BinaryUtils.hpp"
#include "LibdeflateUtils.hpp"
#include "WriterUtils.hpp"

#include <pbcopper/parallel/ThreadPool.h>
#include <rigtorp/SPSCQueue.hpp>

#include <libdeflate.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <exception>
#include <fstream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <cstddef>
#include <cstdint>
#include <cstring>

namespace PacBio {
namespace Samoa {

namespace detail {

using DecompressorPtr = LibdeflateDecompressorPtr;
using CompressorPtr = LibdeflateCompressorPtr;

/// \brief Entry describing one compressed BGZF payload within a batch buffer.
struct CompressedEntry
{
    std::size_t offset;        // offset into batch's contiguous compressed buffer
    std::size_t size;          // DEFLATE payload size (no BGZF framing)
    std::uint32_t isize;       // expected decompressed size from BGZF trailer
    std::uint64_t fileOffset;  // original file position for VirtualOffset tracking
};

/// \brief Batch of decompressed BGZF blocks for producer-consumer pipeline.
struct DecompressedBatch
{
    std::vector<std::byte> data;  // contiguous decompressed output (all blocks concatenated)
};

}  // namespace detail

namespace {

// --- Common BGZF constants ---

constexpr std::size_t GZIP_HEADER_SIZE{10};
constexpr std::uint8_t GZIP_ID1{31};
constexpr std::uint8_t GZIP_ID2{139};
constexpr std::uint8_t GZIP_CM_DEFLATE{8};
constexpr std::uint8_t GZIP_FLG_FEXTRA{4};

constexpr std::uint8_t BGZF_SI1{66};  // 'B'
constexpr std::uint8_t BGZF_SI2{67};  // 'C'
constexpr std::uint32_t GZIP_TRAILER_SIZE{8};

// Shared BGZF header size (used by reader and pipeline)
constexpr std::streamsize BGZF_HEADER_SIZE{18};

// Reader constants
constexpr std::size_t MAX_BGZF_BLOCK_SIZE{65536U + 256U};

// Writer constants
constexpr std::size_t BGZF_MAX_BLOCK_SIZE{65536U};
constexpr std::size_t MAX_UNCOMPRESSED_SIZE{0xFF00U};  // 65280

// Pipeline constants
constexpr std::size_t MAX_COMPRESSED_BLOCK_SIZE{65536U + 256U};
constexpr std::size_t OUTPUT_QUEUE_CAPACITY{4096};
constexpr std::size_t RECORD_BLOCK_SIZE_FIELD{4};
constexpr std::size_t BLOCKS_PER_BATCH{32};

}  // namespace

// =============================================================================
// ParseBgzfBlockHeader / DecompressBgzfBlock / IsBgzfEofMarker
// =============================================================================

std::optional<BgzfBlockInfo> ParseBgzfBlockHeader(std::span<const std::byte> data)
{
    if (std::ssize(data) < 18) {
        return std::nullopt;
    }

    const std::byte* pointer{std::data(data)};

    if ((std::to_integer<std::uint8_t>(pointer[0]) != GZIP_ID1) ||
        (std::to_integer<std::uint8_t>(pointer[1]) != GZIP_ID2) ||
        (std::to_integer<std::uint8_t>(pointer[2]) != GZIP_CM_DEFLATE)) {
        return std::nullopt;
    }

    const std::uint8_t flag{std::to_integer<std::uint8_t>(pointer[3])};
    if ((flag & GZIP_FLG_FEXTRA) == 0) {
        return std::nullopt;
    }

    const std::uint16_t xlen{ReadU16LE(pointer + GZIP_HEADER_SIZE)};
    const std::size_t extraEnd{GZIP_HEADER_SIZE + 2U + xlen};

    if (std::size(data) < extraEnd) {
        return std::nullopt;
    }

    std::size_t offset{GZIP_HEADER_SIZE + 2U};
    while ((offset + 4U) <= extraEnd) {
        const std::uint8_t si1{std::to_integer<std::uint8_t>(pointer[offset])};
        const std::uint8_t si2{std::to_integer<std::uint8_t>(pointer[offset + 1U])};
        const std::uint16_t slen{ReadU16LE(pointer + offset + 2U)};

        if ((si1 == BGZF_SI1) && (si2 == BGZF_SI2) && (slen == 2U)) {
            if ((offset + 6U) > extraEnd) {
                return std::nullopt;
            }

            const std::uint16_t bsize{ReadU16LE(pointer + offset + 4U)};
            const std::uint32_t blockSize{bsize + 1U};
            const std::uint32_t cdataOffset = extraEnd;

            if (blockSize < (cdataOffset + GZIP_TRAILER_SIZE)) {
                return std::nullopt;
            }

            const std::uint32_t cdataSize{blockSize - cdataOffset - GZIP_TRAILER_SIZE};

            return BgzfBlockInfo{
                .blockSize = blockSize,
                .compressedDataOffset = cdataOffset,
                .compressedDataSize = cdataSize,
            };
        }

        offset += 4U + slen;
    }

    return std::nullopt;
}

std::optional<std::size_t> DecompressBgzfBlock(std::span<const std::byte> blockData,
                                               BgzfBlockInfo info, std::span<std::byte> output)
{
    if (std::size(blockData) < info.blockSize) {
        return std::nullopt;
    }

    const std::uint32_t isize{ReadU32LE(std::data(blockData) + info.blockSize - 4U)};
    if (isize == 0U) {
        return 0U;
    }

    if (std::size(output) < isize) {
        return std::nullopt;
    }

    const detail::DecompressorPtr decompressor{libdeflate_alloc_decompressor()};
    if (!decompressor) {
        return std::nullopt;
    }

    std::size_t actualOut{0};
    const libdeflate_result result{libdeflate_deflate_decompress(
        decompressor.get(), std::data(blockData) + info.compressedDataOffset,
        info.compressedDataSize, std::data(output), isize, &actualOut)};
    if (result != LIBDEFLATE_SUCCESS) {
        return std::nullopt;
    }

    return actualOut;
}

bool IsBgzfEofMarker(std::span<const std::byte> data)
{
    if (std::size(data) < std::size(BGZF_EOF_MARKER)) {
        return false;
    }

    return std::ranges::equal(data.first(std::size(BGZF_EOF_MARKER)), BGZF_EOF_MARKER);
}

/// \brief Always-on atomic counters for BGZF reader pipeline introspection.
struct PipelineCounters
{
    // IO stage
    std::atomic<std::uint64_t> bytesRead{0};
    std::atomic<std::uint64_t> blocksRead{0};
    std::atomic<std::uint64_t> bytesDecompressed{0};
    std::atomic<std::uint64_t> ioReadNs{0};
    std::atomic<std::uint64_t> ioStalls{0};

    // Decompression workers (summed across all workers)
    std::atomic<std::uint64_t> decompressNs{0};

    // Consumer thread
    std::atomic<std::uint64_t> recordsProduced{0};
    std::atomic<std::uint64_t> consumerStalls{0};
    std::atomic<std::uint64_t> recordParseNs{0};

    // Reader (caller thread)
    std::atomic<std::uint64_t> recordsConsumed{0};
    std::atomic<std::uint64_t> readerStalls{0};
};

void NotifyReader(std::mutex& readyMutex, std::condition_variable& readyCv)
{
    std::unique_lock lock{readyMutex};
    lock.unlock();
    readyCv.notify_one();
}

void ParseBufferedRecords(std::vector<std::byte>& accumulator, std::size_t& pos,
                          rigtorp::SPSCQueue<RawRecord>& outputQueue,
                          const std::stop_token& stopToken, PipelineCounters& counters,
                          std::mutex& readyMutex, std::condition_variable& readyCv)
{
    const auto parseStart{std::chrono::steady_clock::now()};
    bool pushed{false};

    while (pos + RECORD_BLOCK_SIZE_FIELD <= std::size(accumulator)) {
        const std::uint32_t blockSize{ReadU32LE(std::data(accumulator) + pos)};
        if (blockSize == 0) {
            pos += RECORD_BLOCK_SIZE_FIELD;
            continue;
        }

        const std::size_t totalRecordBytes{RECORD_BLOCK_SIZE_FIELD + blockSize};
        if (pos + totalRecordBytes > std::size(accumulator)) {
            break;  // Need more data
        }

        // Create owning view from record bytes (after block_size field)
        RawRecord view{std::span<const std::byte>{
            std::data(accumulator) + pos + RECORD_BLOCK_SIZE_FIELD, blockSize}};

        // Push to SPSC — spin while queue is full (reader drains on another core).
        while (!outputQueue.try_push(std::move(view))) {
            if (stopToken.stop_requested()) {
                counters.recordParseNs.fetch_add(
                    static_cast<std::uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                                                   std::chrono::steady_clock::now() - parseStart)
                                                   .count()),
                    std::memory_order_relaxed);
                return;
            }
            counters.consumerStalls.fetch_add(1, std::memory_order_relaxed);
        }

        counters.recordsProduced.fetch_add(1, std::memory_order_relaxed);
        pushed = true;
        pos += totalRecordBytes;
    }

    // Compact accumulator: remove consumed bytes
    if (pos > 0) {
        accumulator.erase(std::begin(accumulator),
                          std::begin(accumulator) + static_cast<std::ptrdiff_t>(pos));
        pos = 0;
    }

    counters.recordParseNs.fetch_add(
        static_cast<std::uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                                       std::chrono::steady_clock::now() - parseStart)
                                       .count()),
        std::memory_order_relaxed);

    if (pushed) {
        NotifyReader(readyMutex, readyCv);
    }
}

struct DecompressedBatchConsumer
{
    std::vector<std::byte>& accumulator;
    std::size_t& pos;
    rigtorp::SPSCQueue<RawRecord>& outputQueue;
    std::stop_token stopToken;
    PipelineCounters& counters;
    std::mutex& readyMutex;
    std::condition_variable& readyCv;

    void operator()(detail::DecompressedBatch batch) const
    {
        accumulator.insert(accumulator.end(), batch.data.begin(), batch.data.end());
        ParseBufferedRecords(accumulator, pos, outputQueue, stopToken, counters, readyMutex,
                             readyCv);
    }
};

struct BgzfPipelineState
{
    struct DecompressBatchTask
    {
        std::vector<detail::CompressedEntry> entries;
        std::vector<std::byte> compressedData;
        std::vector<detail::DecompressorPtr>* decompressors;
        PipelineCounters* counters;

        detail::DecompressedBatch operator()(Parallel::ThreadIndex threadIdx) const;
    };

    std::filesystem::path filePath;
    std::size_t numWorkers;
    bool hasEofMarker{false};

    // Sync decompression (for ReadBlock and header parsing)
    std::ifstream syncFile;
    detail::DecompressorPtr syncDecompressor;
    std::array<std::byte, MAX_COMPRESSED_BLOCK_SIZE> syncCompressed;
    std::uint64_t syncBlockOffset{0};

    // Header
    SamHeader header;
    std::vector<std::byte> initialBytes;
    bool headerParsed{false};
    std::uint64_t headerEndFileOffset{0};

    // Pipeline (created after ParseHeader when numWorkers > 0)
    std::unique_ptr<Parallel::ThreadPool<detail::DecompressedBatch>> pool;
    std::jthread ioThread;
    std::jthread consumerThread;
    std::unique_ptr<rigtorp::SPSCQueue<RawRecord>> outputQueue;
    std::mutex readyMutex;
    std::condition_variable readyCv;
    std::atomic<bool> done{false};
    std::atomic<bool> pipelineError{false};
    mutable std::mutex errorMutex;
    std::exception_ptr errorPtr;

    // Per-worker decompressors (indexed by ThreadIndex)
    std::vector<detail::DecompressorPtr> decompressors;

    // Always-on metrics
    PipelineCounters counters;

    explicit BgzfPipelineState(const std::filesystem::path& path, std::size_t workers);
    ~BgzfPipelineState();

    BgzfPipelineState(const BgzfPipelineState&) = delete;
    BgzfPipelineState& operator=(const BgzfPipelineState&) = delete;

    void ParseHeader();
    void StartPipeline();
    void StopPipeline();
    static void RunIoLoop(std::stop_token stopToken, BgzfPipelineState* self);
    static void RunConsumerLoop(std::stop_token stopToken, BgzfPipelineState* self);
    void IoLoop(std::stop_token stopToken);
    void ConsumerLoop(std::stop_token stopToken);
    std::optional<std::size_t> ReadBlockSync(std::span<std::byte> buffer);
    const SamHeader& Header() const;
    std::optional<RawRecord> TryPopOutputRecord();
    bool OutputQueueReady() const;
    std::optional<RawRecord> ReadRecord();
    void Seek(VirtualOffset offset);
    VirtualOffset Tell() const;
    BgzfMetrics GetMetrics() const;
};

// =============================================================================
// BgzfReader
// =============================================================================

enum class BgzfReaderMode : std::uint8_t
{
    BLOCK_ONLY,
    SYNC_BAM,
    PARALLEL_BAM,
};

struct BgzfReader::Impl
{
    BgzfReaderMode mode{BgzfReaderMode::BLOCK_ONLY};

    // Shared sync decompression state (block-only and sync BAM mode)
    std::ifstream syncFile{};
    std::vector<std::byte> compressedBuf{};
    std::uint64_t syncBlockOffset{0};
    bool hasEofMarker{false};
    detail::DecompressorPtr syncDecompressor{};

    // Sync BAM mode state
    SamHeader syncHeader{};
    std::vector<std::byte> recordBuf{};
    std::size_t recordBufSize{0};
    std::size_t recordBufPos{0};
    bool syncEof{false};

    // Parallel BAM mode state
    std::unique_ptr<BgzfPipelineState> pipeline{};

    // Records returned to caller (sync modes only; pipeline tracks internally)
    std::atomic<std::uint64_t> recordsConsumed{0};

    explicit Impl(const std::filesystem::path& path);
    Impl(const std::filesystem::path& path, std::size_t numWorkers);

    Impl(const Impl&) = delete;
    Impl& operator=(const Impl&) = delete;

    void OpenSyncFile(const std::filesystem::path& path);
    std::optional<std::size_t> ReadBlockSync(std::span<std::byte> buffer);
    void ParseHeaderSync();
    bool RefillRecordBuffer();
    std::optional<RawRecord> ReadRecordSync();

    std::optional<std::size_t> ReadBlock(std::span<std::byte> buffer);
    void Seek(VirtualOffset offset);
    VirtualOffset Tell() const;
    bool HasEofMarker() const;
    const SamHeader& Header() const;
    std::optional<RawRecord> ReadRecord();
    BgzfMetrics GetMetrics() const;
};

BgzfReader::BgzfReader(const std::filesystem::path& path) : impl_{std::make_unique<Impl>(path)} {}

BgzfReader::BgzfReader(const std::filesystem::path& path, std::size_t numWorkers)
    : impl_{std::make_unique<Impl>(path, numWorkers)}
{
}

BgzfReader::~BgzfReader() = default;
BgzfReader::BgzfReader(BgzfReader&&) noexcept = default;
BgzfReader& BgzfReader::operator=(BgzfReader&&) noexcept = default;

std::optional<std::size_t> BgzfReader::ReadBlock(std::span<std::byte> buffer)
{
    return impl_->ReadBlock(buffer);
}

void BgzfReader::Seek(VirtualOffset offset) { impl_->Seek(offset); }

bool BgzfReader::HasEofMarker() const { return impl_->HasEofMarker(); }

VirtualOffset BgzfReader::Tell() const { return impl_->Tell(); }

const SamHeader& BgzfReader::Header() const { return impl_->Header(); }

std::optional<RawRecord> BgzfReader::ReadRecord() { return impl_->ReadRecord(); }

BgzfMetrics BgzfReader::GetMetrics() const { return impl_->GetMetrics(); }

// =============================================================================
// BgzfWriter
// =============================================================================

namespace detail {

enum class WriteItemKind : std::uint8_t
{
    DATA,
    CLOSE,
};

struct WriteItem
{
    WriteItemKind kind{WriteItemKind::DATA};
    std::vector<std::byte> data{};
    PendingCallback callback{};

    WriteItem() = default;

    explicit WriteItem(WriteItemKind k) : kind{k} {}

    WriteItem(WriteItemKind k, std::vector<std::byte>&& d, PendingCallback&& cb)
        : kind{k}, data{std::move(d)}, callback{std::move(cb)}
    {
    }
};

struct UncompressedBlock
{
    std::vector<std::byte> data;
    std::vector<PendingCallback> callbacks;
};

struct CompressedBatch
{
    struct Block
    {
        std::vector<std::uint8_t> cdata;
        std::uint32_t crc32{0};
        std::uint32_t isize{0};
        std::vector<PendingCallback> callbacks;
    };

    std::vector<Block> blocks;
};

struct WritePipelineCounters
{
    std::atomic<std::uint64_t> callerStalls{0};
    std::atomic<std::uint64_t> packerStalls{0};
    std::atomic<std::uint64_t> compressNs{0};
    std::atomic<std::uint64_t> bytesCompressed{0};
    std::atomic<std::uint64_t> blocksWritten{0};
    std::atomic<std::uint64_t> ioWriteNs{0};
    std::atomic<std::uint64_t> writerStalls{0};
    std::atomic<std::uint64_t> callbackNs{0};
};

}  // namespace detail

struct BgzfWriter::Impl
{
    struct CompressedBatchWriter
    {
        Impl* self;

        void operator()(detail::CompressedBatch batch) const;
    };

    struct CompressBatchTask
    {
        std::vector<detail::UncompressedBlock> blocks;
        std::vector<detail::CompressorPtr>* compressors;
        detail::WritePipelineCounters* counters;

        detail::CompressedBatch operator()(Parallel::ThreadIndex threadIdx);
    };

    BgzfWriterConfig config;
    std::filesystem::path finalPath;
    std::filesystem::path writePath;
    std::ofstream file{};
    std::uint64_t compressedOffset{0};
    std::unique_ptr<Parallel::ThreadPool<detail::CompressedBatch>> pool;
    std::unique_ptr<rigtorp::SPSCQueue<detail::WriteItem>> inputQueue;
    std::jthread packerThread;
    std::jthread ioWriterThread;
    std::vector<detail::CompressorPtr> compressors;
    IndexCallbackFn indexCallback{};
    mutable std::mutex callbackMutex;
    std::atomic<bool> pipelineError{false};
    std::exception_ptr errorPtr;
    mutable std::mutex errorMutex;
    detail::WritePipelineCounters counters;
    bool closed{false};

    explicit Impl(const std::filesystem::path& path, const BgzfWriterConfig& cfg)
        : config{cfg}, finalPath{path}
    {
        if (config.CompressionLevel < 1 || config.CompressionLevel > 12) {
            throw std::invalid_argument{"BgzfWriter: compression level must be in [1, 12]"};
        }
        if (config.BgzfWorkers == 0) {
            throw std::invalid_argument{"BgzfWriter: BgzfWorkers must be > 0"};
        }
        if (config.InputQueueCapacity == 0) {
            throw std::invalid_argument{"BgzfWriter: InputQueueCapacity must be > 0"};
        }
        if (config.BlocksPerBatch <= 0) {
            throw std::invalid_argument{"BgzfWriter: BlocksPerBatch must be > 0"};
        }

        writePath = detail::ResolveWritePath(finalPath, config.UseTempFile);
        file.open(writePath, std::ios::binary | std::ios::trunc);
        if (!file.is_open()) {
            throw std::runtime_error{"Cannot open file for writing: " + writePath.string()};
        }

        compressors.reserve(config.BgzfWorkers);
        for (std::size_t i{0}; i < config.BgzfWorkers; ++i) {
            compressors.emplace_back(libdeflate_alloc_compressor(config.CompressionLevel));
            if (!compressors.back()) {
                throw std::runtime_error{"Failed to create libdeflate compressor"};
            }
        }

        inputQueue =
            std::make_unique<rigtorp::SPSCQueue<detail::WriteItem>>(config.InputQueueCapacity);
        pool = std::make_unique<Parallel::ThreadPool<detail::CompressedBatch>>(
            Parallel::ThreadPool<detail::CompressedBatch>::Config{
                .NumThreads = config.BgzfWorkers,
                .QueueMultiplier = 3,
                .EnableMetrics = true,
            });

        try {
            // Consumer must run before Submit() starts in producer-consumer mode.
            ioWriterThread = std::jthread{&Impl::RunIoWriterLoop, this};
            packerThread = std::jthread{&Impl::RunPackerLoop, this};
        } catch (...) {
            ShutdownPipelineNoThrow();
            throw;
        }
    }

    Impl(const Impl&) = delete;
    Impl& operator=(const Impl&) = delete;

    static void RunIoWriterLoop(std::stop_token stopToken, Impl* self);
    static void RunPackerLoop(std::stop_token stopToken, Impl* self);

    void FinalizePoolNoThrow() noexcept
    {
        if (!pool) {
            return;
        }
        try {
            pool->Finalize();
        } catch (...) {
        }
    }

    void RequestStopThreads() noexcept
    {
        if (packerThread.joinable()) {
            packerThread.request_stop();
        }
        if (ioWriterThread.joinable()) {
            ioWriterThread.request_stop();
        }
    }

    void JoinThreadsNoThrow() noexcept
    {
        if (packerThread.joinable()) {
            try {
                packerThread.join();
            } catch (...) {
            }
        }
        if (ioWriterThread.joinable()) {
            try {
                ioWriterThread.join();
            } catch (...) {
            }
        }
    }

    void ShutdownPipelineNoThrow() noexcept
    {
        RequestStopThreads();
        // ConsumeWith() can block until finalization; unblock before join.
        FinalizePoolNoThrow();
        JoinThreadsNoThrow();
    }

    void SetError(std::exception_ptr e) noexcept
    {
        const std::lock_guard lock{errorMutex};
        if (pipelineError.load(std::memory_order_relaxed)) {
            return;
        }
        errorPtr = std::move(e);
        pipelineError.store(true, std::memory_order_release);
    }

    void RethrowIfError() const
    {
        if (!pipelineError.load(std::memory_order_acquire)) {
            return;
        }
        const std::lock_guard lock{errorMutex};
        if (errorPtr) {
            std::rethrow_exception(errorPtr);
        }
        throw std::runtime_error{"BgzfWriter: background pipeline failed"};
    }

    void ResetCurrentBlock(detail::UncompressedBlock& currentBlock) const
    {
        currentBlock = detail::UncompressedBlock{};
        currentBlock.data.reserve(MAX_UNCOMPRESSED_SIZE);
    }

    void SubmitBatch(std::vector<detail::UncompressedBlock>& batch, std::size_t batchLimit)
    {
        if (std::empty(batch)) {
            return;
        }
        std::vector<detail::UncompressedBlock> blocks{};
        blocks.swap(batch);
        batch.reserve(batchLimit);

        pool->Submit(CompressBatchTask{std::move(blocks), &compressors, &counters});
    }

    void FlushBatchIfFull(std::vector<detail::UncompressedBlock>& batch, std::size_t batchLimit)
    {
        if (std::size(batch) >= batchLimit) {
            SubmitBatch(batch, batchLimit);
        }
    }

    void SealCurrentBlock(std::vector<detail::UncompressedBlock>& batch,
                          detail::UncompressedBlock& currentBlock, std::size_t batchLimit)
    {
        if (std::empty(currentBlock.data)) {
            return;
        }

        batch.push_back(std::move(currentBlock));
        ResetCurrentBlock(currentBlock);
        FlushBatchIfFull(batch, batchLimit);
    }

    void PackerLoop(std::stop_token stopToken)
    {
        try {
            const std::size_t batchLimit = config.BlocksPerBatch;
            std::vector<detail::UncompressedBlock> batch{};
            batch.reserve(batchLimit);
            std::size_t emptyPollCount{0};

            detail::UncompressedBlock currentBlock{};
            ResetCurrentBlock(currentBlock);

            while (!stopToken.stop_requested()) {
                detail::WriteItem* front{inputQueue->front()};
                if (!front) {
                    counters.packerStalls.fetch_add(1, std::memory_order_relaxed);
                    ++emptyPollCount;
                    if (emptyPollCount < 64) {
                        std::this_thread::yield();
                    } else {
                        std::this_thread::sleep_for(std::chrono::microseconds{50});
                    }
                    continue;
                }
                emptyPollCount = 0;

                detail::WriteItem item{std::move(*front)};
                inputQueue->pop();

                if (item.kind == detail::WriteItemKind::CLOSE) {
                    SealCurrentBlock(batch, currentBlock, batchLimit);
                    SubmitBatch(batch, batchLimit);
                    break;
                }
                if (item.kind != detail::WriteItemKind::DATA || std::empty(item.data)) {
                    continue;
                }

                const bool fitsInSingleBlock{std::size(item.data) <= MAX_UNCOMPRESSED_SIZE};
                if (fitsInSingleBlock && !std::empty(currentBlock.data) &&
                    (std::size(currentBlock.data) + std::size(item.data) > MAX_UNCOMPRESSED_SIZE)) {
                    SealCurrentBlock(batch, currentBlock, batchLimit);
                }

                std::size_t dataPos{0};
                bool callbackAssigned{false};
                while (dataPos < std::size(item.data)) {
                    if (std::size(currentBlock.data) >= MAX_UNCOMPRESSED_SIZE) {
                        SealCurrentBlock(batch, currentBlock, batchLimit);
                    }

                    const std::size_t space{MAX_UNCOMPRESSED_SIZE - std::size(currentBlock.data)};
                    const std::size_t remaining{std::size(item.data) - dataPos};
                    const std::size_t toCopy{std::ranges::min(space, remaining)};

                    if (item.callback.active && !callbackAssigned) {
                        item.callback.withinBlockOffset =
                            static_cast<std::uint16_t>(std::size(currentBlock.data));
                        currentBlock.callbacks.push_back(std::move(item.callback));
                        callbackAssigned = true;
                    }

                    currentBlock.data.insert(
                        std::ranges::end(currentBlock.data),
                        std::ranges::begin(item.data) + static_cast<std::ptrdiff_t>(dataPos),
                        std::ranges::begin(item.data) +
                            static_cast<std::ptrdiff_t>(dataPos + toCopy));
                    dataPos += toCopy;
                }
            }

            pool->Finalize();
        } catch (...) {
            SetError(std::current_exception());
            try {
                pool->Finalize();
            } catch (...) {
            }
        }
    }

    void WriteFrame(const detail::CompressedBatch::Block& block)
    {
        const auto t0{std::chrono::steady_clock::now()};

        const std::uint16_t xlen{6U};
        const std::uint32_t blockSize{10U + 2U + xlen +
                                      static_cast<std::uint32_t>(std::size(block.cdata)) + 8U};
        if (blockSize > 65536U) {
            throw std::runtime_error{"BgzfWriter: block size exceeds BGZF maximum"};
        }
        const std::uint16_t bsize = blockSize - 1U;

        const std::array<std::uint8_t, 18> frameHeader{
            GZIP_ID1,
            GZIP_ID2,
            GZIP_CM_DEFLATE,
            GZIP_FLG_FEXTRA,  // ID1, ID2, CM, FLG
            0U,
            0U,
            0U,
            0U,  // MTIME
            0U,
            0U,  // XFL, OS
            static_cast<std::uint8_t>(xlen & 0xFFU),
            static_cast<std::uint8_t>(xlen >> 8U),
            BGZF_SI1,
            BGZF_SI2,
            2U,
            0U,  // BC subfield ID + length
            static_cast<std::uint8_t>(bsize & 0xFFU),
            static_cast<std::uint8_t>(bsize >> 8U),
        };
        file.write(reinterpret_cast<const char*>(std::data(frameHeader)),
                   static_cast<std::streamsize>(std::size(frameHeader)));
        file.write(reinterpret_cast<const char*>(std::data(block.cdata)),
                   static_cast<std::streamsize>(std::size(block.cdata)));

        const std::array<std::uint8_t, 8> trailer{
            static_cast<std::uint8_t>(block.crc32 & 0xFFU),
            static_cast<std::uint8_t>((block.crc32 >> 8U) & 0xFFU),
            static_cast<std::uint8_t>((block.crc32 >> 16U) & 0xFFU),
            static_cast<std::uint8_t>((block.crc32 >> 24U) & 0xFFU),
            static_cast<std::uint8_t>(block.isize & 0xFFU),
            static_cast<std::uint8_t>((block.isize >> 8U) & 0xFFU),
            static_cast<std::uint8_t>((block.isize >> 16U) & 0xFFU),
            static_cast<std::uint8_t>((block.isize >> 24U) & 0xFFU),
        };
        file.write(reinterpret_cast<const char*>(std::data(trailer)),
                   static_cast<std::streamsize>(std::size(trailer)));

        if (!file.good()) {
            throw std::runtime_error{"BgzfWriter: write failed"};
        }

        compressedOffset += blockSize;
        const auto t1{std::chrono::steady_clock::now()};
        counters.ioWriteNs.fetch_add(
            static_cast<std::uint64_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()),
            std::memory_order_relaxed);
    }

    void FireCallbacks(const detail::CompressedBatch::Block& block, std::uint64_t blockFileOffset)
    {
        IndexCallbackFn callback{};
        {
            const std::lock_guard lock{callbackMutex};
            callback = indexCallback;
        }
        if (!callback || std::empty(block.callbacks)) {
            return;
        }

        const auto t0{std::chrono::steady_clock::now()};
        for (const auto& cb : block.callbacks) {
            const VirtualOffset offset{blockFileOffset, cb.withinBlockOffset};
            callback(static_cast<std::int64_t>(offset.Value()),
                     std::span<const std::byte>{cb.rawData});
        }
        const auto t1{std::chrono::steady_clock::now()};
        counters.callbackNs.fetch_add(
            static_cast<std::uint64_t>(
                std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()),
            std::memory_order_relaxed);
    }

    void WriteCompressedBatch(detail::CompressedBatch batch)
    {
        for (const auto& block : batch.blocks) {
            const std::uint64_t blockFileOffset{compressedOffset};
            WriteFrame(block);
            FireCallbacks(block, blockFileOffset);
            counters.blocksWritten.fetch_add(1, std::memory_order_relaxed);
        }
    }

    void IoWriterLoop(std::stop_token /*stopToken*/)
    {
        try {
            const CompressedBatchWriter writeBatch{this};
            while (true) {
                if (pool->GetMetrics().CurrentResultQueueDepth == 0U) {
                    counters.writerStalls.fetch_add(1, std::memory_order_relaxed);
                }
                const bool keepConsuming = pool->ConsumeWith(writeBatch);
                if (!keepConsuming) {
                    break;
                }
            }
        } catch (...) {
            SetError(std::current_exception());
        }
    }
};

void BgzfWriter::Impl::CompressedBatchWriter::operator()(detail::CompressedBatch batch) const
{
    self->WriteCompressedBatch(std::move(batch));
}

detail::CompressedBatch BgzfWriter::Impl::CompressBatchTask::operator()(
    Parallel::ThreadIndex threadIdx)
{
    detail::CompressedBatch result{};
    result.blocks.reserve(std::size(blocks));

    auto* compressor{(*compressors)[threadIdx.Value()].get()};
    std::vector<std::uint8_t> compressedBuffer(BGZF_MAX_BLOCK_SIZE);

    for (auto& block : blocks) {
        const auto compressStart{std::chrono::steady_clock::now()};
        const std::size_t compressedSize{
            libdeflate_deflate_compress(compressor, std::data(block.data), std::size(block.data),
                                        std::data(compressedBuffer), std::size(compressedBuffer))};
        if (compressedSize == 0U) {
            throw std::runtime_error{"libdeflate_deflate_compress failed"};
        }

        counters->compressNs.fetch_add(
            static_cast<std::uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                                           std::chrono::steady_clock::now() - compressStart)
                                           .count()),
            std::memory_order_relaxed);

        const std::uint32_t crc32{
            libdeflate_crc32(0, std::data(block.data), std::size(block.data))};
        const std::uint32_t isize{static_cast<std::uint32_t>(std::size(block.data))};
        result.blocks.push_back(detail::CompressedBatch::Block{
            .cdata = {std::begin(compressedBuffer), std::begin(compressedBuffer) + compressedSize},
            .crc32 = crc32,
            .isize = isize,
            .callbacks = std::move(block.callbacks),
        });
        counters->bytesCompressed.fetch_add(compressedSize, std::memory_order_relaxed);
    }

    return result;
}

void BgzfWriter::Impl::RunIoWriterLoop(std::stop_token stopToken, Impl* self)
{
    self->IoWriterLoop(stopToken);
}

void BgzfWriter::Impl::RunPackerLoop(std::stop_token stopToken, Impl* self)
{
    self->PackerLoop(stopToken);
}

BgzfWriter::BgzfWriter(const std::filesystem::path& path, const BgzfWriterConfig& config)
    : impl_{std::make_unique<Impl>(path, config)}
{
}

BgzfWriter::~BgzfWriter()
{
    if (impl_ && !impl_->closed) {
        try {
            Close();
        } catch (...) {
        }
    }
}

BgzfWriter::BgzfWriter(BgzfWriter&&) noexcept = default;
BgzfWriter& BgzfWriter::operator=(BgzfWriter&&) noexcept = default;

void BgzfWriter::SetCallback(IndexCallbackFn callback)
{
    const std::lock_guard lock{impl_->callbackMutex};
    impl_->indexCallback = std::move(callback);
}

void BgzfWriter::Write(std::span<const std::byte> data)
{
    if (std::empty(data)) {
        return;
    }
    Write(std::vector<std::byte>{std::ranges::begin(data), std::ranges::end(data)},
          PendingCallback{});
}

void BgzfWriter::Write(std::vector<std::byte>&& data) { Write(std::move(data), PendingCallback{}); }

void BgzfWriter::Write(std::span<const std::byte> data, const PendingCallback& callback)
{
    if (std::empty(data)) {
        return;
    }
    Write(std::vector<std::byte>{std::ranges::begin(data), std::ranges::end(data)},
          PendingCallback{callback});
}

void BgzfWriter::Write(std::vector<std::byte>&& data, PendingCallback callback)
{
    if (impl_->closed) {
        throw std::runtime_error{"BgzfWriter: Write() called after Close()"};
    }
    impl_->RethrowIfError();

    while (!impl_->inputQueue->try_emplace(detail::WriteItemKind::DATA, std::move(data),
                                           std::move(callback))) {
        impl_->RethrowIfError();
        impl_->counters.callerStalls.fetch_add(1, std::memory_order_relaxed);
        std::this_thread::yield();
    }
}

void BgzfWriter::Close()
{
    if (impl_->closed) {
        return;
    }
    // Enter terminal state immediately: after Close() starts, no further writes
    // are valid even if teardown later throws.
    impl_->closed = true;

    try {
        bool closeQueued{false};
        while (!closeQueued) {
            closeQueued = impl_->inputQueue->try_emplace(detail::WriteItemKind::CLOSE);
            if (closeQueued) {
                break;
            }
            if (impl_->pipelineError.load(std::memory_order_acquire)) {
                break;
            }
            std::this_thread::yield();
        }

        if (!closeQueued) {
            impl_->ShutdownPipelineNoThrow();
        } else {
            if (impl_->packerThread.joinable()) {
                impl_->packerThread.join();
            }
            if (impl_->ioWriterThread.joinable()) {
                impl_->ioWriterThread.join();
            }
        }

        impl_->RethrowIfError();

        impl_->file.write(reinterpret_cast<const char*>(std::data(BGZF_EOF_MARKER)),
                          static_cast<std::streamsize>(std::size(BGZF_EOF_MARKER)));
        if (!impl_->file.good()) {
            throw std::runtime_error{"BgzfWriter: write failed during EOF marker"};
        }
        impl_->file.close();
        if (impl_->file.fail()) {
            throw std::runtime_error{"BgzfWriter: close failed"};
        }

        if (impl_->config.UseTempFile) {
            detail::AtomicRename(impl_->writePath, impl_->finalPath, "BgzfWriter");
        }
    } catch (...) {
        impl_->ShutdownPipelineNoThrow();
        if (impl_->file.is_open()) {
            impl_->file.close();
        }
        throw;
    }
}

BgzfWriteMetrics BgzfWriter::GetMetrics() const
{
    BgzfWriteMetrics m{};
    m.CallerStalls = impl_->counters.callerStalls.load(std::memory_order_relaxed);
    m.PackerStalls = impl_->counters.packerStalls.load(std::memory_order_relaxed);
    m.CompressNs = impl_->counters.compressNs.load(std::memory_order_relaxed);
    m.BytesCompressed = impl_->counters.bytesCompressed.load(std::memory_order_relaxed);
    m.BlocksWritten = impl_->counters.blocksWritten.load(std::memory_order_relaxed);
    m.IoWriteNs = impl_->counters.ioWriteNs.load(std::memory_order_relaxed);
    m.WriterStalls = impl_->counters.writerStalls.load(std::memory_order_relaxed);
    m.CallbackNs = impl_->counters.callbackNs.load(std::memory_order_relaxed);

    if (impl_->pool) {
        const auto snap{impl_->pool->GetMetrics()};
        m.Pool.QueueDepth = snap.CurrentQueueDepth;
        m.Pool.PeakQueueDepth = snap.PeakQueueDepth;
        m.Pool.ActiveWorkers = snap.CurrentActiveTasks;
        m.Pool.PeakActiveWorkers = snap.PeakActiveTasks;
        m.Pool.ResultQueueDepth = snap.CurrentResultQueueDepth;
        m.Pool.PeakResultQueueDepth = snap.PeakResultQueueDepth;
    }
    return m;
}

// =============================================================================
// BgzfReader parallel state
// =============================================================================

/// \brief Helper: measure elapsed nanoseconds for a scope.
struct ScopedTimer
{
    std::atomic<std::uint64_t>& target;
    std::chrono::steady_clock::time_point start;

    explicit ScopedTimer(std::atomic<std::uint64_t>& t)
        : target{t}, start{std::chrono::steady_clock::now()}
    {
    }

    ~ScopedTimer()
    {
        const auto elapsed{std::chrono::steady_clock::now() - start};
        target.fetch_add(static_cast<std::uint64_t>(
                             std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count()),
                         std::memory_order_relaxed);
    }

    ScopedTimer(const ScopedTimer&) = delete;
    ScopedTimer& operator=(const ScopedTimer&) = delete;
};

BgzfPipelineState::BgzfPipelineState(const std::filesystem::path& path, std::size_t workers)
    : filePath{path}, numWorkers{workers}
{
    // Check EOF marker
    {
        std::ifstream file{path, std::ios::binary};
        if (!file.is_open()) {
            throw std::runtime_error{"Cannot open file: " + path.string()};
        }
        file.seekg(-static_cast<std::streamoff>(std::size(BGZF_EOF_MARKER)), std::ios::end);
        if (file.good()) {
            std::array<std::byte, 28> tail{};
            file.read(reinterpret_cast<char*>(std::data(tail)),
                      static_cast<std::streamsize>(std::size(tail)));
            hasEofMarker = IsBgzfEofMarker(tail);
        }
    }

    // Open sync file for ReadBlock and header parsing
    syncFile.open(path, std::ios::binary);
    if (!syncFile.is_open()) {
        throw std::runtime_error{"Cannot open file: " + path.string()};
    }
    syncDecompressor = detail::DecompressorPtr{libdeflate_alloc_decompressor()};
    if (!syncDecompressor) {
        throw std::runtime_error{"Failed to create libdeflate decompressor"};
    }
}

BgzfPipelineState::~BgzfPipelineState() { StopPipeline(); }

void BgzfPipelineState::ParseHeader()
{
    // Read BGZF blocks synchronously until we have the complete BAM header
    std::vector<std::byte> headerBuf(MAX_DECOMPRESSED_BLOCK_SIZE * 4);
    std::size_t headerLen{0};
    std::vector<std::byte> blockBuf(MAX_DECOMPRESSED_BLOCK_SIZE);

    while (true) {
        syncBlockOffset = static_cast<std::uint64_t>(syncFile.tellg());
        syncFile.read(reinterpret_cast<char*>(std::data(syncCompressed)), BGZF_HEADER_SIZE);
        if (syncFile.gcount() == 0) {
            break;  // EOF
        }
        if (syncFile.gcount() < BGZF_HEADER_SIZE) {
            throw std::runtime_error{"Truncated BGZF block while parsing BAM header"};
        }

        const std::optional<BgzfBlockInfo> info{ParseBgzfBlockHeader(
            std::span<const std::byte>{syncCompressed}.first(BGZF_HEADER_SIZE))};
        if (!info) {
            throw std::runtime_error{"Invalid BGZF header while parsing BAM header"};
        }

        const std::size_t remaining{info->blockSize - static_cast<std::size_t>(BGZF_HEADER_SIZE)};
        syncFile.read(reinterpret_cast<char*>(std::data(syncCompressed) + BGZF_HEADER_SIZE),
                      static_cast<std::streamsize>(remaining));
        if (static_cast<std::size_t>(syncFile.gcount()) < remaining) {
            throw std::runtime_error{"Truncated BGZF block while parsing BAM header"};
        }

        const std::uint32_t isize{ReadU32LE(std::data(syncCompressed) + info->blockSize - 4U)};
        if (isize == 0U) {
            break;  // EOF marker
        }

        // Decompress
        std::size_t actualOut{0};
        const libdeflate_result result{libdeflate_deflate_decompress(
            syncDecompressor.get(), std::data(syncCompressed) + info->compressedDataOffset,
            info->compressedDataSize, std::data(blockBuf), isize, &actualOut)};
        if (result != LIBDEFLATE_SUCCESS) {
            throw std::runtime_error{"BGZF decompression failed while parsing BAM header"};
        }

        // Append to header buffer
        if (headerLen + actualOut > std::size(headerBuf)) {
            headerBuf.resize(headerLen + actualOut + MAX_DECOMPRESSED_BLOCK_SIZE);
        }
        std::ranges::copy_n(std::data(blockBuf), actualOut, std::data(headerBuf) + headerLen);
        headerLen += actualOut;

        // Check if we have the complete header
        const std::size_t hdrSize{
            ComputeHeaderSize(std::span<const std::byte>{std::data(headerBuf), headerLen})};
        if (hdrSize > 0) {
            auto headerResult{SamHeader::FromBamHeaderBlock(
                std::span<const std::byte>{std::data(headerBuf), hdrSize})};
            if (!headerResult) {
                throw std::runtime_error{std::move(headerResult.error())};
            }
            header = std::move(*headerResult);

            // Save leftover bytes after header
            const std::size_t leftover{headerLen - hdrSize};
            if (leftover > 0) {
                initialBytes.assign(std::data(headerBuf) + hdrSize,
                                    std::data(headerBuf) + headerLen);
            }

            headerEndFileOffset = static_cast<std::uint64_t>(syncFile.tellg());
            headerParsed = true;

            // Start pipeline if parallel mode
            if (numWorkers > 0) {
                StartPipeline();
            }
            return;
        }
    }

    // EOF before complete header
    if (headerLen == 0) {
        throw std::runtime_error{"Empty BAM file: no BGZF blocks"};
    }
    const std::size_t hdrSize{
        ComputeHeaderSize(std::span<const std::byte>{std::data(headerBuf), headerLen})};
    if (hdrSize > 0) {
        auto headerResult{SamHeader::FromBamHeaderBlock(
            std::span<const std::byte>{std::data(headerBuf), hdrSize})};
        if (!headerResult) {
            throw std::runtime_error{std::move(headerResult.error())};
        }
        header = std::move(*headerResult);
        headerParsed = true;
    } else {
        throw std::runtime_error{"Incomplete BAM header"};
    }
}

void BgzfPipelineState::StartPipeline()
{
    // Create per-worker decompressors
    decompressors.clear();
    decompressors.reserve(numWorkers);
    for (std::size_t i{0}; i < numWorkers; ++i) {
        detail::DecompressorPtr d{libdeflate_alloc_decompressor()};
        if (!d) {
            throw std::runtime_error{"Failed to create libdeflate decompressor"};
        }
        decompressors.push_back(std::move(d));
    }

    // Create output queue
    outputQueue = std::make_unique<rigtorp::SPSCQueue<RawRecord>>(OUTPUT_QUEUE_CAPACITY);
    done.store(false, std::memory_order_release);
    pipelineError.store(false, std::memory_order_release);
    errorPtr = nullptr;

    // Create thread pool (producer-consumer mode) with metrics enabled
    pool = std::make_unique<Parallel::ThreadPool<detail::DecompressedBatch>>(
        Parallel::ThreadPool<detail::DecompressedBatch>::Config{
            .NumThreads = numWorkers,
            .QueueMultiplier = 3,
            .EnableMetrics = true,
        });

    // Start consumer thread first — must enter ConsumeWith() before IO thread
    // calls Submit()
    consumerThread = std::jthread{&BgzfPipelineState::RunConsumerLoop, this};

    // Start IO thread
    ioThread = std::jthread{&BgzfPipelineState::RunIoLoop, this};
}

void BgzfPipelineState::RunIoLoop(std::stop_token stopToken, BgzfPipelineState* self)
{
    self->IoLoop(stopToken);
}

void BgzfPipelineState::RunConsumerLoop(std::stop_token stopToken, BgzfPipelineState* self)
{
    self->ConsumerLoop(stopToken);
}

void BgzfPipelineState::StopPipeline()
{
    // 1. Signal IO thread to stop submitting
    if (ioThread.joinable()) {
        ioThread.request_stop();
    }

    // 2. Signal consumer to stop blocking on SPSC push
    if (consumerThread.joinable()) {
        consumerThread.request_stop();
    }

    // 3. Finalize pool — stops workers, unblocks Submit(), makes ConsumeWith()
    // return false
    if (pool) {
        try {
            pool->Finalize();
        } catch (...) {
            // Ignore exceptions during shutdown
        }
    }

    // 4. Join threads (jthread assignment to {} joins then destroys)
    ioThread = {};
    consumerThread = {};

    // 5. Drain SPSC
    if (outputQueue) {
        while (outputQueue->front()) {
            outputQueue->pop();
        }
    }

    // 6. Clean up
    pool.reset();
    outputQueue.reset();
    decompressors.clear();
}

void BgzfPipelineState::IoLoop(std::stop_token stopToken)
{
    try {
        std::ifstream file{filePath, std::ios::binary};
        if (!file.is_open()) {
            throw std::runtime_error{"IO thread: cannot open file: " + filePath.string()};
        }
        file.seekg(static_cast<std::streamoff>(headerEndFileOffset), std::ios::beg);

        std::array<std::byte, MAX_COMPRESSED_BLOCK_SIZE> compressed{};
        std::vector<std::byte> batchBuffer;
        std::vector<detail::CompressedEntry> batchEntries;
        batchEntries.reserve(BLOCKS_PER_BATCH);

        const std::size_t poolCapacity{numWorkers * 3};

        while (!stopToken.stop_requested()) {
            // Accumulate up to BLOCKS_PER_BATCH compressed payloads
            batchBuffer.clear();
            batchEntries.clear();

            std::uint64_t batchCompressedBytes{0};

            for (std::size_t b{0}; b < BLOCKS_PER_BATCH && !stopToken.stop_requested(); ++b) {
                const std::uint64_t blockOffset = file.tellg();

                // Time file reads
                const auto readStart{std::chrono::steady_clock::now()};

                // Read BGZF header
                file.read(reinterpret_cast<char*>(std::data(compressed)), BGZF_HEADER_SIZE);
                if (file.gcount() == 0) {
                    counters.ioReadNs.fetch_add(
                        static_cast<std::uint64_t>(
                            std::chrono::duration_cast<std::chrono::nanoseconds>(
                                std::chrono::steady_clock::now() - readStart)
                                .count()),
                        std::memory_order_relaxed);
                    break;  // EOF
                }
                if (file.gcount() < BGZF_HEADER_SIZE) {
                    counters.ioReadNs.fetch_add(
                        static_cast<std::uint64_t>(
                            std::chrono::duration_cast<std::chrono::nanoseconds>(
                                std::chrono::steady_clock::now() - readStart)
                                .count()),
                        std::memory_order_relaxed);
                    break;  // Truncated
                }

                const std::optional<BgzfBlockInfo> info{ParseBgzfBlockHeader(
                    std::span<const std::byte>{compressed}.first(BGZF_HEADER_SIZE))};
                if (!info) {
                    counters.ioReadNs.fetch_add(
                        static_cast<std::uint64_t>(
                            std::chrono::duration_cast<std::chrono::nanoseconds>(
                                std::chrono::steady_clock::now() - readStart)
                                .count()),
                        std::memory_order_relaxed);
                    break;  // Invalid header
                }

                // Read rest of block
                const std::size_t remaining{info->blockSize -
                                            static_cast<std::size_t>(BGZF_HEADER_SIZE)};
                file.read(reinterpret_cast<char*>(std::data(compressed) + BGZF_HEADER_SIZE),
                          static_cast<std::streamsize>(remaining));

                counters.ioReadNs.fetch_add(
                    static_cast<std::uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                                                   std::chrono::steady_clock::now() - readStart)
                                                   .count()),
                    std::memory_order_relaxed);

                if (static_cast<std::size_t>(file.gcount()) < remaining) {
                    break;  // Truncated
                }

                const std::uint32_t isize{ReadU32LE(std::data(compressed) + info->blockSize - 4U)};
                if (isize == 0U) {
                    break;  // EOF marker
                }

                // Append compressed payload to batch buffer
                const std::size_t batchOffset{std::size(batchBuffer)};
                const std::size_t cdataOffset{info->compressedDataOffset};
                const std::size_t cdataSize{info->compressedDataSize};
                batchBuffer.insert(std::end(batchBuffer), std::data(compressed) + cdataOffset,
                                   std::data(compressed) + cdataOffset + cdataSize);
                batchEntries.push_back(detail::CompressedEntry{
                    .offset = batchOffset,
                    .size = cdataSize,
                    .isize = isize,
                    .fileOffset = blockOffset,
                });
                batchCompressedBytes += info->blockSize;
            }

            if (std::empty(batchEntries)) {
                break;
            }

            // Update IO counters
            counters.bytesRead.fetch_add(batchCompressedBytes, std::memory_order_relaxed);
            counters.blocksRead.fetch_add(std::size(batchEntries), std::memory_order_relaxed);

            // Check if pool queue is at capacity before submitting (stall detection)
            {
                const auto poolMetrics{pool->GetMetrics()};
                if (poolMetrics.CurrentQueueDepth >= poolCapacity) {
                    counters.ioStalls.fetch_add(1, std::memory_order_relaxed);
                }
            }

            // Submit one task for the entire batch
            pool->Submit(DecompressBatchTask{std::move(batchEntries), std::move(batchBuffer),
                                             &decompressors, &counters});
        }
    } catch (...) {
        // Exceptions from Submit() are expected during shutdown (pool finalized).
        // Pool's fail-fast mechanism handles other errors.
    }

    // Signal no more work
    if (pool) {
        try {
            pool->Finalize();
        } catch (...) {
            // Idempotent — may have been finalized by StopPipeline
        }
    }
}

void BgzfPipelineState::ConsumerLoop(std::stop_token stopToken)
{
    try {
        std::vector<std::byte> accumulator{std::move(initialBytes)};
        std::size_t pos{0};

        // Process any records already present from initialBytes (header block
        // may contain record data when header + records fit in one BGZF block).
        ParseBufferedRecords(accumulator, pos, *outputQueue, stopToken, counters, readyMutex,
                             readyCv);

        const DecompressedBatchConsumer consumeBatch{
            accumulator, pos, *outputQueue, stopToken, counters, readyMutex, readyCv};
        while (pool->ConsumeWith(consumeBatch)) {
            if (stopToken.stop_requested()) {
                break;
            }
        }
    } catch (...) {
        {
            const std::lock_guard lock{errorMutex};
            errorPtr = std::current_exception();
        }
        pipelineError.store(true, std::memory_order_release);
    }

    done.store(true, std::memory_order_release);
    NotifyReader(readyMutex, readyCv);
}

std::optional<std::size_t> BgzfPipelineState::ReadBlockSync(std::span<std::byte> buffer)
{
    syncBlockOffset = static_cast<std::uint64_t>(syncFile.tellg());

    syncFile.read(reinterpret_cast<char*>(std::data(syncCompressed)), BGZF_HEADER_SIZE);
    if (syncFile.gcount() == 0) {
        return 0U;
    }
    if (syncFile.gcount() < BGZF_HEADER_SIZE) {
        return std::nullopt;
    }

    const std::optional<BgzfBlockInfo> info{
        ParseBgzfBlockHeader(std::span<const std::byte>{syncCompressed}.first(BGZF_HEADER_SIZE))};
    if (!info) {
        return std::nullopt;
    }

    const std::size_t blockSize{info->blockSize};
    const std::size_t remaining{blockSize - static_cast<std::size_t>(BGZF_HEADER_SIZE)};
    syncFile.read(reinterpret_cast<char*>(std::data(syncCompressed) + BGZF_HEADER_SIZE),
                  static_cast<std::streamsize>(remaining));
    if (static_cast<std::size_t>(syncFile.gcount()) < remaining) {
        return std::nullopt;
    }

    const std::uint32_t isize{ReadU32LE(std::data(syncCompressed) + blockSize - 4U)};
    if (isize == 0U) {
        return 0U;
    }
    if (std::size(buffer) < isize) {
        return std::nullopt;
    }

    std::size_t actualOut{0};
    const libdeflate_result result{libdeflate_deflate_decompress(
        syncDecompressor.get(), std::data(syncCompressed) + info->compressedDataOffset,
        info->compressedDataSize, std::data(buffer), isize, &actualOut)};
    if (result != LIBDEFLATE_SUCCESS) {
        return std::nullopt;
    }

    return actualOut;
}

const SamHeader& BgzfPipelineState::Header() const
{
    if (!headerParsed) {
        throw std::runtime_error{"BAM header not parsed"};
    }
    return header;
}

std::optional<RawRecord> BgzfPipelineState::TryPopOutputRecord()
{
    RawRecord* front{outputQueue->front()};
    if (!front) {
        return std::nullopt;
    }

    RawRecord result{std::move(*front)};
    outputQueue->pop();
    counters.recordsConsumed.fetch_add(1, std::memory_order_relaxed);
    return result;
}

std::optional<RawRecord> BgzfPipelineState::ReadRecord()
{
    if (!headerParsed) {
        throw std::runtime_error{"BAM header not parsed"};
    }
    if (!outputQueue) {
        throw std::runtime_error{"ReadRecord() requires numWorkers > 0"};
    }

    // Fast path: check SPSC queue without locking.
    while (true) {
        if (auto record{TryPopOutputRecord()}) {
            return std::move(*record);
        }

        // Queue empty — check termination before sleeping.
        if (pipelineError.load(std::memory_order_acquire)) {
            if (auto record{TryPopOutputRecord()}) {
                return std::move(*record);
            }
            {
                const std::lock_guard lock{errorMutex};
                std::rethrow_exception(errorPtr);
            }
        }

        if (done.load(std::memory_order_acquire)) {
            if (auto record{TryPopOutputRecord()}) {
                return std::move(*record);
            }
            return std::nullopt;
        }

        // Slow path: wait on condition variable for data/done/error.
        // The mutex lock creates a happens-before with the consumer's
        // lock/unlock in notifyReader(), preventing lost wakeups.
        counters.readerStalls.fetch_add(1, std::memory_order_relaxed);
        std::unique_lock lock{readyMutex};
        while (!OutputQueueReady()) {
            readyCv.wait(lock);
        }
    }
}

detail::DecompressedBatch BgzfPipelineState::DecompressBatchTask::operator()(
    Parallel::ThreadIndex threadIdx) const
{
    const auto decompressStart{std::chrono::steady_clock::now()};

    std::uint64_t totalIsize{0};
    for (const auto& entry : entries) {
        totalIsize += entry.isize;
    }

    detail::DecompressedBatch batch{};
    batch.data.resize(totalIsize);

    std::size_t outOffset{0};
    for (const auto& entry : entries) {
        std::size_t actualOut{0};
        const libdeflate_result result{libdeflate_deflate_decompress(
            (*decompressors)[threadIdx.Value()].get(), std::data(compressedData) + entry.offset,
            entry.size, std::data(batch.data) + outOffset, entry.isize, &actualOut)};
        if (result != LIBDEFLATE_SUCCESS) {
            throw std::runtime_error{"BGZF decompression failed"};
        }
        outOffset += actualOut;
    }
    batch.data.resize(outOffset);

    counters->decompressNs.fetch_add(
        static_cast<std::uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                                       std::chrono::steady_clock::now() - decompressStart)
                                       .count()),
        std::memory_order_relaxed);
    counters->bytesDecompressed.fetch_add(outOffset, std::memory_order_relaxed);

    return batch;
}

bool BgzfPipelineState::OutputQueueReady() const
{
    return outputQueue->front() || done.load(std::memory_order_acquire) ||
           pipelineError.load(std::memory_order_acquire);
}

void BgzfPipelineState::Seek(VirtualOffset offset)
{
    // Stop pipeline if running
    StopPipeline();

    // Reset sync file position
    syncFile.clear();
    // Within-block offset intentionally discarded; callers (BamRawReader::Query)
    // scan past unwanted records within the target block.
    syncFile.seekg(static_cast<std::streamoff>(offset.BlockOffset()), std::ios::beg);
    syncBlockOffset = offset.BlockOffset();

    // Update pipeline start offset
    headerEndFileOffset = offset.BlockOffset();
    initialBytes.clear();

    // Reset pipeline state — safe without lock: pipeline threads joined by StopPipeline()
    done.store(false, std::memory_order_relaxed);
    pipelineError.store(false, std::memory_order_relaxed);
    errorPtr = nullptr;

    // Restart pipeline if in parallel mode and header was parsed
    if ((numWorkers > 0) && headerParsed) {
        StartPipeline();
    }
}

VirtualOffset BgzfPipelineState::Tell() const { return VirtualOffset{syncBlockOffset, 0}; }

BgzfMetrics BgzfPipelineState::GetMetrics() const
{
    BgzfMetrics m{};

    // Throughput counters
    m.BytesRead = counters.bytesRead.load(std::memory_order_relaxed);
    m.BlocksRead = counters.blocksRead.load(std::memory_order_relaxed);
    m.BytesDecompressed = counters.bytesDecompressed.load(std::memory_order_relaxed);

    // ThreadPool metrics
    if (pool) {
        const auto poolSnap{pool->GetMetrics()};
        m.Pool.QueueDepth = poolSnap.CurrentQueueDepth;
        m.Pool.PeakQueueDepth = poolSnap.PeakQueueDepth;
        m.Pool.ActiveWorkers = poolSnap.CurrentActiveTasks;
        m.Pool.PeakActiveWorkers = poolSnap.PeakActiveTasks;
        m.Pool.ResultQueueDepth = poolSnap.CurrentResultQueueDepth;
        m.Pool.PeakResultQueueDepth = poolSnap.PeakResultQueueDepth;
    }

    // SPSC queue
    m.RecordsProduced = counters.recordsProduced.load(std::memory_order_relaxed);
    m.RecordsConsumed = counters.recordsConsumed.load(std::memory_order_relaxed);

    // Stall counters
    m.IoStalls = counters.ioStalls.load(std::memory_order_relaxed);
    m.ConsumerStalls = counters.consumerStalls.load(std::memory_order_relaxed);
    m.ReaderStalls = counters.readerStalls.load(std::memory_order_relaxed);

    // Timing
    m.IoReadNs = counters.ioReadNs.load(std::memory_order_relaxed);
    m.DecompressNs = counters.decompressNs.load(std::memory_order_relaxed);
    m.RecordParseNs = counters.recordParseNs.load(std::memory_order_relaxed);

    return m;
}

BgzfReader::Impl::Impl(const std::filesystem::path& path) : mode{BgzfReaderMode::BLOCK_ONLY}
{
    OpenSyncFile(path);
}

BgzfReader::Impl::Impl(const std::filesystem::path& path, std::size_t numWorkers)
{
    if (numWorkers == 0) {
        mode = BgzfReaderMode::SYNC_BAM;
        OpenSyncFile(path);
        ParseHeaderSync();
        return;
    }

    mode = BgzfReaderMode::PARALLEL_BAM;
    pipeline = std::make_unique<BgzfPipelineState>(path, numWorkers);
    pipeline->ParseHeader();
}

void BgzfReader::Impl::OpenSyncFile(const std::filesystem::path& path)
{
    syncFile.open(path, std::ios::binary);
    if (!syncFile.is_open()) {
        throw std::runtime_error{"Cannot open file: " + path.string()};
    }

    compressedBuf.resize(MAX_BGZF_BLOCK_SIZE);

    syncDecompressor = detail::DecompressorPtr{libdeflate_alloc_decompressor()};
    if (!syncDecompressor) {
        throw std::runtime_error{"Failed to create libdeflate decompressor"};
    }

    syncFile.seekg(-static_cast<std::streamoff>(std::size(BGZF_EOF_MARKER)), std::ios::end);
    if (syncFile.good()) {
        std::array<std::byte, 28> tail{};
        syncFile.read(reinterpret_cast<char*>(std::data(tail)),
                      static_cast<std::streamsize>(std::size(tail)));
        hasEofMarker = IsBgzfEofMarker(tail);
    }

    syncFile.clear();
    syncFile.seekg(0, std::ios::beg);
}

std::optional<std::size_t> BgzfReader::Impl::ReadBlockSync(std::span<std::byte> buffer)
{
    syncBlockOffset = syncFile.tellg();

    syncFile.read(reinterpret_cast<char*>(std::data(compressedBuf)), BGZF_HEADER_SIZE);
    if (syncFile.gcount() == 0) {
        return 0U;
    }
    if (syncFile.gcount() < BGZF_HEADER_SIZE) {
        return std::nullopt;
    }

    const std::optional<BgzfBlockInfo> info{
        ParseBgzfBlockHeader(std::span<const std::byte>{compressedBuf}.first(
            static_cast<std::size_t>(BGZF_HEADER_SIZE)))};
    if (!info) {
        return std::nullopt;
    }

    const std::size_t blockSize{info->blockSize};
    const std::size_t remaining{blockSize - BGZF_HEADER_SIZE};
    if (std::size(compressedBuf) < blockSize) {
        compressedBuf.resize(blockSize);
    }

    syncFile.read(reinterpret_cast<char*>(std::data(compressedBuf) + BGZF_HEADER_SIZE),
                  static_cast<std::streamsize>(remaining));
    if (static_cast<std::size_t>(syncFile.gcount()) < remaining) {
        return std::nullopt;
    }

    const std::uint32_t isize{ReadU32LE(std::data(compressedBuf) + blockSize - 4U)};
    if (isize == 0U) {
        return 0U;
    }

    if (std::size(buffer) < static_cast<std::size_t>(isize)) {
        return std::nullopt;
    }

    std::size_t actualOut{0};
    const libdeflate_result result{libdeflate_deflate_decompress(
        syncDecompressor.get(), std::data(compressedBuf) + info->compressedDataOffset,
        info->compressedDataSize, std::data(buffer), isize, &actualOut)};
    if (result != LIBDEFLATE_SUCCESS) {
        return std::nullopt;
    }

    return actualOut;
}

void BgzfReader::Impl::ParseHeaderSync()
{
    std::vector<std::byte> headerBuf(MAX_DECOMPRESSED_BLOCK_SIZE * 4);
    std::size_t headerLen{0};
    std::vector<std::byte> blockBuf(MAX_DECOMPRESSED_BLOCK_SIZE);

    while (true) {
        const std::optional<std::size_t> bytesRead{ReadBlockSync(std::span<std::byte>{blockBuf})};

        if (!bytesRead) {
            throw std::runtime_error{"Failed to read BGZF block while parsing BAM header"};
        }
        if (*bytesRead == 0) {
            break;
        }

        if (headerLen + *bytesRead > std::size(headerBuf)) {
            headerBuf.resize(headerLen + *bytesRead + MAX_DECOMPRESSED_BLOCK_SIZE);
        }

        std::ranges::copy_n(std::data(blockBuf), *bytesRead, std::data(headerBuf) + headerLen);
        headerLen += *bytesRead;

        const std::size_t hdrSize{
            ComputeHeaderSize(std::span<const std::byte>{std::data(headerBuf), headerLen})};
        if (hdrSize > 0) {
            auto headerResult{SamHeader::FromBamHeaderBlock(
                std::span<const std::byte>{std::data(headerBuf), hdrSize})};
            if (!headerResult) {
                throw std::runtime_error{std::move(headerResult.error())};
            }
            syncHeader = std::move(*headerResult);

            const std::size_t remaining{headerLen - hdrSize};
            recordBuf.resize(std::ranges::max(remaining, MAX_DECOMPRESSED_BLOCK_SIZE));
            if (remaining > 0) {
                std::ranges::copy_n(std::data(headerBuf) + hdrSize, remaining,
                                    std::data(recordBuf));
            }
            recordBufSize = remaining;
            recordBufPos = 0;
            syncEof = false;
            return;
        }
    }

    if (headerLen == 0) {
        throw std::runtime_error{"Empty BAM file: no BGZF blocks"};
    }

    const std::size_t hdrSize{
        ComputeHeaderSize(std::span<const std::byte>{std::data(headerBuf), headerLen})};
    if (hdrSize > 0) {
        auto headerResult{SamHeader::FromBamHeaderBlock(
            std::span<const std::byte>{std::data(headerBuf), hdrSize})};
        if (!headerResult) {
            throw std::runtime_error{std::move(headerResult.error())};
        }
        syncHeader = std::move(*headerResult);
        recordBuf.resize(MAX_DECOMPRESSED_BLOCK_SIZE);
        recordBufSize = 0;
        recordBufPos = 0;
        syncEof = true;
        return;
    }

    throw std::runtime_error{"Incomplete BAM header"};
}

bool BgzfReader::Impl::RefillRecordBuffer()
{
    const std::size_t remaining{recordBufSize - recordBufPos};
    if ((recordBufPos > 0) && ((remaining == 0) || (recordBufPos > recordBufSize / 2))) {
        if (remaining > 0) {
            std::memmove(std::data(recordBuf), std::data(recordBuf) + recordBufPos, remaining);
        }
        recordBufPos = 0;
        recordBufSize = remaining;
    }

    if (syncEof) {
        return remaining > 0;
    }

    const std::size_t needed{recordBufSize + MAX_DECOMPRESSED_BLOCK_SIZE};
    if (std::size(recordBuf) < needed) {
        recordBuf.resize(needed);
    }

    const std::optional<std::size_t> bytesRead{ReadBlockSync(
        std::span<std::byte>{recordBuf}.subspan(recordBufSize, MAX_DECOMPRESSED_BLOCK_SIZE))};

    if (!bytesRead) {
        syncEof = true;
        return remaining > 0;
    }
    if (*bytesRead == 0) {
        syncEof = true;
        return remaining > 0;
    }

    recordBufSize += *bytesRead;
    return true;
}

std::optional<RawRecord> BgzfReader::Impl::ReadRecordSync()
{
    while ((recordBufSize - recordBufPos) < RECORD_BLOCK_SIZE_FIELD) {
        if (!RefillRecordBuffer()) {
            return std::nullopt;
        }
    }

    const std::byte* pos{std::data(recordBuf) + recordBufPos};
    const std::uint32_t blockSize{ReadU32LE(pos)};

    if (blockSize == 0) {
        return std::nullopt;
    }

    const std::size_t totalRecordBytes{RECORD_BLOCK_SIZE_FIELD + blockSize};

    while ((recordBufSize - recordBufPos) < totalRecordBytes) {
        if (!RefillRecordBuffer()) {
            return std::nullopt;
        }
    }

    pos = std::data(recordBuf) + recordBufPos + RECORD_BLOCK_SIZE_FIELD;
    recordBufPos += totalRecordBytes;

    return RawRecord{std::span<const std::byte>{pos, blockSize}};
}

std::optional<std::size_t> BgzfReader::Impl::ReadBlock(std::span<std::byte> buffer)
{
    if (mode == BgzfReaderMode::PARALLEL_BAM) {
        return pipeline->ReadBlockSync(buffer);
    }
    return ReadBlockSync(buffer);
}

void BgzfReader::Impl::Seek(VirtualOffset offset)
{
    if (mode == BgzfReaderMode::PARALLEL_BAM) {
        pipeline->Seek(offset);
        return;
    }

    syncFile.clear();
    syncFile.seekg(static_cast<std::streamoff>(offset.BlockOffset()), std::ios::beg);
    syncBlockOffset = offset.BlockOffset();

    if (mode != BgzfReaderMode::SYNC_BAM) {
        return;
    }

    recordBufPos = 0;
    recordBufSize = 0;
    syncEof = false;

    const std::uint16_t withinBlock{offset.WithinBlockOffset()};
    if (withinBlock == 0) {
        return;
    }

    if (std::size(recordBuf) < MAX_DECOMPRESSED_BLOCK_SIZE) {
        recordBuf.resize(MAX_DECOMPRESSED_BLOCK_SIZE);
    }

    const std::optional<std::size_t> bytesRead{
        ReadBlockSync(std::span<std::byte>{recordBuf}.first(MAX_DECOMPRESSED_BLOCK_SIZE))};
    if (!bytesRead || (*bytesRead == 0)) {
        syncEof = true;
        return;
    }

    recordBufSize = *bytesRead;
    if (recordBufSize >= withinBlock) {
        recordBufPos = withinBlock;
    }
}

VirtualOffset BgzfReader::Impl::Tell() const
{
    if (mode == BgzfReaderMode::PARALLEL_BAM) {
        return pipeline->Tell();
    }
    return VirtualOffset{syncBlockOffset, std::uint16_t{0}};
}

bool BgzfReader::Impl::HasEofMarker() const
{
    if (mode == BgzfReaderMode::PARALLEL_BAM) {
        return pipeline->hasEofMarker;
    }
    return hasEofMarker;
}

const SamHeader& BgzfReader::Impl::Header() const
{
    if (mode == BgzfReaderMode::BLOCK_ONLY) {
        throw std::runtime_error{"Header() is unavailable in block-only mode"};
    }
    if (mode == BgzfReaderMode::PARALLEL_BAM) {
        return pipeline->Header();
    }
    return syncHeader;
}

std::optional<RawRecord> BgzfReader::Impl::ReadRecord()
{
    if (mode == BgzfReaderMode::BLOCK_ONLY) {
        throw std::runtime_error{"ReadRecord() is unavailable in block-only mode"};
    }

    if (mode == BgzfReaderMode::PARALLEL_BAM) {
        return pipeline->ReadRecord();
    }

    auto rec{ReadRecordSync()};
    if (rec) {
        recordsConsumed.fetch_add(1, std::memory_order_relaxed);
    }
    return rec;
}

BgzfMetrics BgzfReader::Impl::GetMetrics() const
{
    if (mode == BgzfReaderMode::PARALLEL_BAM) {
        return pipeline->GetMetrics();
    }

    BgzfMetrics m{};
    m.RecordsConsumed = recordsConsumed.load(std::memory_order_relaxed);
    return m;
}

}  // namespace Samoa
}  // namespace PacBio
