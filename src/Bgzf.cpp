#include <pbsamoa/core/Bgzf.hpp>

#include "BinaryUtils.hpp"

#include <pbcopper/parallel/ThreadPool.h>
#include <pbcopper/third-party/rigtorp/SPSCQueue.hpp>

#include <libdeflate.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <fstream>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {
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
constexpr std::size_t MAX_DECOMPRESSED_BLOCK_SIZE{65536U};
constexpr std::size_t OUTPUT_QUEUE_CAPACITY{4096};
constexpr std::size_t RECORD_BLOCK_SIZE_FIELD{4};
constexpr std::size_t BLOCKS_PER_BATCH{32};

// --- Shared helpers ---

struct LibdeflateDecompressorDeleter
{
    void operator()(libdeflate_decompressor* decompressor) const
    {
        libdeflate_free_decompressor(decompressor);
    }
};

using DecompressorPtr = std::unique_ptr<libdeflate_decompressor, LibdeflateDecompressorDeleter>;

struct LibdeflateCompressorDeleter
{
    void operator()(libdeflate_compressor* compressor) const
    {
        libdeflate_free_compressor(compressor);
    }
};

using CompressorPtr = std::unique_ptr<libdeflate_compressor, LibdeflateCompressorDeleter>;

// --- Pipeline helpers ---

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

    const DecompressorPtr decompressor{libdeflate_alloc_decompressor()};
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

// =============================================================================
// BgzfReader
// =============================================================================

struct BgzfReader::Impl
{
    std::ifstream file{};
    std::vector<std::byte> compressedBuf{};
    std::uint64_t blockFileOffset{0};
    bool hasEofMarker{false};
    DecompressorPtr decompressor{};

    explicit Impl(const std::filesystem::path& path)
    {
        file.open(path, std::ios::binary);
        if (!file.is_open()) {
            throw std::runtime_error{"Cannot open file: " + path.string()};
        }

        compressedBuf.resize(MAX_BGZF_BLOCK_SIZE);

        decompressor = DecompressorPtr{libdeflate_alloc_decompressor()};
        if (!decompressor) {
            throw std::runtime_error{"Failed to create libdeflate decompressor"};
        }

        file.seekg(-static_cast<std::streamoff>(std::size(BGZF_EOF_MARKER)), std::ios::end);
        if (file.good()) {
            std::array<std::byte, 28> tail{};
            file.read(reinterpret_cast<char*>(std::data(tail)),
                      static_cast<std::streamsize>(std::size(tail)));
            hasEofMarker = IsBgzfEofMarker(tail);
        }

        file.clear();
        file.seekg(0, std::ios::beg);
    }

    ~Impl() = default;
    Impl(const Impl&) = delete;
    Impl& operator=(const Impl&) = delete;

    std::optional<std::size_t> ReadBlock(std::span<std::byte> buffer)
    {
        blockFileOffset = file.tellg();

        file.read(reinterpret_cast<char*>(std::data(compressedBuf)), BGZF_HEADER_SIZE);
        if (file.gcount() == 0) {
            return 0U;
        }
        if (file.gcount() < BGZF_HEADER_SIZE) {
            return std::nullopt;
        }

        const std::optional<BgzfBlockInfo> info{
            ParseBgzfBlockHeader(std::span<const std::byte>{compressedBuf}.first(
                static_cast<std::size_t>(BGZF_HEADER_SIZE)))};
        if (!info.has_value()) {
            return std::nullopt;
        }

        const std::size_t blockSize{info->blockSize};
        const std::size_t headerSize{static_cast<std::size_t>(BGZF_HEADER_SIZE)};
        const std::size_t remaining{blockSize - headerSize};
        if (std::size(compressedBuf) < blockSize) {
            compressedBuf.resize(blockSize);
        }

        file.read(reinterpret_cast<char*>(std::data(compressedBuf) + headerSize),
                  static_cast<std::streamsize>(remaining));
        if (static_cast<std::size_t>(file.gcount()) < remaining) {
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
            decompressor.get(), std::data(compressedBuf) + info->compressedDataOffset,
            info->compressedDataSize, std::data(buffer), isize, &actualOut)};
        if (result != LIBDEFLATE_SUCCESS) {
            return std::nullopt;
        }

        return actualOut;
    }
};

BgzfReader::BgzfReader(const std::filesystem::path& path) : impl_{std::make_unique<Impl>(path)} {}

BgzfReader::~BgzfReader() = default;
BgzfReader::BgzfReader(BgzfReader&&) noexcept = default;
BgzfReader& BgzfReader::operator=(BgzfReader&&) noexcept = default;

std::optional<std::size_t> BgzfReader::ReadBlock(std::span<std::byte> buffer)
{
    return impl_->ReadBlock(buffer);
}

void BgzfReader::Seek(VirtualOffset offset)
{
    impl_->file.clear();
    impl_->file.seekg(static_cast<std::streamoff>(offset.BlockOffset()), std::ios::beg);
    impl_->blockFileOffset = offset.BlockOffset();
}

bool BgzfReader::HasEofMarker() const { return impl_->hasEofMarker; }

VirtualOffset BgzfReader::Tell() const
{
    return VirtualOffset{impl_->blockFileOffset, std::uint16_t{0}};
}

// =============================================================================
// BgzfWriter
// =============================================================================

struct BgzfWriter::Impl
{
    std::ofstream file{};
    std::vector<std::byte> uncompressedBuf{};
    std::vector<std::uint8_t> compressedBuf{};
    CompressorPtr compressor{};
    std::uint64_t compressedOffset_{0};
    bool closed{false};

    explicit Impl(const std::filesystem::path& path, int level)
        : compressor{libdeflate_alloc_compressor(level)}
    {
        file.open(path, std::ios::binary);
        if (!file.is_open()) {
            throw std::runtime_error{"Cannot open file for writing: " + path.string()};
        }

        if (!compressor) {
            throw std::runtime_error{"Failed to create libdeflate compressor"};
        }

        uncompressedBuf.reserve(MAX_UNCOMPRESSED_SIZE);
        compressedBuf.resize(BGZF_MAX_BLOCK_SIZE);
    }

    ~Impl() = default;
    Impl(const Impl&) = delete;
    Impl& operator=(const Impl&) = delete;
};

BgzfWriter::BgzfWriter(const std::filesystem::path& path, int compressionLevel)
    : impl_{std::make_unique<Impl>(path, compressionLevel)}
{
}

BgzfWriter::~BgzfWriter()
{
    if ((impl_ != nullptr) && (!impl_->closed)) {
        Close();
    }
}

BgzfWriter::BgzfWriter(BgzfWriter&&) noexcept = default;
BgzfWriter& BgzfWriter::operator=(BgzfWriter&&) noexcept = default;

void BgzfWriter::Write(std::span<const std::byte> data)
{
    while (!std::empty(data)) {
        const std::size_t space{MAX_UNCOMPRESSED_SIZE - std::size(impl_->uncompressedBuf)};
        const std::size_t toCopy{std::ranges::min(space, std::size(data))};
        const std::span<const std::byte> chunk{data.first(toCopy)};

        impl_->uncompressedBuf.insert(std::ranges::end(impl_->uncompressedBuf),
                                      std::ranges::begin(chunk), std::ranges::end(chunk));
        data = data.subspan(toCopy);

        if (std::size(impl_->uncompressedBuf) >= MAX_UNCOMPRESSED_SIZE) {
            FlushBlock();
        }
    }
}

void BgzfWriter::Close()
{
    if (impl_->closed) {
        return;
    }

    if (!std::empty(impl_->uncompressedBuf)) {
        FlushBlock();
    }

    impl_->file.write(reinterpret_cast<const char*>(std::data(BGZF_EOF_MARKER)),
                      static_cast<std::streamsize>(std::size(BGZF_EOF_MARKER)));
    if (!impl_->file.good()) {
        throw std::runtime_error{"BgzfWriter: write failed during EOF marker"};
    }

    impl_->file.close();
    if (impl_->file.fail()) {
        throw std::runtime_error{"BgzfWriter: close failed"};
    }
    impl_->closed = true;
}

void BgzfWriter::FlushBlock()
{
    const std::vector<std::byte>& uncompressedBuffer{impl_->uncompressedBuf};
    std::vector<std::uint8_t>& compressedBuffer{impl_->compressedBuf};

    const std::size_t compressedSize{libdeflate_deflate_compress(
        impl_->compressor.get(), std::data(uncompressedBuffer), std::size(uncompressedBuffer),
        std::data(compressedBuffer), std::size(compressedBuffer))};
    if (compressedSize == 0U) {
        throw std::runtime_error{"BGZF compression failed"};
    }

    const std::uint32_t crc{
        libdeflate_crc32(0, std::data(uncompressedBuffer), std::size(uncompressedBuffer))};
    const std::uint32_t isize{static_cast<std::uint32_t>(std::size(uncompressedBuffer))};

    const std::uint16_t xlen{6U};
    const std::uint32_t blockSize{10U + 2U + xlen + static_cast<std::uint32_t>(compressedSize) +
                                  8U};
    const std::uint16_t bsize{static_cast<std::uint16_t>(blockSize - 1U)};

    std::ofstream& file{impl_->file};

    // Combined: BGZF header (12B) + BC extra field (6B) = 18B — one syscall
    const std::array<std::uint8_t, 18> frameHeader{
        31U,
        139U,
        8U,
        4U,  // ID1, ID2, CM, FLG
        0U,
        0U,
        0U,
        0U,  // MTIME
        0U,
        0U,  // XFL, OS
        static_cast<std::uint8_t>(xlen & 0xFFU),
        static_cast<std::uint8_t>(xlen >> 8U),
        66U,
        67U,
        2U,
        0U,  // BC subfield ID + length
        static_cast<std::uint8_t>(bsize & 0xFFU),
        static_cast<std::uint8_t>(bsize >> 8U),
    };
    file.write(reinterpret_cast<const char*>(std::data(frameHeader)),
               static_cast<std::streamsize>(std::size(frameHeader)));

    file.write(reinterpret_cast<const char*>(std::data(compressedBuffer)),
               static_cast<std::streamsize>(compressedSize));

    // Combined: CRC32 (4B) + ISIZE (4B) = 8B — one syscall
    const std::array<std::uint8_t, 8> trailer{
        static_cast<std::uint8_t>(crc & 0xFFU),
        static_cast<std::uint8_t>((crc >> 8U) & 0xFFU),
        static_cast<std::uint8_t>((crc >> 16U) & 0xFFU),
        static_cast<std::uint8_t>((crc >> 24U) & 0xFFU),
        static_cast<std::uint8_t>(isize & 0xFFU),
        static_cast<std::uint8_t>((isize >> 8U) & 0xFFU),
        static_cast<std::uint8_t>((isize >> 16U) & 0xFFU),
        static_cast<std::uint8_t>((isize >> 24U) & 0xFFU),
    };
    file.write(reinterpret_cast<const char*>(std::data(trailer)),
               static_cast<std::streamsize>(std::size(trailer)));

    if (!file.good()) {
        throw std::runtime_error{"BgzfWriter: write failed during BGZF block flush"};
    }

    impl_->uncompressedBuf.clear();
    impl_->compressedOffset_ = impl_->file.tellp();
}

VirtualOffset BgzfWriter::Tell() const
{
    return VirtualOffset{impl_->compressedOffset_,
                         static_cast<std::uint16_t>(std::size(impl_->uncompressedBuf))};
}

// =============================================================================
// BgzfPipeline
// =============================================================================

/// \brief Always-on atomic counters for pipeline introspection.
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

struct BgzfPipeline::Impl
{
    std::filesystem::path filePath;
    std::size_t numWorkers;
    bool hasEofMarker{false};

    // Sync decompression (for ReadBlock and header parsing)
    std::ifstream syncFile;
    DecompressorPtr syncDecompressor;
    std::array<std::byte, MAX_COMPRESSED_BLOCK_SIZE> syncCompressed;
    std::uint64_t syncBlockOffset{0};

    // Header
    SamHeader header;
    std::vector<std::byte> initialBytes;
    bool headerParsed{false};
    std::uint64_t headerEndFileOffset{0};

    // Pipeline (created after ParseHeader when numWorkers > 0)
    std::unique_ptr<Parallel::ThreadPool<DecompressedBatch>> pool;
    std::jthread ioThread;
    std::jthread consumerThread;
    std::unique_ptr<rigtorp::SPSCQueue<RawRecord>> outputQueue;
    std::mutex readyMutex;
    std::condition_variable readyCv;
    std::atomic<bool> done{false};
    std::atomic<bool> pipelineError{false};
    std::exception_ptr errorPtr;

    // Per-worker decompressors (indexed by ThreadIndex)
    std::vector<DecompressorPtr> decompressors;

    // Always-on metrics
    PipelineCounters counters;

    explicit Impl(const std::filesystem::path& path, std::size_t workers);
    ~Impl();

    Impl(const Impl&) = delete;
    Impl& operator=(const Impl&) = delete;

    void ParseHeader();
    void StartPipeline();
    void StopPipeline();
    void IoLoop(std::stop_token stopToken);
    void ConsumerLoop(std::stop_token stopToken);
    std::optional<std::size_t> ReadBlockSync(std::span<std::byte> buffer);
};

BgzfPipeline::Impl::Impl(const std::filesystem::path& path, std::size_t workers)
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
    syncDecompressor = DecompressorPtr{libdeflate_alloc_decompressor()};
    if (!syncDecompressor) {
        throw std::runtime_error{"Failed to create libdeflate decompressor"};
    }
}

BgzfPipeline::Impl::~Impl() { StopPipeline(); }

void BgzfPipeline::Impl::ParseHeader()
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
        if (!info.has_value()) {
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
        const std::size_t hdrSize{ComputeHeaderSize(std::data(headerBuf), headerLen)};
        if (hdrSize > 0) {
            header = SamHeader::FromBamHeaderBlock(
                std::span<const std::byte>{std::data(headerBuf), hdrSize});

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
    const std::size_t hdrSize{ComputeHeaderSize(std::data(headerBuf), headerLen)};
    if (hdrSize > 0) {
        header = SamHeader::FromBamHeaderBlock(
            std::span<const std::byte>{std::data(headerBuf), hdrSize});
        headerParsed = true;
    } else {
        throw std::runtime_error{"Incomplete BAM header"};
    }
}

void BgzfPipeline::Impl::StartPipeline()
{
    // Create per-worker decompressors
    decompressors.clear();
    decompressors.reserve(numWorkers);
    for (std::size_t i{0}; i < numWorkers; ++i) {
        DecompressorPtr d{libdeflate_alloc_decompressor()};
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
    pool = std::make_unique<Parallel::ThreadPool<DecompressedBatch>>(
        Parallel::ThreadPool<DecompressedBatch>::Config{
            .NumThreads = numWorkers,
            .QueueMultiplier = 3,
            .EnableMetrics = true,
        });

    // Start consumer thread first — must enter ConsumeWith() before IO thread calls Submit()
    consumerThread = std::jthread{[this](std::stop_token st) { ConsumerLoop(st); }};

    // Start IO thread
    ioThread = std::jthread{[this](std::stop_token st) { IoLoop(st); }};
}

void BgzfPipeline::Impl::StopPipeline()
{
    // 1. Signal IO thread to stop submitting
    if (ioThread.joinable()) {
        ioThread.request_stop();
    }

    // 2. Signal consumer to stop blocking on SPSC push
    if (consumerThread.joinable()) {
        consumerThread.request_stop();
    }

    // 3. Finalize pool — stops workers, unblocks Submit(), makes ConsumeWith() return false
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

void BgzfPipeline::Impl::IoLoop(std::stop_token stopToken)
{
    try {
        std::ifstream file{filePath, std::ios::binary};
        if (!file.is_open()) {
            throw std::runtime_error{"IO thread: cannot open file: " + filePath.string()};
        }
        file.seekg(static_cast<std::streamoff>(headerEndFileOffset), std::ios::beg);

        std::array<std::byte, MAX_COMPRESSED_BLOCK_SIZE> compressed{};
        std::vector<std::byte> batchBuffer;
        std::vector<CompressedEntry> batchEntries;
        batchEntries.reserve(BLOCKS_PER_BATCH);

        const std::size_t poolCapacity{numWorkers * 3};

        while (!stopToken.stop_requested()) {
            // Accumulate up to BLOCKS_PER_BATCH compressed payloads
            batchBuffer.clear();
            batchEntries.clear();

            std::uint64_t batchCompressedBytes{0};

            for (std::size_t b{0}; b < BLOCKS_PER_BATCH && !stopToken.stop_requested(); ++b) {
                const std::uint64_t blockOffset{static_cast<std::uint64_t>(file.tellg())};

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
                if (!info.has_value()) {
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
                batchEntries.push_back(CompressedEntry{
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
            pool->Submit([entries = batchEntries, cdata = std::move(batchBuffer),
                          &decomps = this->decompressors, &ctr = this->counters](
                             Parallel::ThreadIndex threadIdx) -> DecompressedBatch {
                const auto decompStart{std::chrono::steady_clock::now()};

                // Compute total expected output size
                std::uint64_t totalIsize{0};
                for (const auto& entry : entries) {
                    totalIsize += entry.isize;
                }

                DecompressedBatch batch;
                batch.data.resize(totalIsize);

                std::size_t outOffset{0};
                for (const auto& entry : entries) {
                    std::size_t actualOut{0};
                    const libdeflate_result result{libdeflate_deflate_decompress(
                        decomps[threadIdx.Value()].get(), std::data(cdata) + entry.offset,
                        entry.size, std::data(batch.data) + outOffset, entry.isize, &actualOut)};

                    if (result != LIBDEFLATE_SUCCESS) {
                        throw std::runtime_error{"BGZF decompression failed"};
                    }
                    outOffset += actualOut;
                }
                batch.data.resize(outOffset);

                // Update decompression counters
                ctr.decompressNs.fetch_add(
                    static_cast<std::uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                                                   std::chrono::steady_clock::now() - decompStart)
                                                   .count()),
                    std::memory_order_relaxed);
                ctr.bytesDecompressed.fetch_add(outOffset, std::memory_order_relaxed);

                return batch;
            });
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

void BgzfPipeline::Impl::ConsumerLoop(std::stop_token stopToken)
{
    // Signal the reader via condition variable.  Lock/unlock the mutex to
    // establish happens-before with the reader's cv.wait() predicate.
    const auto notifyReader = [this]() {
        {
            const std::lock_guard lock{readyMutex};
        }
        readyCv.notify_one();
    };

    try {
        std::vector<std::byte> accumulator{std::move(initialBytes)};
        std::size_t pos{0};

        // Parse complete records from accumulator and push to SPSC queue.
        // Also compacts consumed bytes from the front.  Notifies once per
        // call (batch-level) rather than per record.
        const auto parseRecords = [&]() {
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
                while (!outputQueue->try_push(std::move(view))) {
                    if (stopToken.stop_requested()) {
                        counters.recordParseNs.fetch_add(
                            static_cast<std::uint64_t>(
                                std::chrono::duration_cast<std::chrono::nanoseconds>(
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
                notifyReader();
            }
        };

        // Process any records already present from initialBytes (header block
        // may contain record data when header + records fit in one BGZF block).
        parseRecords();

        while (pool->ConsumeWith([&](DecompressedBatch batch) {
            // Append decompressed data to accumulator
            accumulator.insert(std::end(accumulator), std::data(batch.data),
                               std::data(batch.data) + std::size(batch.data));
            parseRecords();
        })) {
            if (stopToken.stop_requested()) {
                break;
            }
        }
    } catch (...) {
        errorPtr = std::current_exception();
        pipelineError.store(true, std::memory_order_release);
    }

    done.store(true, std::memory_order_release);
    notifyReader();
}

std::optional<std::size_t> BgzfPipeline::Impl::ReadBlockSync(std::span<std::byte> buffer)
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
    if (!info.has_value()) {
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

// --- BgzfPipeline public API ---

BgzfPipeline::BgzfPipeline(const std::filesystem::path& path, std::size_t numWorkers)
    : impl_{std::make_unique<Impl>(path, numWorkers)}
{
}

BgzfPipeline::~BgzfPipeline() = default;

void BgzfPipeline::ParseHeader() { impl_->ParseHeader(); }

const SamHeader& BgzfPipeline::Header() const
{
    if (!impl_->headerParsed) {
        throw std::runtime_error{"ParseHeader() must be called before accessing Header()"};
    }
    return impl_->header;
}

std::optional<RawRecord> BgzfPipeline::ReadRecord()
{
    if (!impl_->headerParsed) {
        throw std::runtime_error{"ParseHeader() must be called before ReadRecord()"};
    }
    if (!impl_->outputQueue) {
        throw std::runtime_error{"ReadRecord() requires numWorkers > 0"};
    }

    // Fast path: check SPSC queue without locking.
    while (true) {
        RawRecord* front{impl_->outputQueue->front()};
        if (front != nullptr) {
            RawRecord result{std::move(*front)};
            impl_->outputQueue->pop();
            impl_->counters.recordsConsumed.fetch_add(1, std::memory_order_relaxed);
            return result;
        }

        // Queue empty — check termination before sleeping.
        if (impl_->pipelineError.load(std::memory_order_acquire)) {
            front = impl_->outputQueue->front();
            if (front != nullptr) {
                RawRecord result{std::move(*front)};
                impl_->outputQueue->pop();
                impl_->counters.recordsConsumed.fetch_add(1, std::memory_order_relaxed);
                return result;
            }
            std::rethrow_exception(impl_->errorPtr);
        }

        if (impl_->done.load(std::memory_order_acquire)) {
            front = impl_->outputQueue->front();
            if (front != nullptr) {
                RawRecord result{std::move(*front)};
                impl_->outputQueue->pop();
                impl_->counters.recordsConsumed.fetch_add(1, std::memory_order_relaxed);
                return result;
            }
            return std::nullopt;
        }

        // Slow path: wait on condition variable for data/done/error.
        // The mutex lock creates a happens-before with the consumer's
        // lock/unlock in notifyReader(), preventing lost wakeups.
        impl_->counters.readerStalls.fetch_add(1, std::memory_order_relaxed);
        std::unique_lock lock{impl_->readyMutex};
        impl_->readyCv.wait(lock, [this]() {
            return impl_->outputQueue->front() != nullptr ||
                   impl_->done.load(std::memory_order_acquire) ||
                   impl_->pipelineError.load(std::memory_order_acquire);
        });
    }
}

std::optional<std::size_t> BgzfPipeline::ReadBlock(std::span<std::byte> buffer)
{
    return impl_->ReadBlockSync(buffer);
}

void BgzfPipeline::Seek(VirtualOffset offset)
{
    // Stop pipeline if running
    impl_->StopPipeline();

    // Reset sync file position
    impl_->syncFile.clear();
    impl_->syncFile.seekg(static_cast<std::streamoff>(offset.BlockOffset()), std::ios::beg);
    impl_->syncBlockOffset = offset.BlockOffset();

    // Update pipeline start offset
    impl_->headerEndFileOffset = offset.BlockOffset();
    impl_->initialBytes.clear();

    // Reset pipeline state
    impl_->done.store(false, std::memory_order_relaxed);
    impl_->pipelineError.store(false, std::memory_order_relaxed);
    impl_->errorPtr = nullptr;

    // Restart pipeline if in parallel mode and header was parsed
    if ((impl_->numWorkers > 0) && impl_->headerParsed) {
        impl_->StartPipeline();
    }
}

VirtualOffset BgzfPipeline::Tell() const { return VirtualOffset{impl_->syncBlockOffset, 0}; }

bool BgzfPipeline::HasEofMarker() const { return impl_->hasEofMarker; }

BgzfMetrics BgzfPipeline::GetMetrics() const
{
    BgzfMetrics m{};

    // Throughput counters
    m.BytesRead = impl_->counters.bytesRead.load(std::memory_order_relaxed);
    m.BlocksRead = impl_->counters.blocksRead.load(std::memory_order_relaxed);
    m.BytesDecompressed = impl_->counters.bytesDecompressed.load(std::memory_order_relaxed);

    // ThreadPool metrics
    if (impl_->pool) {
        const auto poolSnap{impl_->pool->GetMetrics()};
        m.PoolQueueDepth = poolSnap.CurrentQueueDepth;
        m.PoolPeakQueueDepth = poolSnap.PeakQueueDepth;
        m.PoolActiveWorkers = poolSnap.CurrentActiveTasks;
        m.PoolPeakActiveWorkers = poolSnap.PeakActiveTasks;
        m.PoolResultQueueDepth = poolSnap.CurrentResultQueueDepth;
        m.PoolPeakResultQueueDepth = poolSnap.PeakResultQueueDepth;
    }

    // SPSC queue
    m.RecordsProduced = impl_->counters.recordsProduced.load(std::memory_order_relaxed);
    m.RecordsConsumed = impl_->counters.recordsConsumed.load(std::memory_order_relaxed);

    // Stall counters
    m.IoStalls = impl_->counters.ioStalls.load(std::memory_order_relaxed);
    m.ConsumerStalls = impl_->counters.consumerStalls.load(std::memory_order_relaxed);
    m.ReaderStalls = impl_->counters.readerStalls.load(std::memory_order_relaxed);

    // Timing
    m.IoReadNs = impl_->counters.ioReadNs.load(std::memory_order_relaxed);
    m.DecompressNs = impl_->counters.decompressNs.load(std::memory_order_relaxed);
    m.RecordParseNs = impl_->counters.recordParseNs.load(std::memory_order_relaxed);

    return m;
}

}  // namespace Samoa
}  // namespace PacBio
