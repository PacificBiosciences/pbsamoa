#include <pbsamoa/index/BaiIndex.hpp>

#include "BinaryUtils.hpp"
#include "LibdeflateUtils.hpp"

#include <pbsamoa/core/Bgzf.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/io/BamRawReader.hpp>

#include <pbcopper/parallel/FireAndForgetIndexed.h>

#include <libdeflate.h>

#include <algorithm>
#include <array>
#include <filesystem>
#include <format>
#include <fstream>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <cassert>
#include <cstdint>
#include <cstring>

namespace PacBio {
namespace Samoa {

namespace {

constexpr std::size_t BAI_LINEAR_INDEX_WINDOW{16384};  // 2^14 = 16 kbp
constexpr std::uint32_t BAI_METADATA_BIN{37450U};
constexpr std::array<std::byte, 4> BAI_MAGIC{
    std::byte{'B'},
    std::byte{'A'},
    std::byte{'I'},
    std::byte{'\1'},
};
constexpr std::size_t BAM_FIXED_FIELDS_SIZE{32};
constexpr std::size_t BGZF_BLOCK_HEADER_SIZE{18U};
constexpr std::size_t MAX_COMPRESSED_BLOCK_SIZE{65536U};  // BSIZE is 16-bit: blockSize <= 65536
// Compressed blocks inflated per worker before the consumer drains the window.
// Bounds memory to numWorkers * WINDOW_BLOCKS_PER_WORKER * (64 KiB + 64 KiB).
constexpr std::size_t WINDOW_BLOCKS_PER_WORKER{4U};

struct BlockSegment
{
    std::uint64_t blockOffset;
    std::uint16_t withinBlockStart;
    std::size_t byteCount;
};

void ReadExact(std::ifstream& in, std::span<std::byte> bytes, std::string_view fieldName)
{
    in.read(reinterpret_cast<char*>(std::data(bytes)),
            static_cast<std::streamsize>(std::size(bytes)));
    if (!in) {
        throw std::runtime_error{std::format("Truncated BAI file reading {}", fieldName)};
    }
}

template <typename T>
    requires std::is_arithmetic_v<T>
T ReadLEFromFile(std::ifstream& in, std::string_view fieldName)
{
    std::array<std::byte, sizeof(T)> bytes{};
    ReadExact(in, bytes, fieldName);
    T value{};
    std::memcpy(&value, std::data(bytes), sizeof(T));
    return value;
}

std::size_t CheckedNonNegativeCount(std::int32_t value, std::string_view fieldName)
{
    if (value < 0) {
        throw std::runtime_error{
            std::format("Invalid BAI file: {} must be >= 0, got {}", fieldName, value)};
    }
    return static_cast<std::size_t>(value);
}

template <typename T>
    requires std::is_arithmetic_v<T>
void WriteLEToFile(std::ofstream& out, T v)
{
    std::array<std::byte, sizeof(T)> bytes{};
    std::memcpy(std::data(bytes), &v, sizeof(T));
    out.write(reinterpret_cast<const char*>(std::data(bytes)),
              static_cast<std::streamsize>(std::size(bytes)));
}

template <typename T>
    requires std::is_integral_v<T>
std::int32_t CheckedInt32(T value, std::string_view fieldName)
{
    using Int32Limits = std::numeric_limits<std::int32_t>;
    if constexpr (std::is_signed_v<T>) {
        if ((value < Int32Limits::min()) || (value > Int32Limits::max())) {
            throw std::runtime_error{std::format("Invalid BAI value for {}: {}", fieldName, value)};
        }
    } else if (value > static_cast<std::uint64_t>(Int32Limits::max())) {
        throw std::runtime_error{std::format("Invalid BAI value for {}: {}", fieldName, value)};
    }
    return static_cast<std::int32_t>(value);
}

/// \brief Merge overlapping/adjacent chunks sorted by Begin.
std::vector<Chunk> MergeChunks(std::vector<Chunk> chunks)
{
    if (std::empty(chunks)) {
        return {};
    }

    std::ranges::sort(chunks, {}, &Chunk::Begin);

    std::vector<Chunk> merged;
    merged.reserve(std::size(chunks));
    merged.push_back(chunks[0]);

    for (std::size_t i{1}; i < std::size(chunks); ++i) {
        const Chunk& chunk{chunks[i]};
        Chunk& last{merged.back()};
        if (chunk.Begin <= last.End) {
            if (chunk.End > last.End) {
                last.End = chunk.End;
            }
        } else {
            merged.push_back(chunk);
        }
    }

    return merged;
}

bool ConsumesReference(std::uint8_t opCode)
{
    return (opCode == 0) || (opCode == 2) || (opCode == 3) || (opCode == 7) || (opCode == 8);
}

VirtualOffset VirtualOffsetAt(std::span<const BlockSegment> segments, std::size_t accumPos)
{
    std::size_t pos{0};
    for (const BlockSegment& segment : segments) {
        if (accumPos < (pos + segment.byteCount)) {
            const std::size_t within{segment.withinBlockStart + (accumPos - pos)};
            assert(within <= std::numeric_limits<std::uint16_t>::max());
            return VirtualOffset{segment.blockOffset, static_cast<std::uint16_t>(within)};
        }
        pos += segment.byteCount;
    }

    if (!std::empty(segments)) {
        const BlockSegment& last{segments.back()};
        const std::size_t pastEnd{last.withinBlockStart + last.byteCount};
        assert(pastEnd <= std::numeric_limits<std::uint16_t>::max());
        return VirtualOffset{last.blockOffset, static_cast<std::uint16_t>(pastEnd)};
    }

    return VirtualOffset{};
}

void CompactAccumulator(std::vector<std::byte>& recordAccum, std::vector<BlockSegment>& segments,
                        std::size_t& accumPos)
{
    std::size_t consumed{0};
    std::size_t segmentIndex{0};
    while ((segmentIndex < std::size(segments)) &&
           ((consumed + segments[segmentIndex].byteCount) <= accumPos)) {
        consumed += segments[segmentIndex].byteCount;
        ++segmentIndex;
    }

    if (segmentIndex > 0) {
        if ((consumed < accumPos) && (segmentIndex < std::size(segments))) {
            const std::size_t skipInSegment{accumPos - consumed};
            const std::size_t newWithin{segments[segmentIndex].withinBlockStart + skipInSegment};
            assert(newWithin <= std::numeric_limits<std::uint16_t>::max());
            segments[segmentIndex].withinBlockStart = static_cast<std::uint16_t>(newWithin);
            segments[segmentIndex].byteCount -= skipInSegment;
        }
        segments.erase(std::ranges::begin(segments),
                       std::ranges::begin(segments) + static_cast<std::ptrdiff_t>(segmentIndex));
    }

    recordAccum.erase(std::ranges::begin(recordAccum),
                      std::ranges::begin(recordAccum) + static_cast<std::ptrdiff_t>(accumPos));
    accumPos = 0;
}

/// \brief One decompressed BGZF block tagged with the compressed file offset of
///        its start. The data span is valid only until the next Next() call.
struct DecodedBlock
{
    std::uint64_t coffset{0};
    std::span<const std::byte> data{};
};

/// \brief Ordered source of decompressed BGZF blocks. Build() consumes blocks
///        through this seam, so the (coffset, bytes) contract — and therefore the
///        virtual-offset arithmetic in VirtualOffsetAt() — is identical for the
///        serial and parallel implementations; only the decompression differs.
class BlockSource
{
public:
    BlockSource() = default;
    BlockSource(const BlockSource&) = delete;
    BlockSource& operator=(const BlockSource&) = delete;
    virtual ~BlockSource() = default;

    /// \brief Next non-empty block in file order, or nullopt at clean EOF.
    /// \throws std::runtime_error on a truncated or corrupt block.
    virtual std::optional<DecodedBlock> Next() = 0;
};

/// \brief Single-threaded source wrapping BgzfReader's synchronous ReadBlock().
///        Tell() reports the start offset of the block just read, matching the
///        coffset contract; this preserves the original Build() behaviour exactly.
class SerialBlockSource final : public BlockSource
{
public:
    explicit SerialBlockSource(const std::filesystem::path& path)
        : reader_{path}, buffer_(MAX_DECOMPRESSED_BLOCK_SIZE)
    {
    }

    std::optional<DecodedBlock> Next() override
    {
        const std::optional<std::size_t> bytesRead{
            reader_.ReadBlock(std::span<std::byte>{buffer_})};
        if (!bytesRead || (*bytesRead == 0)) {
            return std::nullopt;
        }
        return DecodedBlock{reader_.Tell().BlockOffset(),
                            std::span<const std::byte>{buffer_}.first(*bytesRead)};
    }

private:
    BgzfReader reader_;
    std::vector<std::byte> buffer_;
};

/// \brief Multi-threaded source: frames compressed BGZF blocks sequentially (cheap
///        I/O), then inflates a bounded window of them in parallel. Blocks are
///        emitted in file order with their coffsets, so the resulting index is
///        byte-identical to the serial path for any worker count (deterministic).
class ParallelBlockSource final : public BlockSource
{
    struct Slot
    {
        std::uint64_t coffset{0};
        std::vector<std::byte> compressed;
        BgzfBlockInfo info{};
        std::uint32_t isize{0};
        std::uint32_t expectedCrc{0};
        std::vector<std::byte> out;
        std::size_t outSize{0};
    };

public:
    ParallelBlockSource(const std::filesystem::path& path, std::size_t numWorkers)
        : numWorkers_{std::max<std::size_t>(numWorkers, 1)}
        , file_{path, std::ios::binary}
        , slots_(numWorkers_ * WINDOW_BLOCKS_PER_WORKER)
    {
        if (!file_) {
            throw std::runtime_error{"Cannot open BAM file: " + path.string()};
        }
        decompressors_.reserve(numWorkers_);
        for (std::size_t i{0}; i < numWorkers_; ++i) {
            LibdeflateDecompressorPtr decompressor{libdeflate_alloc_decompressor()};
            if (!decompressor) {
                throw std::runtime_error{"Failed to allocate libdeflate decompressor"};
            }
            decompressors_.push_back(std::move(decompressor));
        }
        for (Slot& slot : slots_) {
            slot.compressed.resize(MAX_COMPRESSED_BLOCK_SIZE);
            slot.out.resize(MAX_DECOMPRESSED_BLOCK_SIZE);
        }
    }

    std::optional<DecodedBlock> Next() override
    {
        if (cursor_ >= ready_) {
            FillWindow();
            if (ready_ == 0) {
                return std::nullopt;
            }
        }
        const Slot& slot{slots_[cursor_]};
        ++cursor_;
        return DecodedBlock{slot.coffset, std::span<const std::byte>{slot.out}.first(slot.outSize)};
    }

private:
    // Read the next compressed block into slot. Returns false at physical EOF.
    // Empty blocks (ISIZE == 0: the EOF marker or appended-file boundaries) are
    // skipped here, mirroring BgzfReader::ReadBlockSync.
    bool ReadCompressedBlock(Slot& slot)
    {
        while (true) {
            const std::uint64_t coffset{static_cast<std::uint64_t>(file_.tellg())};
            file_.read(reinterpret_cast<char*>(std::data(slot.compressed)),
                       static_cast<std::streamsize>(BGZF_BLOCK_HEADER_SIZE));
            const std::streamsize headerRead{file_.gcount()};
            if (headerRead == 0) {
                return false;  // physical EOF
            }
            if (headerRead < static_cast<std::streamsize>(BGZF_BLOCK_HEADER_SIZE)) {
                throw std::runtime_error{"Truncated BGZF block header during BAI build"};
            }

            const std::optional<BgzfBlockInfo> info{ParseBgzfBlockHeader(
                std::span<const std::byte>{slot.compressed}.first(BGZF_BLOCK_HEADER_SIZE))};
            if (!info) {
                throw std::runtime_error{"Invalid BGZF block header during BAI build"};
            }

            const std::size_t blockSize{info->blockSize};
            if ((blockSize < BGZF_BLOCK_HEADER_SIZE) || (blockSize > MAX_COMPRESSED_BLOCK_SIZE)) {
                throw std::runtime_error{"Invalid BGZF block size during BAI build"};
            }

            const std::streamsize remaining{static_cast<std::streamsize>(blockSize) -
                                            static_cast<std::streamsize>(BGZF_BLOCK_HEADER_SIZE)};
            file_.read(reinterpret_cast<char*>(std::data(slot.compressed) + BGZF_BLOCK_HEADER_SIZE),
                       remaining);
            if (file_.gcount() < remaining) {
                throw std::runtime_error{"Truncated BGZF block during BAI build"};
            }

            const std::uint32_t isize{ReadU32LE(std::data(slot.compressed) + blockSize - 4U)};
            if (isize == 0U) {
                continue;  // empty block — skip
            }

            slot.coffset = coffset;
            slot.info = *info;
            slot.isize = isize;
            slot.expectedCrc = ReadU32LE(std::data(slot.compressed) + blockSize - 8U);
            return true;
        }
    }

    void InflateSlot(Parallel::FireAndForgetIndexed::Index workerIdx, std::size_t slotIdx)
    {
        Slot& slot{slots_[slotIdx]};
        std::size_t actualOut{0};
        const libdeflate_result result{libdeflate_deflate_decompress(
            decompressors_[static_cast<std::size_t>(workerIdx)].get(),
            std::data(slot.compressed) + slot.info.compressedDataOffset,
            slot.info.compressedDataSize, std::data(slot.out), slot.isize, &actualOut)};
        if (result != LIBDEFLATE_SUCCESS) {
            throw std::runtime_error{"BGZF decompression failed during BAI build"};
        }
        if (libdeflate_crc32(0, std::data(slot.out), actualOut) != slot.expectedCrc) {
            throw std::runtime_error{"BGZF CRC32 mismatch during BAI build"};
        }
        slot.outSize = actualOut;
    }

    void FillWindow()
    {
        cursor_ = 0;
        ready_ = 0;
        while (ready_ < std::size(slots_)) {
            if (!ReadCompressedBlock(slots_[ready_])) {
                break;
            }
            ++ready_;
        }
        if (ready_ == 0) {
            return;
        }

        // Inflate the framed blocks across numWorkers_ threads. Each worker owns a
        // private decompressor (indexed by its worker id) and writes a disjoint
        // slot, so there are no data races. Finalize() joins every task and
        // rethrows the first exception (corruption fails loud).
        Parallel::FireAndForgetIndexed pool{static_cast<std::int32_t>(numWorkers_)};
        for (std::size_t i{0}; i < ready_; ++i) {
            pool.ProduceWith([this](Parallel::FireAndForgetIndexed::Index workerIdx,
                                    std::size_t slotIdx) { InflateSlot(workerIdx, slotIdx); },
                             i);
        }
        pool.Finalize();
    }

    std::size_t numWorkers_;
    std::ifstream file_;
    std::vector<LibdeflateDecompressorPtr> decompressors_;
    std::vector<Slot> slots_;
    std::size_t ready_{0};
    std::size_t cursor_{0};
};

std::unique_ptr<BlockSource> MakeBlockSource(const std::filesystem::path& path,
                                             std::size_t numWorkers)
{
    if (numWorkers > 1) {
        return std::make_unique<ParallelBlockSource>(path, numWorkers);
    }
    return std::make_unique<SerialBlockSource>(path);
}

bool EnsureAccumulatedBytes(BlockSource& source, std::vector<std::byte>& recordAccum,
                            std::vector<BlockSegment>& segments, std::uint64_t& currentBlockOffset,
                            std::size_t accumPos, std::size_t requiredBytes, bool& eof)
{
    while ((std::size(recordAccum) - accumPos) < requiredBytes) {
        if (eof) {
            return false;
        }
        const std::optional<DecodedBlock> block{source.Next()};
        if (!block) {
            eof = true;
            return false;
        }
        currentBlockOffset = block->coffset;
        recordAccum.insert(std::ranges::end(recordAccum), std::ranges::begin(block->data),
                           std::ranges::end(block->data));
        segments.push_back(BlockSegment{currentBlockOffset, 0, std::size(block->data)});
    }
    return true;
}

}  // namespace

// --- BaiIndex accessors ---

std::int32_t BaiIndex::NumReferences() const { return std::ssize(references_); }

const ReferenceIndex& BaiIndex::Reference(std::int32_t refId) const
{
    return references_.at(refId);
}

std::uint64_t BaiIndex::MappedCount() const { return mappedCount_; }

std::uint64_t BaiIndex::UnmappedCount() const { return unmappedCount_; }

// --- FromFile ---

BaiIndex BaiIndex::FromFile(const std::filesystem::path& path)
{
    std::ifstream file{path, std::ios::binary};
    if (!file.is_open()) {
        throw std::runtime_error{"Cannot open BAI file: " + path.string()};
    }

    // Read magic
    std::array<std::byte, 4> magic{};
    ReadExact(file, magic, "magic");
    if (magic != BAI_MAGIC) {
        throw std::runtime_error{"Invalid BAI magic in: " + path.string()};
    }

    // Read n_ref
    const std::int32_t nRef{ReadLEFromFile<std::int32_t>(file, "n_ref")};
    const std::size_t nRefCount{CheckedNonNegativeCount(nRef, "n_ref")};

    BaiIndex index;
    index.references_.resize(nRefCount);

    for (std::size_t r{0}; r < nRefCount; ++r) {
        ReferenceIndex& ref{index.references_[r]};

        // n_bin
        const std::int32_t nBin{ReadLEFromFile<std::int32_t>(file, "n_bin")};
        const std::size_t nBinCount{CheckedNonNegativeCount(nBin, "n_bin")};

        for (std::size_t b{0}; b < nBinCount; ++b) {
            const std::uint32_t binNumber{ReadLEFromFile<std::uint32_t>(file, "bin number")};
            const std::int32_t nChunks{ReadLEFromFile<std::int32_t>(file, "n_chunks")};
            const std::size_t nChunkCount{CheckedNonNegativeCount(nChunks, "n_chunks")};

            std::vector<Chunk> chunks;
            chunks.reserve(nChunkCount);
            for (std::size_t c{0}; c < nChunkCount; ++c) {
                const std::uint64_t chunkBeg{ReadLEFromFile<std::uint64_t>(file, "chunk begin")};
                const std::uint64_t chunkEnd{ReadLEFromFile<std::uint64_t>(file, "chunk end")};
                chunks.push_back(Chunk{VirtualOffset{chunkBeg}, VirtualOffset{chunkEnd}});
            }

            if (binNumber == BAI_METADATA_BIN) {
                if (nChunkCount != 2) {
                    throw std::runtime_error{std::format(
                        "Invalid BAI metadata bin: expected 2 chunks, got {}", nChunkCount)};
                }
                index.mappedCount_ += chunks[1].Begin.Value();
                index.unmappedCount_ += chunks[1].End.Value();
            }

            ref.bins[binNumber] = std::move(chunks);
        }

        // n_intv
        const std::int32_t nIntv{ReadLEFromFile<std::int32_t>(file, "n_intv")};
        const std::size_t nIntvCount{CheckedNonNegativeCount(nIntv, "n_intv")};

        ref.linearIndex.resize(nIntvCount);
        for (std::size_t i{0}; i < nIntvCount; ++i) {
            const std::uint64_t offset{ReadLEFromFile<std::uint64_t>(file, "linear index offset")};
            ref.linearIndex[i] = VirtualOffset{offset};
        }
    }

    // Optional: read n_no_coor (unmapped count) at end of file
    {
        std::array<std::byte, sizeof(std::uint64_t)> bytes{};
        file.read(reinterpret_cast<char*>(std::data(bytes)),
                  static_cast<std::streamsize>(std::size(bytes)));
        const std::streamsize bytesRead{file.gcount()};
        if (bytesRead != 0) {
            if (bytesRead != static_cast<std::streamsize>(std::size(bytes))) {
                throw std::runtime_error{"Truncated BAI file reading n_no_coor"};
            }
            std::uint64_t value{};
            std::memcpy(&value, std::data(bytes), sizeof(value));
            index.unmappedCount_ += value;
            index.noCoorCount_ += value;
        }
    }

    return index;
}

// --- ToFile ---

void BaiIndex::ToFile(const std::filesystem::path& path) const
{
    std::ofstream file{path, std::ios::binary};
    if (!file.is_open()) {
        throw std::runtime_error{"Cannot create BAI file: " + path.string()};
    }
    file.exceptions(std::ios::failbit | std::ios::badbit);

    try {
        // Magic
        file.write("BAI\1", 4);

        // n_ref
        const std::int32_t nRef{CheckedInt32(std::ssize(references_), "n_ref")};
        WriteLEToFile(file, nRef);

        for (const ReferenceIndex& ref : references_) {
            // n_bin
            const std::int32_t nBin{CheckedInt32(std::ssize(ref.bins), "n_bin")};
            WriteLEToFile(file, nBin);

            for (const auto& [binNumber, chunks] : ref.bins) {
                WriteLEToFile(file, binNumber);
                const std::int32_t nChunks{CheckedInt32(std::ssize(chunks), "n_chunk")};
                WriteLEToFile(file, nChunks);
                for (const Chunk& chunk : chunks) {
                    WriteLEToFile(file, chunk.Begin.Value());
                    WriteLEToFile(file, chunk.End.Value());
                }
            }

            // n_intv
            const std::int32_t nIntv{CheckedInt32(std::ssize(ref.linearIndex), "n_intv")};
            WriteLEToFile(file, nIntv);
            for (const VirtualOffset& offset : ref.linearIndex) {
                WriteLEToFile(file, offset.Value());
            }
        }

        // n_no_coor: only unplaced (refId < 0) reads, per SAMv1 §5.2 and hts.c.
        WriteLEToFile(file, noCoorCount_);
    } catch (...) {
        file.close();
        std::filesystem::remove(path);
        throw;
    }
}

// --- Query ---

std::vector<Chunk> BaiIndex::Query(std::int32_t refId, std::int32_t beg, std::int32_t end) const
{
    if ((refId < 0) || (refId >= std::ssize(references_)) || (end <= beg) || (end <= 0)) {
        return {};
    }

    const std::int32_t normalizedBeg{std::max<std::int32_t>(beg, 0)};
    const std::int32_t normalizedEnd{end};

    const ReferenceIndex& ref{references_[refId]};

    // 1. Get overlapping bins
    const std::vector<std::uint16_t> overlappingBins{Reg2Bins(normalizedBeg, normalizedEnd)};

    // 2. Collect chunks from matching bins
    std::vector<Chunk> candidates;
    for (const std::uint16_t binNum : overlappingBins) {
        const auto it = ref.bins.find(binNum);
        if (it != ref.bins.end()) {
            candidates.insert(std::ranges::end(candidates), std::ranges::begin(it->second),
                              std::ranges::end(it->second));
        }
    }

    if (std::empty(candidates)) {
        return {};
    }

    // 3. Linear index pruning: find minimum offset for the query start window
    const std::size_t linearIdx = normalizedBeg / BAI_LINEAR_INDEX_WINDOW;
    // Discard chunks that end before the minimum offset
    if (linearIdx < std::size(ref.linearIndex)) {
        const VirtualOffset threshold{ref.linearIndex[linearIdx]};
        std::erase_if(candidates,
                      [threshold](const Chunk& chunk) { return chunk.End <= threshold; });
    }

    if (std::empty(candidates)) {
        return {};
    }

    // 4. Sort and merge
    return MergeChunks(candidates);
}

// --- Build ---

BaiIndex BaiIndex::Build(const std::filesystem::path& bamPath, std::size_t numWorkers)
{
    const std::unique_ptr<BlockSource> source{MakeBlockSource(bamPath, numWorkers)};
    BaiIndex index;

    // Step 1: Read and skip the BAM header
    std::vector<std::byte> headerBuf(MAX_DECOMPRESSED_BLOCK_SIZE * 4);
    std::size_t headerLen{0};

    // Track BGZF block offsets during header parsing
    std::uint64_t firstRecordBlockOffset{0};
    std::size_t firstRecordWithinBlock{0};
    std::size_t parsedHeaderSize{0};
    std::int32_t nRef{0};

    for (;;) {
        const std::optional<DecodedBlock> block{source->Next()};
        if (!block) {
            throw std::runtime_error{"Truncated BAM header"};
        }

        const std::uint64_t blockOffset{block->coffset};
        const std::size_t bytesRead{std::size(block->data)};

        if (headerLen + bytesRead > std::size(headerBuf)) {
            headerBuf.resize(headerLen + bytesRead + MAX_DECOMPRESSED_BLOCK_SIZE);
        }

        std::ranges::copy_n(std::data(block->data), bytesRead, std::data(headerBuf) + headerLen);
        headerLen += bytesRead;

        const std::size_t hdrSize{
            ComputeHeaderSize(std::span<const std::byte>{std::data(headerBuf), headerLen})};
        if (hdrSize > 0) {
            parsedHeaderSize = hdrSize;

            // Validate sort order — BAI indexing requires coordinate-sorted input
            const std::span<const std::byte> headerBytes{std::data(headerBuf), hdrSize};
            const auto headerResult{SamHeader::FromBamHeaderBlock(headerBytes)};
            if (headerResult && (headerResult->SortOrder() != "coordinate")) {
                throw std::runtime_error{
                    std::format("BaiIndex::Build requires coordinate-sorted BAM, got SO:{}",
                                headerResult->SortOrder())};
            }

            // Get n_ref from the header
            const std::uint32_t lText{ReadU32LE(std::data(headerBuf) + 4)};
            nRef = CheckedInt32(ReadU32LE(std::data(headerBuf) + 8 + lText), "n_ref");

            // Figure out where records start
            // The header ends at hdrSize bytes within the decompressed stream.
            // We need to figure out which block boundary that falls on.
            // Since we've been accumulating blocks, the remaining bytes after
            // the header are in the current block.
            const std::size_t remaining{headerLen - hdrSize};
            if (remaining > 0) {
                // Records start within this last block we read
                firstRecordBlockOffset = blockOffset;
                firstRecordWithinBlock = bytesRead - remaining;
            } else {
                // Records start at the next block. firstRecordBlockOffset is unused in this
                // case (headerLeftover == 0 → no initial segment), so the just-read block's
                // offset is a harmless placeholder.
                firstRecordBlockOffset = blockOffset;
                firstRecordWithinBlock = 0;
            }
            break;
        }
    }

    index.references_.resize(nRef);

    // Step 2: Re-open and seek to the first record position
    // We have leftover data from the header block. Rather than re-parsing,
    // let's process records from the decompressed data we already have.

    // We'll process BGZF blocks and track records within them.
    // We need to know: for each record, its virtual offset (block_offset << 16 |
    // within_block_offset)

    // Re-use the data we already have: the last block read during header parsing
    // may contain record data starting at firstRecordWithinBlock

    std::uint64_t currentBlockOffset{firstRecordBlockOffset};

    // Buffer for accumulating record data that spans blocks
    std::vector<std::byte> recordAccum;
    recordAccum.reserve(MAX_DECOMPRESSED_BLOCK_SIZE * 2);

    // Initialize with leftover data from header parsing
    const std::size_t headerLeftover{headerLen - parsedHeaderSize};
    if (headerLeftover > 0) {
        recordAccum.insert(std::ranges::end(recordAccum), std::data(headerBuf) + parsedHeaderSize,
                           std::data(headerBuf) + parsedHeaderSize + headerLeftover);
    }

    std::vector<BlockSegment> segments;
    if (headerLeftover > 0) {
        segments.push_back(BlockSegment{firstRecordBlockOffset,
                                        static_cast<std::uint16_t>(firstRecordWithinBlock),
                                        headerLeftover});
    }

    std::uint64_t mappedCount{0};
    std::uint64_t placedUnmappedCount{0};  // FLAG 0x4 reads that still carry a refId >= 0
    std::uint64_t noCoorCount{0};          // unplaced reads (refId < 0) → n_no_coor
    bool eof{false};
    std::optional<std::pair<std::int32_t, std::int32_t>> lastMappedCoordinate;
    bool sawUnmappedTailRecord{false};

    // Per-reference statistics for the optional metadata pseudo-bin (37450), enabling
    // samtools idxstats. begVo/endVo bracket every read placed on the reference.
    struct RefMeta
    {
        std::uint64_t nMapped{0};
        std::uint64_t nUnmapped{0};
        std::uint64_t begVo{std::numeric_limits<std::uint64_t>::max()};
        std::uint64_t endVo{0};
    };

    std::vector<RefMeta> refMeta(static_cast<std::size_t>(nRef));

    std::size_t accumPos{0};

    while (true) {
        // Ensure we have at least 4 bytes for block_size
        if (!EnsureAccumulatedBytes(*source, recordAccum, segments, currentBlockOffset, accumPos, 4,
                                    eof)) {
            break;
        }

        const std::uint32_t blockSize{ReadU32LE(std::data(recordAccum) + accumPos)};
        if (blockSize == 0) {
            break;
        }

        const std::size_t totalRecordBytes{4 + blockSize};

        // Ensure we have the complete record
        if (!EnsureAccumulatedBytes(*source, recordAccum, segments, currentBlockOffset, accumPos,
                                    totalRecordBytes, eof)) {
            break;
        }

        // Record virtual offset = position of the block_size field
        const VirtualOffset recordVo{VirtualOffsetAt(segments, accumPos)};
        // End virtual offset = position just past this record
        const VirtualOffset recordEndVo{VirtualOffsetAt(segments, accumPos + totalRecordBytes)};

        // Parse minimal record fields from the BAM binary data (after block_size)
        const std::byte* const rec{std::data(recordAccum) + accumPos + 4};
        const std::int32_t refId{ReadI32LE(rec)};
        const std::int32_t pos{ReadI32LE(rec + 4)};

        // bin_mq_nl at offset 8: bin(16) | mapq(8) | nl(8)
        const std::uint32_t binMqNl{ReadU32LE(rec + 8)};
        // flag_nc at offset 12: flag(16) | nc(16)
        const std::uint32_t flagNc{ReadU32LE(rec + 12)};
        const std::uint16_t flag = flagNc >> 16;
        const std::uint16_t nCigarOp = flagNc & 0xFFFF;
        const std::uint8_t nameLen = binMqNl & 0xFF;

        // Compute reference length from CIGAR
        std::int64_t refLen{0};
        if (nCigarOp > 0) {
            const std::size_t cigarOffset{BAM_FIXED_FIELDS_SIZE + nameLen};
            for (std::uint16_t ci{0}; ci < nCigarOp; ++ci) {
                const std::uint32_t cigarVal{ReadU32LE(rec + cigarOffset + ci * 4)};
                const std::uint8_t opCode{static_cast<std::uint8_t>(cigarVal & 0xFU)};
                const std::uint32_t opLen{cigarVal >> 4};
                if (ConsumesReference(opCode)) {
                    refLen += opLen;
                }
            }
        }

        const bool isUnmapped{(flag & 0x4) != 0};

        if (refId < 0) {
            // Unplaced read (no coordinate): the true n_no_coor population, which always
            // sorts after every placed read. htslib counts only these in n_no_coor.
            sawUnmappedTailRecord = true;
            ++noCoorCount;
        } else if (refId < nRef) {
            // Placed read (refId >= 0): either mapped, or placed-unmapped (FLAG 0x4 carrying
            // a coordinate). Both participate in coordinate ordering.
            if (sawUnmappedTailRecord) {
                throw std::runtime_error{
                    std::format("BaiIndex::Build requires coordinate-sorted BAM "
                                "payload; saw placed "
                                "refId={} pos={} after unplaced tail in {}",
                                refId, pos, bamPath.string())};
            }
            const std::pair<std::int32_t, std::int32_t> currentCoordinate{refId, pos};
            if (lastMappedCoordinate && currentCoordinate < *lastMappedCoordinate) {
                throw std::runtime_error{
                    std::format("BaiIndex::Build requires coordinate-sorted BAM "
                                "payload; saw refId={} pos={} "
                                "after refId={} pos={} in {}",
                                refId, pos, lastMappedCoordinate->first,
                                lastMappedCoordinate->second, bamPath.string())};
            }
            lastMappedCoordinate = currentCoordinate;

            RefMeta& meta{refMeta[static_cast<std::size_t>(refId)]};
            meta.begVo = std::min(meta.begVo, recordVo.Value());
            meta.endVo = std::max(meta.endVo, recordEndVo.Value());

            if (isUnmapped) {
                // Placed but unmapped: counted as per-reference n_unmapped (not n_no_coor),
                // and not added to a bin since it has no alignment span.
                ++placedUnmappedCount;
                ++meta.nUnmapped;
            } else {
                ++mappedCount;
                ++meta.nMapped;
                ReferenceIndex& refIdx{index.references_[refId]};

                // Compute bin
                const std::int32_t endPos{
                    CheckedInt32(static_cast<std::int64_t>(pos) + refLen, "alignment end")};
                const std::int32_t nonEmptyEnd{NonEmptyAlignmentEnd(pos, endPos)};
                const std::uint16_t bin{Reg2Bin(pos, nonEmptyEnd)};

                // Add chunk to bin
                refIdx.bins[bin].push_back(Chunk{recordVo, recordEndVo});

                // Update linear index
                const std::int32_t begWindow = pos / BAI_LINEAR_INDEX_WINDOW;
                const std::int32_t endWindow{
                    static_cast<std::int32_t>((nonEmptyEnd - 1) / BAI_LINEAR_INDEX_WINDOW)};
                const std::size_t maxWindow = endWindow + 1;
                if (maxWindow > std::size(refIdx.linearIndex)) {
                    refIdx.linearIndex.resize(maxWindow);
                }
                for (std::int32_t w{begWindow}; w <= endWindow; ++w) {
                    VirtualOffset& entry{refIdx.linearIndex[w]};
                    if (entry.Value() == 0 || recordVo < entry) {
                        entry = recordVo;
                    }
                }
            }
        }

        accumPos += totalRecordBytes;

        // Compact the accumulator periodically to avoid unbounded growth
        if (accumPos > MAX_DECOMPRESSED_BLOCK_SIZE) {
            CompactAccumulator(recordAccum, segments, accumPos);
        }
    }

    // Merge chunks within each bin
    for (ReferenceIndex& ref : index.references_) {
        for (auto& [binNum, chunks] : ref.bins) {
            chunks = MergeChunks(chunks);
        }
    }

    // Emit the metadata pseudo-bin (37450) per reference. Added after the merge above so
    // its two stat "chunks" — {begVo,endVo} and {n_mapped,n_unmapped} — are never merged.
    for (std::int32_t r{0}; r < nRef; ++r) {
        const RefMeta& meta{refMeta[static_cast<std::size_t>(r)]};
        if ((meta.nMapped + meta.nUnmapped) == 0) {
            continue;
        }
        std::vector<Chunk>& metaBin{index.references_[r].bins[BAI_METADATA_BIN]};
        metaBin.push_back(Chunk{VirtualOffset{meta.begVo}, VirtualOffset{meta.endVo}});
        metaBin.push_back(Chunk{VirtualOffset{meta.nMapped}, VirtualOffset{meta.nUnmapped}});
    }

    index.mappedCount_ = mappedCount;
    index.unmappedCount_ = placedUnmappedCount + noCoorCount;
    index.noCoorCount_ = noCoorCount;

    return index;
}

}  // namespace Samoa
}  // namespace PacBio
