#include <pbsamoa/index/BaiIndex.hpp>

#include "BinaryUtils.hpp"

#include <pbsamoa/core/Bgzf.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/io/BamRawReader.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <filesystem>
#include <format>
#include <fstream>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>

#include <cassert>
#include <cstdint>

namespace PacBio {
namespace Samoa {

namespace {

constexpr std::size_t BAI_LINEAR_INDEX_WINDOW{16384};  // 2^14 = 16 kbp

void ReadExact(std::ifstream& in, std::span<std::byte> bytes, std::string_view fieldName)
{
    in.read(reinterpret_cast<char*>(std::data(bytes)),
            static_cast<std::streamsize>(std::size(bytes)));
    if (!in) {
        throw std::runtime_error{std::format("Truncated BAI file reading {}", fieldName)};
    }
}

std::uint32_t ReadU32LEFromFile(std::ifstream& in, std::string_view fieldName)
{
    std::array<std::byte, 4> bytes{};
    ReadExact(in, bytes, fieldName);
    return ReadU32LE(std::data(bytes));
}

std::int32_t ReadI32LEFromFile(std::ifstream& in, std::string_view fieldName)
{
    return std::bit_cast<std::int32_t>(ReadU32LEFromFile(in, fieldName));
}

std::uint64_t ReadU64LEFromFile(std::ifstream& in, std::string_view fieldName)
{
    std::array<std::byte, 8> bytes{};
    ReadExact(in, bytes, fieldName);
    return ReadU64LE(std::data(bytes));
}

std::optional<std::uint64_t> ReadOptionalU64LEFromFile(std::ifstream& in,
                                                       std::string_view fieldName)
{
    std::array<std::byte, 8> bytes{};
    in.read(reinterpret_cast<char*>(std::data(bytes)),
            static_cast<std::streamsize>(std::size(bytes)));
    const std::streamsize bytesRead{in.gcount()};
    if (bytesRead == 0) {
        return std::nullopt;
    }
    if (bytesRead != static_cast<std::streamsize>(std::size(bytes))) {
        throw std::runtime_error{std::format("Truncated BAI file reading {}", fieldName)};
    }
    return ReadU64LE(std::data(bytes));
}

std::size_t CheckedNonNegativeCount(std::int32_t value, std::string_view fieldName)
{
    if (value < 0) {
        throw std::runtime_error{
            std::format("Invalid BAI file: {} must be >= 0, got {}", fieldName, value)};
    }
    return static_cast<std::size_t>(value);
}

void WriteU32LE(std::ofstream& out, std::uint32_t v)
{
    const std::array<std::byte, 4> bytes{
        static_cast<std::byte>(v & 0xFFU),
        static_cast<std::byte>((v >> 8U) & 0xFFU),
        static_cast<std::byte>((v >> 16U) & 0xFFU),
        static_cast<std::byte>((v >> 24U) & 0xFFU),
    };
    out.write(reinterpret_cast<const char*>(std::data(bytes)),
              static_cast<std::streamsize>(std::size(bytes)));
}

void WriteI32LE(std::ofstream& out, std::int32_t v)
{
    WriteU32LE(out, std::bit_cast<std::uint32_t>(v));
}

void WriteU64LE(std::ofstream& out, std::uint64_t v)
{
    const std::array<std::byte, 8> bytes{
        static_cast<std::byte>(v & 0xFFULL),
        static_cast<std::byte>((v >> 8ULL) & 0xFFULL),
        static_cast<std::byte>((v >> 16ULL) & 0xFFULL),
        static_cast<std::byte>((v >> 24ULL) & 0xFFULL),
        static_cast<std::byte>((v >> 32ULL) & 0xFFULL),
        static_cast<std::byte>((v >> 40ULL) & 0xFFULL),
        static_cast<std::byte>((v >> 48ULL) & 0xFFULL),
        static_cast<std::byte>((v >> 56ULL) & 0xFFULL),
    };
    out.write(reinterpret_cast<const char*>(std::data(bytes)),
              static_cast<std::streamsize>(std::size(bytes)));
}

/// \brief Merge overlapping/adjacent chunks. Input must be sorted by Begin.
std::vector<Chunk> MergeChunks(std::vector<Chunk>& chunks)
{
    if (std::empty(chunks)) {
        return {};
    }

    std::ranges::sort(chunks, {}, &Chunk::Begin);

    std::vector<Chunk> merged;
    merged.push_back(chunks[0]);

    for (std::size_t i{1}; i < std::size(chunks); ++i) {
        Chunk& last{merged.back()};
        if (chunks[i].Begin <= last.End) {
            if (chunks[i].End > last.End) {
                last.End = chunks[i].End;
            }
        } else {
            merged.push_back(chunks[i]);
        }
    }

    return merged;
}

}  // namespace

// --- Chunk ---

bool Chunk::Overlaps(const Chunk& other) const
{
    return (Begin <= other.End) && (other.Begin <= End);
}

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
    if (magic[0] != static_cast<std::byte>('B') || magic[1] != static_cast<std::byte>('A') ||
        magic[2] != static_cast<std::byte>('I') || magic[3] != static_cast<std::byte>('\1')) {
        throw std::runtime_error{"Invalid BAI magic in: " + path.string()};
    }

    // Read n_ref
    const std::int32_t nRef{ReadI32LEFromFile(file, "n_ref")};
    const std::size_t nRefCount{CheckedNonNegativeCount(nRef, "n_ref")};

    BaiIndex index;
    index.references_.resize(nRefCount);

    for (std::size_t r{0}; r < nRefCount; ++r) {
        ReferenceIndex& ref{index.references_[r]};

        // n_bin
        const std::int32_t nBin{ReadI32LEFromFile(file, "n_bin")};
        const std::size_t nBinCount{CheckedNonNegativeCount(nBin, "n_bin")};

        for (std::size_t b{0}; b < nBinCount; ++b) {
            const std::uint32_t binNumber{ReadU32LEFromFile(file, "bin number")};
            const std::int32_t nChunks{ReadI32LEFromFile(file, "n_chunks")};
            const std::size_t nChunkCount{CheckedNonNegativeCount(nChunks, "n_chunks")};

            std::vector<Chunk> chunks;
            chunks.reserve(nChunkCount);
            for (std::size_t c{0}; c < nChunkCount; ++c) {
                const std::uint64_t chunkBeg{ReadU64LEFromFile(file, "chunk begin")};
                const std::uint64_t chunkEnd{ReadU64LEFromFile(file, "chunk end")};
                chunks.push_back(Chunk{VirtualOffset{chunkBeg}, VirtualOffset{chunkEnd}});
            }

            ref.bins[binNumber] = std::move(chunks);
        }

        // n_intv
        const std::int32_t nIntv{ReadI32LEFromFile(file, "n_intv")};
        const std::size_t nIntvCount{CheckedNonNegativeCount(nIntv, "n_intv")};

        ref.linearIndex.resize(nIntvCount);
        for (std::size_t i{0}; i < nIntvCount; ++i) {
            const std::uint64_t offset{ReadU64LEFromFile(file, "linear index offset")};
            ref.linearIndex[i] = VirtualOffset{offset};
        }
    }

    // Optional: read n_no_coor (unmapped count) at end of file
    const std::optional<std::uint64_t> nNoCoor{ReadOptionalU64LEFromFile(file, "n_no_coor")};
    if (nNoCoor.has_value()) {
        index.unmappedCount_ = *nNoCoor;
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
        const std::int32_t nRef = std::ssize(references_);
        WriteI32LE(file, nRef);

        for (const ReferenceIndex& ref : references_) {
            // n_bin
            const std::int32_t nBin = std::ssize(ref.bins);
            WriteI32LE(file, nBin);

            for (const auto& [binNumber, chunks] : ref.bins) {
                WriteU32LE(file, binNumber);
                const std::int32_t nChunks = std::ssize(chunks);
                WriteI32LE(file, nChunks);
                for (const Chunk& chunk : chunks) {
                    WriteU64LE(file, chunk.Begin.Value());
                    WriteU64LE(file, chunk.End.Value());
                }
            }

            // n_intv
            const std::int32_t nIntv = std::ssize(ref.linearIndex);
            WriteI32LE(file, nIntv);
            for (const VirtualOffset& offset : ref.linearIndex) {
                WriteU64LE(file, offset.Value());
            }
        }

        // n_no_coor
        WriteU64LE(file, unmappedCount_);
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
    const std::int32_t normalizedEnd{std::max<std::int32_t>(end, 0)};
    if (normalizedEnd <= normalizedBeg) {
        return {};
    }

    const ReferenceIndex& ref{references_[refId]};

    // 1. Get overlapping bins
    const std::vector<std::uint16_t> overlappingBins{Reg2Bins(normalizedBeg, normalizedEnd)};

    // 2. Collect chunks from matching bins
    std::vector<Chunk> candidates;
    for (const std::uint16_t binNum : overlappingBins) {
        const auto it = ref.bins.find(binNum);
        if (it != std::ranges::end(ref.bins)) {
            candidates.insert(std::ranges::end(candidates), std::ranges::begin(it->second),
                              std::ranges::end(it->second));
        }
    }

    if (std::empty(candidates)) {
        return {};
    }

    // 3. Linear index pruning: find minimum offset for the query start window
    const std::size_t linearIdx = normalizedBeg / BAI_LINEAR_INDEX_WINDOW;
    VirtualOffset minOffset;
    if (linearIdx < std::size(ref.linearIndex)) {
        minOffset = ref.linearIndex[linearIdx];
    }

    // Discard chunks that end before the minimum offset
    std::erase_if(candidates, [&minOffset](const Chunk& c) { return c.End < minOffset; });

    if (std::empty(candidates)) {
        return {};
    }

    // 4. Sort and merge
    return MergeChunks(candidates);
}

// --- Build ---

BaiIndex BaiIndex::Build(const std::filesystem::path& bamPath)
{
    BgzfReader bgzf{bamPath};
    BaiIndex index;

    // Step 1: Read and skip the BAM header
    std::vector<std::byte> headerBuf(MAX_DECOMPRESSED_BLOCK_SIZE * 4);
    std::size_t headerLen{0};
    std::vector<std::byte> blockBuf(MAX_DECOMPRESSED_BLOCK_SIZE);

    // Track BGZF block offsets during header parsing
    std::uint64_t firstRecordBlockOffset{0};
    std::size_t firstRecordWithinBlock{0};
    bool headerParsed{false};
    std::int32_t nRef{0};

    while (!headerParsed) {
        const std::optional<std::size_t> bytesRead{bgzf.ReadBlock(std::span<std::byte>{blockBuf})};

        if (!bytesRead.has_value() || *bytesRead == 0) {
            throw std::runtime_error{"Truncated BAM header"};
        }

        const std::uint64_t blockOffset{bgzf.Tell().BlockOffset()};

        if (headerLen + *bytesRead > std::size(headerBuf)) {
            headerBuf.resize(headerLen + *bytesRead + MAX_DECOMPRESSED_BLOCK_SIZE);
        }

        std::ranges::copy_n(std::data(blockBuf), *bytesRead, std::data(headerBuf) + headerLen);
        headerLen += *bytesRead;

        const std::size_t hdrSize{
            ComputeHeaderSize(std::span<const std::byte>{std::data(headerBuf), headerLen})};
        if (hdrSize > 0) {
            // Validate sort order — BAI indexing requires coordinate-sorted input
            auto headerResult{SamHeader::FromBamHeaderBlock(
                std::span<const std::byte>{std::data(headerBuf), hdrSize})};
            if (headerResult.has_value() && headerResult->SortOrder() != "coordinate") {
                throw std::runtime_error{
                    std::format("BaiIndex::Build requires coordinate-sorted BAM, got SO:{}",
                                headerResult->SortOrder())};
            }

            // Get n_ref from the header
            const std::uint32_t lText{ReadU32LE(std::data(headerBuf) + 4)};
            nRef = ReadU32LE(std::data(headerBuf) + 8 + lText);

            // Figure out where records start
            // The header ends at hdrSize bytes within the decompressed stream.
            // We need to figure out which block boundary that falls on.
            // Since we've been accumulating blocks, the remaining bytes after
            // the header are in the current block.
            const std::size_t remaining{headerLen - hdrSize};
            if (remaining > 0) {
                // Records start within this last block we read
                firstRecordBlockOffset = blockOffset;
                firstRecordWithinBlock = *bytesRead - remaining;
            } else {
                // Records start at the next block
                firstRecordBlockOffset = bgzf.Tell().BlockOffset();
                firstRecordWithinBlock = 0;
            }
            headerParsed = true;
        }
    }

    index.references_.resize(nRef);

    // Step 2: Re-open and seek to the first record position
    // We have leftover data from the header block. Rather than re-parsing,
    // let's process records from the decompressed data we already have.

    // We'll process BGZF blocks and track records within them.
    // We need to know: for each record, its virtual offset (block_offset << 16 | within_block_offset)

    // Re-use the data we already have: the last block read during header parsing
    // may contain record data starting at firstRecordWithinBlock

    std::uint64_t currentBlockOffset{firstRecordBlockOffset};

    // Buffer for accumulating record data that spans blocks
    std::vector<std::byte> recordAccum;
    recordAccum.reserve(MAX_DECOMPRESSED_BLOCK_SIZE * 2);

    // Initialize with leftover data from header parsing
    const std::size_t headerLeftover{
        headerLen - ComputeHeaderSize(std::span<const std::byte>{std::data(headerBuf), headerLen})};
    if (headerLeftover > 0) {
        const std::size_t hdrSize{
            ComputeHeaderSize(std::span<const std::byte>{std::data(headerBuf), headerLen})};
        recordAccum.insert(std::ranges::end(recordAccum), std::data(headerBuf) + hdrSize,
                           std::data(headerBuf) + hdrSize + headerLeftover);
    }

    // We need to track which parts of recordAccum came from which blocks
    // to compute virtual offsets correctly.
    // Approach: track (blockOffset, withinBlockStartPos, byteCount) for each segment in accumulator
    struct BlockSegment
    {
        std::uint64_t blockOffset;
        std::uint16_t withinBlockStart;
        std::size_t byteCount;
    };

    std::vector<BlockSegment> segments;
    if (headerLeftover > 0) {
        segments.push_back(BlockSegment{firstRecordBlockOffset,
                                        static_cast<std::uint16_t>(firstRecordWithinBlock),
                                        headerLeftover});
    }

    std::uint64_t mappedCount{0};
    std::uint64_t unmappedCount{0};
    bool eof{false};

    auto readMoreData = [&]() -> bool {
        if (eof) {
            return false;
        }
        const std::optional<std::size_t> bytesRead{bgzf.ReadBlock(std::span<std::byte>{blockBuf})};
        currentBlockOffset = bgzf.Tell().BlockOffset();
        if (!bytesRead.has_value() || *bytesRead == 0) {
            eof = true;
            return false;
        }
        recordAccum.insert(std::ranges::end(recordAccum), std::data(blockBuf),
                           std::data(blockBuf) + *bytesRead);
        segments.push_back(BlockSegment{currentBlockOffset, 0, *bytesRead});
        return true;
    };

    // Helper: compute virtual offset for a position within recordAccum
    auto virtualOffsetAt = [&](std::size_t accumPos) -> VirtualOffset {
        std::size_t pos{0};
        for (const BlockSegment& seg : segments) {
            if (accumPos < pos + seg.byteCount) {
                const std::size_t within{seg.withinBlockStart + (accumPos - pos)};
                assert(within <= std::numeric_limits<std::uint16_t>::max());
                return VirtualOffset{seg.blockOffset, static_cast<std::uint16_t>(within)};
            }
            pos += seg.byteCount;
        }
        // Past end — return the next block position
        if (!std::empty(segments)) {
            const BlockSegment& last{segments.back()};
            const std::size_t pastEnd{last.withinBlockStart + last.byteCount};
            assert(pastEnd <= std::numeric_limits<std::uint16_t>::max());
            return VirtualOffset{last.blockOffset, static_cast<std::uint16_t>(pastEnd)};
        }
        return VirtualOffset{};
    };

    std::size_t accumPos{0};

    const auto ensureBytes = [&](std::size_t n) -> bool {
        while ((std::size(recordAccum) - accumPos) < n) {
            if (!readMoreData()) {
                return false;
            }
        }
        return true;
    };

    while (true) {
        // Ensure we have at least 4 bytes for block_size
        if (!ensureBytes(4)) {
            break;
        }

        const std::uint32_t blockSize{ReadU32LE(std::data(recordAccum) + accumPos)};
        if (blockSize == 0) {
            break;
        }

        const std::size_t totalRecordBytes{4 + blockSize};

        // Ensure we have the complete record
        if (!ensureBytes(totalRecordBytes)) {
            break;
        }

        // Record virtual offset = position of the block_size field
        const VirtualOffset recordVo{virtualOffsetAt(accumPos)};
        // End virtual offset = position just past this record
        const VirtualOffset recordEndVo{virtualOffsetAt(accumPos + totalRecordBytes)};

        // Parse minimal record fields from the BAM binary data (after block_size)
        const std::byte* rec{std::data(recordAccum) + accumPos + 4};
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
            const std::size_t cigarOffset{32U + nameLen};
            for (std::uint16_t ci{0}; ci < nCigarOp; ++ci) {
                const std::uint32_t cigarVal{ReadU32LE(rec + cigarOffset + ci * 4)};
                const std::uint8_t opCode = cigarVal & 0xFU;
                const std::uint32_t opLen{cigarVal >> 4};
                // Ops that consume reference: M(0), D(2), N(3), =(7), X(8)
                if ((opCode == 0) || (opCode == 2) || (opCode == 3) || (opCode == 7) ||
                    (opCode == 8)) {
                    refLen += opLen;
                }
            }
        }

        const bool isUnmapped{(flag & 0x4) != 0};

        if (isUnmapped) {
            ++unmappedCount;
        } else if ((refId >= 0) && (refId < nRef)) {
            ++mappedCount;
            ReferenceIndex& refIdx{index.references_[refId]};

            // Compute bin
            const std::int32_t endPos = pos + refLen;
            const std::uint16_t bin{Reg2Bin(pos, endPos > pos ? endPos : pos + 1)};

            // Add chunk to bin
            refIdx.bins[bin].push_back(Chunk{recordVo, recordEndVo});

            // Update linear index
            const std::int32_t begWindow = pos / BAI_LINEAR_INDEX_WINDOW;
            const std::int32_t endWindow =
                ((endPos > pos) ? (endPos - 1) : pos) / BAI_LINEAR_INDEX_WINDOW;
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

        accumPos += totalRecordBytes;

        // Compact the accumulator periodically to avoid unbounded growth
        if (accumPos > MAX_DECOMPRESSED_BLOCK_SIZE) {
            // Remove consumed segments
            std::size_t consumed{0};
            std::ptrdiff_t segIdx{0};
            while (segIdx < std::ssize(segments) &&
                   (consumed + segments[segIdx].byteCount) <= accumPos) {
                consumed += segments[segIdx].byteCount;
                ++segIdx;
            }
            if (segIdx > 0) {
                // Partial segment at boundary
                if ((consumed < accumPos) && (segIdx < std::ssize(segments))) {
                    // The segment at segIdx straddles accumPos
                    const std::size_t skipInSeg{accumPos - consumed};
                    const std::size_t newWithin{segments[segIdx].withinBlockStart + skipInSeg};
                    assert(newWithin <= std::numeric_limits<std::uint16_t>::max());
                    segments[segIdx].withinBlockStart = static_cast<std::uint16_t>(newWithin);
                    segments[segIdx].byteCount -= skipInSeg;
                } else if (consumed > accumPos) {
                    // This shouldn't happen if segments are contiguous
                } else {
                    // consumed == accumPos, clean boundary
                }
                segments.erase(std::ranges::begin(segments), std::ranges::begin(segments) + segIdx);
            }
            recordAccum.erase(
                std::ranges::begin(recordAccum),
                std::ranges::begin(recordAccum) + static_cast<std::ptrdiff_t>(accumPos));
            accumPos = 0;
        }
    }

    // Merge chunks within each bin
    for (ReferenceIndex& ref : index.references_) {
        for (auto& [binNum, chunks] : ref.bins) {
            chunks = MergeChunks(chunks);
        }
    }

    index.mappedCount_ = mappedCount;
    index.unmappedCount_ = unmappedCount;

    return index;
}

}  // namespace Samoa
}  // namespace PacBio
