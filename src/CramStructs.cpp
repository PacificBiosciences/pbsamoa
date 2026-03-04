#include <pbsamoa/cram/CramCompression.hpp>
#include <pbsamoa/cram/CramStructs.hpp>

#include "BinaryUtils.hpp"
#include "CramInternal.hpp"

#include <libdeflate.h>

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace PacBio {
namespace Samoa {

namespace {

std::uint32_t ComputeCrc32(std::span<const std::byte> data)
{
    return libdeflate_crc32(0, data.data(), std::size(data));
}

std::size_t Itf8EncodedSize(const std::int32_t value)
{
    const auto uval = static_cast<std::uint32_t>(value);
    if (uval < 0x80) {
        return 1;
    }
    if (uval < 0x4000) {
        return 2;
    }
    if (uval < 0x200000) {
        return 3;
    }
    if (uval < 0x10000000) {
        return 4;
    }
    return 5;
}

std::size_t SerializedBlockSize(const CramBlock& block)
{
    const auto dataSize = std::size(block.Data);
    if (dataSize > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
        throw std::runtime_error("block: compressed size exceeds ITF8 range");
    }
    const auto compressedSize = static_cast<std::int32_t>(dataSize);
    return 2 + Itf8EncodedSize(block.ContentId) + Itf8EncodedSize(compressedSize) +
           Itf8EncodedSize(block.RawSize) + dataSize + 4;
}

void ValidateAndAdvance(std::int64_t& accumulator, const std::int64_t delta, const char* context)
{
    const auto next = accumulator + delta;
    if (next < 0 || next > std::numeric_limits<std::int32_t>::max()) {
        throw std::runtime_error(std::string{"CRAM size overflow while "} + context);
    }
    accumulator = next;
}

}  // namespace

// ---------------------------------------------------------------------------
// ITF-8
// ---------------------------------------------------------------------------

std::int32_t ReadItf8(std::span<const std::byte> data, std::size_t& bytesRead)
{
    if (std::empty(data)) {
        throw std::runtime_error("ReadItf8: empty data");
    }

    const auto b0 = static_cast<std::uint8_t>(data[0]);

    if ((b0 & 0x80) == 0) {
        // 0xxxxxxx
        bytesRead = 1;
        return static_cast<std::int32_t>(b0);
    }
    if ((b0 & 0xC0) == 0x80) {
        // 10xxxxxx + 1 byte
        if (std::size(data) < 2) {
            throw std::runtime_error("ReadItf8: truncated 2-byte value");
        }
        bytesRead = 2;
        const auto value = ((static_cast<std::uint32_t>(b0 & 0x3F) << 8) |
                            static_cast<std::uint32_t>(static_cast<std::uint8_t>(data[1])));
        return static_cast<std::int32_t>(value);
    }
    if ((b0 & 0xE0) == 0xC0) {
        // 110xxxxx + 2 bytes
        if (std::size(data) < 3) {
            throw std::runtime_error("ReadItf8: truncated 3-byte value");
        }
        bytesRead = 3;
        const auto value = ((static_cast<std::uint32_t>(b0 & 0x1F) << 16) |
                            (static_cast<std::uint32_t>(static_cast<std::uint8_t>(data[1])) << 8) |
                            static_cast<std::uint32_t>(static_cast<std::uint8_t>(data[2])));
        return static_cast<std::int32_t>(value);
    }
    if ((b0 & 0xF0) == 0xE0) {
        // 1110xxxx + 3 bytes
        if (std::size(data) < 4) {
            throw std::runtime_error("ReadItf8: truncated 4-byte value");
        }
        bytesRead = 4;
        const auto value = ((static_cast<std::uint32_t>(b0 & 0x0F) << 24) |
                            (static_cast<std::uint32_t>(static_cast<std::uint8_t>(data[1])) << 16) |
                            (static_cast<std::uint32_t>(static_cast<std::uint8_t>(data[2])) << 8) |
                            static_cast<std::uint32_t>(static_cast<std::uint8_t>(data[3])));
        return static_cast<std::int32_t>(value);
    }
    // 1111xxxx + 4 bytes (only lower 4 bits of last byte used)
    if (std::size(data) < 5) {
        throw std::runtime_error("ReadItf8: truncated 5-byte value");
    }
    bytesRead = 5;
    const auto value = ((static_cast<std::uint32_t>(b0 & 0x0F) << 28) |
                        (static_cast<std::uint32_t>(static_cast<std::uint8_t>(data[1])) << 20) |
                        (static_cast<std::uint32_t>(static_cast<std::uint8_t>(data[2])) << 12) |
                        (static_cast<std::uint32_t>(static_cast<std::uint8_t>(data[3])) << 4) |
                        (static_cast<std::uint32_t>(static_cast<std::uint8_t>(data[4])) & 0x0F));
    return static_cast<std::int32_t>(value);
}

std::int64_t ReadLtf8(std::span<const std::byte> data, std::size_t& bytesRead)
{
    if (std::empty(data)) {
        throw std::runtime_error("ReadLtf8: empty data");
    }

    const auto b0 = static_cast<std::uint8_t>(data[0]);

    if ((b0 & 0x80) == 0) {
        bytesRead = 1;
        return static_cast<std::int64_t>(b0);
    }
    if ((b0 & 0xC0) == 0x80) {
        if (std::size(data) < 2) {
            throw std::runtime_error("ReadLtf8: truncated");
        }
        bytesRead = 2;
        return (static_cast<std::int64_t>(b0 & 0x3F) << 8) |
               static_cast<std::int64_t>(static_cast<std::uint8_t>(data[1]));
    }
    if ((b0 & 0xE0) == 0xC0) {
        if (std::size(data) < 3) {
            throw std::runtime_error("ReadLtf8: truncated");
        }
        bytesRead = 3;
        return (static_cast<std::int64_t>(b0 & 0x1F) << 16) |
               (static_cast<std::int64_t>(static_cast<std::uint8_t>(data[1])) << 8) |
               static_cast<std::int64_t>(static_cast<std::uint8_t>(data[2]));
    }
    if ((b0 & 0xF0) == 0xE0) {
        if (std::size(data) < 4) {
            throw std::runtime_error("ReadLtf8: truncated");
        }
        bytesRead = 4;
        return (static_cast<std::int64_t>(b0 & 0x0F) << 24) |
               (static_cast<std::int64_t>(static_cast<std::uint8_t>(data[1])) << 16) |
               (static_cast<std::int64_t>(static_cast<std::uint8_t>(data[2])) << 8) |
               static_cast<std::int64_t>(static_cast<std::uint8_t>(data[3]));
    }
    if ((b0 & 0xF8) == 0xF0) {
        if (std::size(data) < 5) {
            throw std::runtime_error("ReadLtf8: truncated");
        }
        bytesRead = 5;
        return (static_cast<std::int64_t>(b0 & 0x07) << 32) |
               (static_cast<std::int64_t>(static_cast<std::uint8_t>(data[1])) << 24) |
               (static_cast<std::int64_t>(static_cast<std::uint8_t>(data[2])) << 16) |
               (static_cast<std::int64_t>(static_cast<std::uint8_t>(data[3])) << 8) |
               static_cast<std::int64_t>(static_cast<std::uint8_t>(data[4]));
    }
    if ((b0 & 0xFC) == 0xF8) {
        if (std::size(data) < 6) {
            throw std::runtime_error("ReadLtf8: truncated");
        }
        bytesRead = 6;
        return (static_cast<std::int64_t>(b0 & 0x03) << 40) |
               (static_cast<std::int64_t>(static_cast<std::uint8_t>(data[1])) << 32) |
               (static_cast<std::int64_t>(static_cast<std::uint8_t>(data[2])) << 24) |
               (static_cast<std::int64_t>(static_cast<std::uint8_t>(data[3])) << 16) |
               (static_cast<std::int64_t>(static_cast<std::uint8_t>(data[4])) << 8) |
               static_cast<std::int64_t>(static_cast<std::uint8_t>(data[5]));
    }
    if ((b0 & 0xFE) == 0xFC) {
        if (std::size(data) < 7) {
            throw std::runtime_error("ReadLtf8: truncated");
        }
        bytesRead = 7;
        return (static_cast<std::int64_t>(b0 & 0x01) << 48) |
               (static_cast<std::int64_t>(static_cast<std::uint8_t>(data[1])) << 40) |
               (static_cast<std::int64_t>(static_cast<std::uint8_t>(data[2])) << 32) |
               (static_cast<std::int64_t>(static_cast<std::uint8_t>(data[3])) << 24) |
               (static_cast<std::int64_t>(static_cast<std::uint8_t>(data[4])) << 16) |
               (static_cast<std::int64_t>(static_cast<std::uint8_t>(data[5])) << 8) |
               static_cast<std::int64_t>(static_cast<std::uint8_t>(data[6]));
    }
    if (b0 == 0xFE) {
        if (std::size(data) < 8) {
            throw std::runtime_error("ReadLtf8: truncated");
        }
        bytesRead = 8;
        return (static_cast<std::int64_t>(static_cast<std::uint8_t>(data[1])) << 48) |
               (static_cast<std::int64_t>(static_cast<std::uint8_t>(data[2])) << 40) |
               (static_cast<std::int64_t>(static_cast<std::uint8_t>(data[3])) << 32) |
               (static_cast<std::int64_t>(static_cast<std::uint8_t>(data[4])) << 24) |
               (static_cast<std::int64_t>(static_cast<std::uint8_t>(data[5])) << 16) |
               (static_cast<std::int64_t>(static_cast<std::uint8_t>(data[6])) << 8) |
               static_cast<std::int64_t>(static_cast<std::uint8_t>(data[7]));
    }
    // 0xFF + 8 bytes
    if (std::size(data) < 9) {
        throw std::runtime_error("ReadLtf8: truncated");
    }
    bytesRead = 9;
    std::int64_t result = 0;
    for (int i = 1; i <= 8; ++i) {
        result = (result << 8) | static_cast<std::int64_t>(static_cast<std::uint8_t>(data[i]));
    }
    return result;
}

std::size_t WriteItf8(std::vector<std::byte>& out, std::int32_t value)
{
    const auto uval = static_cast<std::uint32_t>(value);

    if (uval < 0x80) {
        out.push_back(static_cast<std::byte>(uval));
        return 1;
    }
    if (uval < 0x4000) {
        out.push_back(static_cast<std::byte>(0x80 | (uval >> 8)));
        out.push_back(static_cast<std::byte>(uval & 0xFF));
        return 2;
    }
    if (uval < 0x200000) {
        out.push_back(static_cast<std::byte>(0xC0 | (uval >> 16)));
        out.push_back(static_cast<std::byte>((uval >> 8) & 0xFF));
        out.push_back(static_cast<std::byte>(uval & 0xFF));
        return 3;
    }
    if (uval < 0x10000000) {
        out.push_back(static_cast<std::byte>(0xE0 | (uval >> 24)));
        out.push_back(static_cast<std::byte>((uval >> 16) & 0xFF));
        out.push_back(static_cast<std::byte>((uval >> 8) & 0xFF));
        out.push_back(static_cast<std::byte>(uval & 0xFF));
        return 4;
    }
    out.push_back(static_cast<std::byte>(0xF0 | ((uval >> 28) & 0x0F)));
    out.push_back(static_cast<std::byte>((uval >> 20) & 0xFF));
    out.push_back(static_cast<std::byte>((uval >> 12) & 0xFF));
    out.push_back(static_cast<std::byte>((uval >> 4) & 0xFF));
    out.push_back(static_cast<std::byte>(uval & 0x0F));
    return 5;
}

std::size_t WriteLtf8(std::vector<std::byte>& out, std::int64_t value)
{
    const std::uint64_t uval = value;

    if (uval < 0x80) {
        out.push_back(static_cast<std::byte>(uval));
        return 1;
    }
    if (uval < 0x4000) {
        out.push_back(static_cast<std::byte>(0x80 | (uval >> 8)));
        out.push_back(static_cast<std::byte>(uval & 0xFF));
        return 2;
    }
    if (uval < 0x200000) {
        out.push_back(static_cast<std::byte>(0xC0 | (uval >> 16)));
        out.push_back(static_cast<std::byte>((uval >> 8) & 0xFF));
        out.push_back(static_cast<std::byte>(uval & 0xFF));
        return 3;
    }
    if (uval < 0x10000000) {
        out.push_back(static_cast<std::byte>(0xE0 | (uval >> 24)));
        out.push_back(static_cast<std::byte>((uval >> 16) & 0xFF));
        out.push_back(static_cast<std::byte>((uval >> 8) & 0xFF));
        out.push_back(static_cast<std::byte>(uval & 0xFF));
        return 4;
    }
    if (uval < 0x800000000ULL) {
        out.push_back(static_cast<std::byte>(0xF0 | ((uval >> 32) & 0x07)));
        out.push_back(static_cast<std::byte>((uval >> 24) & 0xFF));
        out.push_back(static_cast<std::byte>((uval >> 16) & 0xFF));
        out.push_back(static_cast<std::byte>((uval >> 8) & 0xFF));
        out.push_back(static_cast<std::byte>(uval & 0xFF));
        return 5;
    }
    if (uval < 0x40000000000ULL) {
        out.push_back(static_cast<std::byte>(0xF8 | ((uval >> 40) & 0x03)));
        out.push_back(static_cast<std::byte>((uval >> 32) & 0xFF));
        out.push_back(static_cast<std::byte>((uval >> 24) & 0xFF));
        out.push_back(static_cast<std::byte>((uval >> 16) & 0xFF));
        out.push_back(static_cast<std::byte>((uval >> 8) & 0xFF));
        out.push_back(static_cast<std::byte>(uval & 0xFF));
        return 6;
    }
    if (uval < 0x2000000000000ULL) {
        out.push_back(static_cast<std::byte>(0xFC | ((uval >> 48) & 0x01)));
        out.push_back(static_cast<std::byte>((uval >> 40) & 0xFF));
        out.push_back(static_cast<std::byte>((uval >> 32) & 0xFF));
        out.push_back(static_cast<std::byte>((uval >> 24) & 0xFF));
        out.push_back(static_cast<std::byte>((uval >> 16) & 0xFF));
        out.push_back(static_cast<std::byte>((uval >> 8) & 0xFF));
        out.push_back(static_cast<std::byte>(uval & 0xFF));
        return 7;
    }
    if (uval < 0x100000000000000ULL) {
        out.push_back(static_cast<std::byte>(0xFE));
        out.push_back(static_cast<std::byte>((uval >> 48) & 0xFF));
        out.push_back(static_cast<std::byte>((uval >> 40) & 0xFF));
        out.push_back(static_cast<std::byte>((uval >> 32) & 0xFF));
        out.push_back(static_cast<std::byte>((uval >> 24) & 0xFF));
        out.push_back(static_cast<std::byte>((uval >> 16) & 0xFF));
        out.push_back(static_cast<std::byte>((uval >> 8) & 0xFF));
        out.push_back(static_cast<std::byte>(uval & 0xFF));
        return 8;
    }
    out.push_back(static_cast<std::byte>(0xFF));
    for (int i = 56; i >= 0; i -= 8) {
        out.push_back(static_cast<std::byte>((uval >> i) & 0xFF));
    }
    return 9;
}

// ---------------------------------------------------------------------------
// File Definition
// ---------------------------------------------------------------------------

CramFileDefinition ParseFileDefinition(std::span<const std::byte> data)
{
    if (std::size(data) < 26) {
        throw std::runtime_error("CRAM file definition too short");
    }
    if (!std::ranges::equal(CRAM_MAGIC, data.first(std::size(CRAM_MAGIC)))) {
        throw std::runtime_error("not a CRAM file: bad magic");
    }

    CramFileDefinition def;
    def.MajorVersion = static_cast<std::uint8_t>(data[4]);
    def.MinorVersion = static_cast<std::uint8_t>(data[5]);
    std::copy_n(data.data() + 6, 20, def.FileId.data());
    return def;
}

std::vector<std::byte> SerializeFileDefinition(const CramFileDefinition& def)
{
    std::vector<std::byte> out;
    out.reserve(26);
    out.insert(std::end(out), std::begin(CRAM_MAGIC), std::end(CRAM_MAGIC));
    out.push_back(static_cast<std::byte>(def.MajorVersion));
    out.push_back(static_cast<std::byte>(def.MinorVersion));
    out.insert(std::end(out), std::begin(def.FileId), std::end(def.FileId));
    return out;
}

// ---------------------------------------------------------------------------
// Container Header
// ---------------------------------------------------------------------------

CramContainerHeader ParseContainerHeader(std::span<const std::byte> data, std::size_t& bytesRead)
{
    if (std::size(data) < 4) {
        throw std::runtime_error("container header too short");
    }

    CramContainerHeader hdr;
    std::size_t pos = 0;

    hdr.Length = ReadI32LE(data.data() + pos);
    pos += 4;

    std::size_t n = 0;
    hdr.RefSeqId = ReadItf8(data.subspan(pos), n);
    pos += n;
    hdr.StartPos = ReadItf8(data.subspan(pos), n);
    pos += n;
    hdr.AlignmentSpan = ReadItf8(data.subspan(pos), n);
    pos += n;
    hdr.NumRecords = ReadItf8(data.subspan(pos), n);
    pos += n;
    hdr.RecordCounter = ReadLtf8(data.subspan(pos), n);
    pos += n;
    hdr.Bases = ReadLtf8(data.subspan(pos), n);
    pos += n;
    hdr.NumBlocks = ReadItf8(data.subspan(pos), n);
    pos += n;
    if (hdr.NumBlocks < 0) {
        throw std::runtime_error("container header: negative num_blocks");
    }

    // landmarks array
    const auto landmarkCount = ReadItf8(data.subspan(pos), n);
    pos += n;
    if (landmarkCount < 0) {
        throw std::runtime_error("container header: negative landmark count");
    }
    hdr.Landmarks.resize(landmarkCount);
    for (std::int32_t i = 0; i < landmarkCount; ++i) {
        hdr.Landmarks[i] = ReadItf8(data.subspan(pos), n);
        pos += n;
    }

    // CRC32
    if (pos + 4 > std::size(data)) {
        throw std::runtime_error("container header: truncated CRC32");
    }
    hdr.Crc32 = ReadU32LE(data.data() + pos);
    pos += 4;

    const auto computedCrc = ComputeCrc32(data.subspan(0, pos - 4));
    if (computedCrc != hdr.Crc32) {
        throw std::runtime_error("container header CRC32 mismatch");
    }

    bytesRead = pos;
    return hdr;
}

std::vector<std::byte> SerializeContainerHeader(const CramContainerHeader& header)
{
    std::vector<std::byte> out;
    out.reserve(64);

    WriteI32LE(out, header.Length);
    WriteItf8(out, header.RefSeqId);
    WriteItf8(out, header.StartPos);
    WriteItf8(out, header.AlignmentSpan);
    WriteItf8(out, header.NumRecords);
    WriteLtf8(out, header.RecordCounter);
    WriteLtf8(out, header.Bases);
    WriteItf8(out, header.NumBlocks);

    WriteItf8(out, static_cast<std::int32_t>(std::size(header.Landmarks)));
    for (const auto lm : header.Landmarks) {
        WriteItf8(out, lm);
    }

    // CRC32 of everything so far
    const auto crc = ComputeCrc32(out);
    WriteU32LE(out, crc);

    return out;
}

// ---------------------------------------------------------------------------
// Block
// ---------------------------------------------------------------------------

CramBlock ParseBlock(std::span<const std::byte> data, std::size_t& bytesRead)
{
    if (std::size(data) < 5) {
        throw std::runtime_error("block too short");
    }

    CramBlock block;
    std::size_t pos = 0;

    const auto methodByte = static_cast<std::uint8_t>(data[pos]);
    ++pos;
    block.Method = static_cast<CramBlockMethod>(methodByte);
    const auto contentTypeByte = static_cast<std::uint8_t>(data[pos]);
    ++pos;
    block.ContentType = static_cast<CramBlockContentType>(contentTypeByte);

    std::size_t n = 0;
    block.ContentId = ReadItf8(data.subspan(pos), n);
    pos += n;
    block.CompressedSize = ReadItf8(data.subspan(pos), n);
    pos += n;
    block.RawSize = ReadItf8(data.subspan(pos), n);
    pos += n;
    if (block.CompressedSize < 0 || block.RawSize < 0) {
        throw std::runtime_error("block: negative compressed/raw size");
    }

    if (std::size(data) - pos < 4) {
        throw std::runtime_error("block: truncated CRC32");
    }
    const auto payloadMax = std::size(data) - pos - 4;

    // Block data (compressed)
    if (static_cast<std::size_t>(block.CompressedSize) > payloadMax) {
        throw std::runtime_error("block: truncated data");
    }
    block.Data.assign(data.data() + pos, data.data() + pos + block.CompressedSize);
    pos += block.CompressedSize;

    // CRC32
    block.Crc32 = ReadU32LE(data.data() + pos);
    pos += 4;

    const auto computedCrc = ComputeCrc32(data.subspan(0, pos - 4));
    if (computedCrc != block.Crc32) {
        throw std::runtime_error("block CRC32 mismatch");
    }

    bytesRead = pos;
    return block;
}

std::vector<std::byte> SerializeBlock(const CramBlock& block)
{
    const auto dataSize = std::size(block.Data);
    if (dataSize > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())) {
        throw std::runtime_error("block: compressed size exceeds ITF8 range");
    }

    std::vector<std::byte> out;
    out.reserve(32 + dataSize);
    out.resize(2);
    out[0] = static_cast<std::byte>(std::to_underlying(block.Method));
    out[1] = static_cast<std::byte>(std::to_underlying(block.ContentType));
    WriteItf8(out, block.ContentId);
    WriteItf8(out, static_cast<std::int32_t>(dataSize));  // compressed size
    WriteItf8(out, block.RawSize);

    out.insert(std::end(out), std::begin(block.Data), std::end(block.Data));

    // CRC32 of everything in this block
    const auto crc = ComputeCrc32(out);
    WriteU32LE(out, crc);

    return out;
}

// ---------------------------------------------------------------------------
// Slice Header
// ---------------------------------------------------------------------------

CramSliceHeader ParseSliceHeader(std::span<const std::byte> data)
{
    CramSliceHeader hdr;
    std::size_t pos = 0;
    std::size_t n = 0;

    hdr.RefSeqId = ReadItf8(data.subspan(pos), n);
    pos += n;
    hdr.AlignmentStart = ReadItf8(data.subspan(pos), n);
    pos += n;
    hdr.AlignmentSpan = ReadItf8(data.subspan(pos), n);
    pos += n;
    hdr.NumRecords = ReadItf8(data.subspan(pos), n);
    pos += n;
    hdr.RecordCounter = ReadLtf8(data.subspan(pos), n);
    pos += n;
    hdr.NumBlocks = ReadItf8(data.subspan(pos), n);
    pos += n;

    const auto numContentIds = ReadItf8(data.subspan(pos), n);
    pos += n;
    hdr.BlockContentIds.resize(numContentIds);
    for (std::int32_t i = 0; i < numContentIds; ++i) {
        hdr.BlockContentIds[i] = ReadItf8(data.subspan(pos), n);
        pos += n;
    }

    hdr.EmbeddedRefBlockId = ReadItf8(data.subspan(pos), n);
    pos += n;

    // MD5 (16 bytes)
    if (pos + 16 > std::size(data)) {
        throw std::runtime_error("slice header: truncated MD5");
    }
    std::copy_n(data.data() + pos, 16, hdr.RefMd5.data());
    pos += 16;

    if (pos < std::size(data)) {
        hdr.OptionalTags.assign(data.begin() + static_cast<std::ptrdiff_t>(pos), data.end());
    }

    return hdr;
}

std::vector<std::byte> SerializeSliceHeader(const CramSliceHeader& header)
{
    std::vector<std::byte> out;
    out.reserve(64);

    WriteItf8(out, header.RefSeqId);
    WriteItf8(out, header.AlignmentStart);
    WriteItf8(out, header.AlignmentSpan);
    WriteItf8(out, header.NumRecords);
    WriteLtf8(out, header.RecordCounter);
    WriteItf8(out, header.NumBlocks);

    WriteItf8(out, static_cast<std::int32_t>(std::size(header.BlockContentIds)));
    for (const auto id : header.BlockContentIds) {
        WriteItf8(out, id);
    }

    WriteItf8(out, header.EmbeddedRefBlockId);

    out.insert(std::end(out), std::begin(header.RefMd5), std::end(header.RefMd5));
    out.insert(std::end(out), std::begin(header.OptionalTags), std::end(header.OptionalTags));

    return out;
}

// ---------------------------------------------------------------------------
// Compression Header
// ---------------------------------------------------------------------------

namespace {

CramEncodingDescriptor ParseEncodingDescriptor(std::span<const std::byte> data, std::size_t& pos,
                                               std::size_t endPos)
{
    if (pos > std::size(data) || endPos > std::size(data) || pos > endPos) {
        throw std::runtime_error("encoding descriptor: invalid bounds");
    }
    CramEncodingDescriptor desc;
    std::size_t n = 0;

    if (pos >= endPos) {
        throw std::runtime_error("encoding descriptor: truncated codec id");
    }
    desc.CodecId = static_cast<CramCodecId>(ReadItf8(data.subspan(pos), n));
    pos += n;
    if (pos > endPos) {
        throw std::runtime_error("encoding descriptor: truncated codec id");
    }

    if (pos >= endPos) {
        throw std::runtime_error("encoding descriptor: truncated parameter length");
    }
    const auto paramLen = ReadItf8(data.subspan(pos), n);
    pos += n;
    if (paramLen < 0) {
        throw std::runtime_error("encoding descriptor: negative parameter length");
    }
    if (pos > endPos || static_cast<std::size_t>(paramLen) > (endPos - pos)) {
        throw std::runtime_error("encoding descriptor: truncated parameters");
    }

    desc.Parameters.assign(data.data() + pos, data.data() + pos + paramLen);
    pos += paramLen;

    return desc;
}

void SerializeEncodingDescriptor(std::vector<std::byte>& out, const CramEncodingDescriptor& desc)
{
    WriteItf8(out, std::to_underlying(desc.CodecId));
    WriteItf8(out, static_cast<std::int32_t>(std::size(desc.Parameters)));
    out.insert(std::end(out), std::begin(desc.Parameters), std::end(desc.Parameters));
}

}  // namespace

CramCompressionHeader ParseCompressionHeader(std::span<const std::byte> data)
{
    CramCompressionHeader header;
    std::size_t pos = 0;
    std::size_t n = 0;
    const auto nextByte = [&data, &pos]() {
        const std::byte value{data[pos]};
        ++pos;
        return value;
    };

    // --- Preservation map ---
    if (std::empty(data)) {
        throw std::runtime_error("compression header: truncated preservation map size");
    }
    const auto pmapSize = ReadItf8(data.subspan(pos), n);
    pos += n;
    if (pmapSize < 0 || static_cast<std::size_t>(pmapSize) > (std::size(data) - pos)) {
        throw std::runtime_error("compression header: invalid preservation map size");
    }
    const auto pmapEnd = pos + pmapSize;

    if (pos >= pmapEnd) {
        throw std::runtime_error("compression header: truncated preservation map count");
    }
    const auto pmapCount = ReadItf8(data.subspan(pos), n);
    pos += n;
    if (pmapCount < 0) {
        throw std::runtime_error("compression header: negative preservation map count");
    }

    for (std::int32_t i = 0; i < pmapCount; ++i) {
        if (pos + 2 > pmapEnd) {
            throw std::runtime_error("compression header: truncated preservation map");
        }
        const char k0 = static_cast<char>(nextByte());
        const char k1 = static_cast<char>(nextByte());

        if (k0 == 'R' && k1 == 'N') {
            if (pos >= pmapEnd) {
                throw std::runtime_error("compression header: truncated RN value");
            }
            header.PreservationMap.ReadNamesIncluded = (static_cast<std::uint8_t>(nextByte()) != 0);
        } else if (k0 == 'A' && k1 == 'P') {
            if (pos >= pmapEnd) {
                throw std::runtime_error("compression header: truncated AP value");
            }
            header.PreservationMap.ApDelta = (static_cast<std::uint8_t>(nextByte()) != 0);
        } else if (k0 == 'R' && k1 == 'R') {
            if (pos >= pmapEnd) {
                throw std::runtime_error("compression header: truncated RR value");
            }
            header.PreservationMap.ReferenceRequired = (static_cast<std::uint8_t>(nextByte()) != 0);
        } else if (k0 == 'S' && k1 == 'M') {
            if (pmapEnd - pos < 5) {
                throw std::runtime_error("compression header: truncated SM value");
            }
            std::copy_n(data.data() + pos, 5, header.PreservationMap.SubstitutionMatrix.data());
            pos += 5;
        } else if (k0 == 'T' && k1 == 'D') {
            if (pos >= pmapEnd) {
                throw std::runtime_error("compression header: truncated TD length");
            }
            const auto tdLen = ReadItf8(data.subspan(pos), n);
            pos += n;
            if (tdLen < 0 || static_cast<std::size_t>(tdLen) > (pmapEnd - pos)) {
                throw std::runtime_error("compression header: truncated TD payload");
            }
            header.PreservationMap.TagIdsDictionary.assign(data.data() + pos,
                                                           data.data() + pos + tdLen);
            pos += tdLen;
        }
    }
    if (pos > pmapEnd) {
        throw std::runtime_error("compression header: preservation map overrun");
    }
    pos = pmapEnd;

    // --- Data series encoding map ---
    if (pos >= std::size(data)) {
        throw std::runtime_error("compression header: truncated data series map size");
    }
    const auto dsMapSize = ReadItf8(data.subspan(pos), n);
    pos += n;
    if (dsMapSize < 0 || static_cast<std::size_t>(dsMapSize) > (std::size(data) - pos)) {
        throw std::runtime_error("compression header: invalid data series map size");
    }
    const auto dsMapEnd = pos + dsMapSize;

    if (pos >= dsMapEnd) {
        throw std::runtime_error("compression header: truncated data series count");
    }
    const auto dsCount = ReadItf8(data.subspan(pos), n);
    pos += n;
    if (dsCount < 0) {
        throw std::runtime_error("compression header: negative data series count");
    }

    for (std::int32_t i = 0; i < dsCount; ++i) {
        if (pos + 2 > dsMapEnd) {
            throw std::runtime_error("compression header: truncated data series map");
        }
        const auto c0 = static_cast<std::uint8_t>(nextByte());
        const auto c1 = static_cast<std::uint8_t>(nextByte());
        const auto key = static_cast<CramDataSeries>((c0 << 8) | c1);

        auto desc = ParseEncodingDescriptor(data, pos, dsMapEnd);
        header.DataSeriesEncodings.emplace_back(key, std::move(desc));
    }
    if (pos > dsMapEnd) {
        throw std::runtime_error("compression header: data series map overrun");
    }
    pos = dsMapEnd;

    // --- Tag encoding map ---
    if (pos >= std::size(data)) {
        throw std::runtime_error("compression header: truncated tag map size");
    }
    const auto tagMapSize = ReadItf8(data.subspan(pos), n);
    pos += n;
    if (tagMapSize < 0 || static_cast<std::size_t>(tagMapSize) > (std::size(data) - pos)) {
        throw std::runtime_error("compression header: invalid tag map size");
    }
    const auto tagMapEnd = pos + tagMapSize;

    if (pos >= tagMapEnd) {
        throw std::runtime_error("compression header: truncated tag count");
    }
    const auto tagCount = ReadItf8(data.subspan(pos), n);
    pos += n;
    if (tagCount < 0) {
        throw std::runtime_error("compression header: negative tag count");
    }

    for (std::int32_t i = 0; i < tagCount; ++i) {
        if (pos >= tagMapEnd) {
            throw std::runtime_error("compression header: truncated tag key");
        }
        const auto tagKey = ReadItf8(data.subspan(pos), n);
        pos += n;
        auto desc = ParseEncodingDescriptor(data, pos, tagMapEnd);
        header.TagEncodings.emplace_back(tagKey, std::move(desc));
    }
    if (pos > tagMapEnd) {
        throw std::runtime_error("compression header: tag map overrun");
    }
    pos = tagMapEnd;

    return header;
}

std::vector<std::byte> SerializeCompressionHeader(const CramCompressionHeader& header)
{
    std::vector<std::byte> out;
    out.reserve(256);

    // --- Preservation map ---
    std::vector<std::byte> pmapBody;
    std::int32_t pmapCount = 0;

    // RN
    pmapBody.push_back(static_cast<std::byte>('R'));
    pmapBody.push_back(static_cast<std::byte>('N'));
    pmapBody.push_back(header.PreservationMap.ReadNamesIncluded ? std::byte{1} : std::byte{0});
    ++pmapCount;

    // AP
    pmapBody.push_back(static_cast<std::byte>('A'));
    pmapBody.push_back(static_cast<std::byte>('P'));
    pmapBody.push_back(header.PreservationMap.ApDelta ? std::byte{1} : std::byte{0});
    ++pmapCount;

    // RR
    pmapBody.push_back(static_cast<std::byte>('R'));
    pmapBody.push_back(static_cast<std::byte>('R'));
    pmapBody.push_back(header.PreservationMap.ReferenceRequired ? std::byte{1} : std::byte{0});
    ++pmapCount;

    // SM
    pmapBody.push_back(static_cast<std::byte>('S'));
    pmapBody.push_back(static_cast<std::byte>('M'));
    pmapBody.insert(std::end(pmapBody), std::begin(header.PreservationMap.SubstitutionMatrix),
                    std::end(header.PreservationMap.SubstitutionMatrix));
    ++pmapCount;

    // TD
    pmapBody.push_back(static_cast<std::byte>('T'));
    pmapBody.push_back(static_cast<std::byte>('D'));
    WriteItf8(pmapBody,
              static_cast<std::int32_t>(std::size(header.PreservationMap.TagIdsDictionary)));
    pmapBody.insert(std::end(pmapBody), std::begin(header.PreservationMap.TagIdsDictionary),
                    std::end(header.PreservationMap.TagIdsDictionary));
    ++pmapCount;

    // Write pmap: size + count + body
    std::vector<std::byte> pmapSized;
    WriteItf8(pmapSized, pmapCount);
    pmapSized.insert(std::end(pmapSized), std::begin(pmapBody), std::end(pmapBody));
    WriteItf8(out, static_cast<std::int32_t>(std::size(pmapSized)));
    out.insert(std::end(out), std::begin(pmapSized), std::end(pmapSized));

    // --- Data series encoding map ---
    std::vector<std::byte> dsBody;
    WriteItf8(dsBody, static_cast<std::int32_t>(std::size(header.DataSeriesEncodings)));
    for (const auto& [key, desc] : header.DataSeriesEncodings) {
        const auto k = std::to_underlying(key);
        dsBody.push_back(static_cast<std::byte>((k >> 8) & 0xFF));
        dsBody.push_back(static_cast<std::byte>(k & 0xFF));
        SerializeEncodingDescriptor(dsBody, desc);
    }
    WriteItf8(out, static_cast<std::int32_t>(std::size(dsBody)));
    out.insert(std::end(out), std::begin(dsBody), std::end(dsBody));

    // --- Tag encoding map ---
    std::vector<std::byte> tagBody;
    WriteItf8(tagBody, static_cast<std::int32_t>(std::size(header.TagEncodings)));
    for (const auto& [key, desc] : header.TagEncodings) {
        WriteItf8(tagBody, key);
        SerializeEncodingDescriptor(tagBody, desc);
    }
    WriteItf8(out, static_cast<std::int32_t>(std::size(tagBody)));
    out.insert(std::end(out), std::begin(tagBody), std::end(tagBody));

    return out;
}

std::vector<std::byte> SerializeSlice(const CramSlice& slice)
{
    const auto sliceHdrData = SerializeSliceHeader(slice.Header);
    CramBlock sliceHdrBlock;
    sliceHdrBlock.Method = CramBlockMethod::RAW;
    sliceHdrBlock.ContentType = CramBlockContentType::SLICE_HEADER;
    sliceHdrBlock.ContentId = 0;
    sliceHdrBlock.RawSize = static_cast<std::int32_t>(std::size(sliceHdrData));
    sliceHdrBlock.CompressedSize = sliceHdrBlock.RawSize;
    sliceHdrBlock.Data = sliceHdrData;

    std::vector<std::byte> out;
    out.reserve(SerializedSliceSize(slice));

    const auto hdrBytes = SerializeBlock(sliceHdrBlock);
    out.insert(std::end(out), std::begin(hdrBytes), std::end(hdrBytes));

    const auto coreBytes = SerializeBlock(slice.CoreBlock);
    out.insert(std::end(out), std::begin(coreBytes), std::end(coreBytes));

    for (const auto& extBlock : slice.ExternalBlocks) {
        const auto extBytes = SerializeBlock(extBlock);
        out.insert(std::end(out), std::begin(extBytes), std::end(extBytes));
    }

    return out;
}

std::size_t SerializedSliceSize(const CramSlice& slice)
{
    const auto sliceHdrData = SerializeSliceHeader(slice.Header);
    CramBlock sliceHdrBlock;
    sliceHdrBlock.Method = CramBlockMethod::RAW;
    sliceHdrBlock.ContentType = CramBlockContentType::SLICE_HEADER;
    sliceHdrBlock.ContentId = 0;
    sliceHdrBlock.RawSize = static_cast<std::int32_t>(std::size(sliceHdrData));
    sliceHdrBlock.CompressedSize = sliceHdrBlock.RawSize;
    sliceHdrBlock.Data = std::move(sliceHdrData);

    std::size_t total = SerializedBlockSize(sliceHdrBlock);
    total += SerializedBlockSize(slice.CoreBlock);
    for (const auto& extBlock : slice.ExternalBlocks) {
        total += SerializedBlockSize(extBlock);
    }
    return total;
}

std::vector<std::byte> SerializeContainer(const CramContainer& container)
{
    const auto compHdrData = SerializeCompressionHeader(container.CompressionHeader);
    CramBlock compHdrBlock;
    compHdrBlock.Method = CramBlockMethod::RAW;
    compHdrBlock.ContentType = CramBlockContentType::COMPRESSION_HEADER;
    compHdrBlock.ContentId = 0;
    compHdrBlock.RawSize = static_cast<std::int32_t>(std::size(compHdrData));
    compHdrBlock.CompressedSize = compHdrBlock.RawSize;
    compHdrBlock.Data = compHdrData;
    const auto compHdrBytes = SerializeBlock(compHdrBlock);

    CramContainerHeader containerHeader = container.Header;
    containerHeader.Landmarks.clear();
    containerHeader.Landmarks.reserve(std::size(container.Slices));

    std::int64_t payloadSize = static_cast<std::int64_t>(std::size(compHdrBytes));
    for (const auto& slice : container.Slices) {
        containerHeader.Landmarks.push_back(static_cast<std::int32_t>(payloadSize));
        const auto sliceSize = static_cast<std::int64_t>(SerializedSliceSize(slice));
        ValidateAndAdvance(payloadSize, sliceSize, "serializing container slices");
    }
    containerHeader.Length = static_cast<std::int32_t>(payloadSize);

    std::int64_t totalBlocks = 1;  // compression header block
    for (const auto& slice : container.Slices) {
        ValidateAndAdvance(totalBlocks, 2, "counting slice header/core blocks");
        ValidateAndAdvance(totalBlocks, static_cast<std::int64_t>(std::size(slice.ExternalBlocks)),
                           "counting external blocks");
    }
    containerHeader.NumBlocks = static_cast<std::int32_t>(totalBlocks);

    const auto containerHeaderBytes = SerializeContainerHeader(containerHeader);
    std::vector<std::byte> out;
    out.reserve(std::size(containerHeaderBytes) + static_cast<std::size_t>(payloadSize));
    out.insert(std::end(out), std::begin(containerHeaderBytes), std::end(containerHeaderBytes));
    out.insert(std::end(out), std::begin(compHdrBytes), std::end(compHdrBytes));
    for (const auto& slice : container.Slices) {
        const auto sliceBytes = SerializeSlice(slice);
        out.insert(std::end(out), std::begin(sliceBytes), std::end(sliceBytes));
    }
    return out;
}

CramContainer ParseContainer(const CramContainerHeader& header, std::span<const std::byte> payload)
{
    if (header.Length < 0) {
        throw std::runtime_error("ParseContainer: negative container length");
    }
    if (static_cast<std::size_t>(header.Length) != std::size(payload)) {
        throw std::runtime_error("ParseContainer: payload size does not match container header");
    }
    if (header.NumBlocks < 0) {
        throw std::runtime_error("ParseContainer: negative container block count");
    }

    CramContainer container;
    container.Header = header;

    if (header.NumBlocks == 0) {
        if (!std::empty(payload)) {
            throw std::runtime_error("ParseContainer: non-empty payload with zero block count");
        }
        return container;
    }

    std::size_t pos = 0;
    std::size_t bytesRead = 0;
    std::size_t blocksParsed = 0;

    CramBlock compressionHeaderBlock =
        ParseBlock(std::span<const std::byte>{payload}.subspan(pos), bytesRead);
    pos += bytesRead;
    ++blocksParsed;
    if (compressionHeaderBlock.ContentType != CramBlockContentType::COMPRESSION_HEADER) {
        throw std::runtime_error("ParseContainer: first block is not COMPRESSION_HEADER");
    }
    DecompressCramBlock(compressionHeaderBlock);
    container.CompressionHeader = ParseCompressionHeader(compressionHeaderBlock.Data);

    std::vector<std::int32_t> observedLandmarks;
    observedLandmarks.reserve(std::size(header.Landmarks));

    while (blocksParsed < static_cast<std::size_t>(header.NumBlocks)) {
        if (pos >= std::size(payload)) {
            throw std::runtime_error("ParseContainer: truncated slice list");
        }

        observedLandmarks.push_back(static_cast<std::int32_t>(pos));

        CramBlock sliceHeaderBlock =
            ParseBlock(std::span<const std::byte>{payload}.subspan(pos), bytesRead);
        pos += bytesRead;
        ++blocksParsed;
        if (sliceHeaderBlock.ContentType != CramBlockContentType::SLICE_HEADER) {
            throw std::runtime_error("ParseContainer: expected SLICE_HEADER block");
        }

        DecompressCramBlock(sliceHeaderBlock);
        CramSlice slice;
        slice.Header = ParseSliceHeader(sliceHeaderBlock.Data);

        if (slice.Header.NumBlocks < 0) {
            throw std::runtime_error("ParseContainer: negative slice block count");
        }
        if (blocksParsed + static_cast<std::size_t>(slice.Header.NumBlocks) >
            static_cast<std::size_t>(header.NumBlocks)) {
            throw std::runtime_error("ParseContainer: slice block count exceeds container payload");
        }

        bool sawCore = false;
        for (std::int32_t i = 0; i < slice.Header.NumBlocks; ++i) {
            CramBlock block =
                ParseBlock(std::span<const std::byte>{payload}.subspan(pos), bytesRead);
            pos += bytesRead;
            ++blocksParsed;

            switch (block.ContentType) {
                case CramBlockContentType::CORE_DATA:
                    if (sawCore) {
                        throw std::runtime_error(
                            "ParseContainer: multiple CORE_DATA blocks in slice");
                    }
                    slice.CoreBlock = std::move(block);
                    sawCore = true;
                    break;
                case CramBlockContentType::EXTERNAL_DATA:
                    slice.ExternalBlocks.push_back(std::move(block));
                    break;
                default:
                    throw std::runtime_error(
                        "ParseContainer: unsupported slice data block content type");
            }
        }

        if (!sawCore) {
            throw std::runtime_error("ParseContainer: slice missing CORE_DATA block");
        }

        container.Slices.push_back(std::move(slice));
    }

    if (pos != std::size(payload)) {
        throw std::runtime_error("ParseContainer: trailing bytes after declared block list");
    }

    if (!std::empty(header.Landmarks)) {
        if (std::size(header.Landmarks) != std::size(observedLandmarks)) {
            throw std::runtime_error("ParseContainer: landmark count mismatch");
        }
        for (std::size_t i = 0; i < std::size(observedLandmarks); ++i) {
            if (header.Landmarks[i] != observedLandmarks[i]) {
                throw std::runtime_error("ParseContainer: landmark offset mismatch");
            }
        }
    }

    return container;
}

// ---------------------------------------------------------------------------
// EOF marker
// ---------------------------------------------------------------------------

bool IsCramEofMarker(std::span<const std::byte> data)
{
    if (std::size(data) < std::size(CRAM_EOF_MARKER)) {
        return false;
    }
    return std::ranges::equal(CRAM_EOF_MARKER, data.first(std::size(CRAM_EOF_MARKER)));
}

}  // namespace Samoa
}  // namespace PacBio
