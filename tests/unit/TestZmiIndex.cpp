#include "../../src/tools/zmi-index/RecordParser.hpp"

#include <pbsamoa/core/Bgzf.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <gtest/gtest.h>

#include <array>
#include <cstdint>
#include <span>
#include <vector>

#include <cstddef>
#include <cstring>

namespace PacBio {
namespace Samoa {
namespace ZmiIndex {
namespace {

// Helpers for building synthetic BAM bytes in tests.

void AppendU32LE(std::vector<std::byte>& out, std::uint32_t value)
{
    for (int i{0}; i < 4; ++i) {
        out.push_back(static_cast<std::byte>((value >> (i * 8)) & 0xFFU));
    }
}

void AppendI32LE(std::vector<std::byte>& out, std::int32_t value)
{
    AppendU32LE(out, static_cast<std::uint32_t>(value));
}

void AppendU16LE(std::vector<std::byte>& out, std::uint16_t value)
{
    out.push_back(static_cast<std::byte>(value & 0xFFU));
    out.push_back(static_cast<std::byte>((value >> 8U) & 0xFFU));
}

void AppendString(std::vector<std::byte>& out, std::string_view text, bool nulTerminate)
{
    for (const char c : text) {
        out.push_back(static_cast<std::byte>(c));
    }
    if (nulTerminate) {
        out.push_back(std::byte{0});
    }
}

void AppendRgZ(std::vector<std::byte>& aux, std::string_view value)
{
    aux.push_back(static_cast<std::byte>('R'));
    aux.push_back(static_cast<std::byte>('G'));
    aux.push_back(static_cast<std::byte>('Z'));
    AppendString(aux, value, true);
}

void AppendXiInt(std::vector<std::byte>& aux, char tag1, char tag2, std::int32_t value)
{
    aux.push_back(static_cast<std::byte>(tag1));
    aux.push_back(static_cast<std::byte>(tag2));
    aux.push_back(static_cast<std::byte>('i'));
    AppendI32LE(aux, value);
}

// Build a minimal BAM record body (bytes after the 4-byte block_size prefix).
// l_seq=0, no CIGAR, no qual. Aux is provided externally.
std::vector<std::byte> BuildRecordBody(std::string_view name, std::span<const std::byte> aux)
{
    std::vector<std::byte> body;

    // Fixed fields (32 bytes), zeroed except l_read_name
    body.resize(32U, std::byte{0});
    const std::uint8_t nameLen{static_cast<std::uint8_t>(name.size() + 1U)};  // includes NUL
    body[8] = static_cast<std::byte>(nameLen);
    // refID=-1, pos=-1, nextRefID=-1, nextPos=-1
    for (const std::size_t off : {0U, 4U, 20U, 24U}) {
        for (std::size_t i{0}; i < 4U; ++i) {
            body[off + i] = std::byte{0xFFU};
        }
    }

    // Read name (NUL-terminated)
    AppendString(body, name, true);

    // No CIGAR, no seq, no qual (l_seq=0, n_cigar_op=0)
    // Aux
    body.insert(body.end(), aux.begin(), aux.end());
    return body;
}

// Wrap a record body with the 4-byte block_size prefix.
std::vector<std::byte> WithBlockSize(std::span<const std::byte> body)
{
    std::vector<std::byte> result;
    AppendU32LE(result, static_cast<std::uint32_t>(body.size()));
    result.insert(result.end(), body.begin(), body.end());
    return result;
}

// Build a minimal BAM header bytes (decompressed): "BAM\1", l_text=0, n_ref=0.
std::vector<std::byte> BuildEmptyBamHeader()
{
    std::vector<std::byte> header;
    header.push_back(static_cast<std::byte>('B'));
    header.push_back(static_cast<std::byte>('A'));
    header.push_back(static_cast<std::byte>('M'));
    header.push_back(static_cast<std::byte>(0x01));
    AppendU32LE(header, 0U);  // l_text
    AppendI32LE(header, 0);   // n_ref
    return header;
}

}  // namespace

// --- AuxValueLength ---

TEST(ZmiIndexAux, FixedSizeTypes)
{
    using detail::AuxValueLength;
    std::array<std::byte, 8> dummy{};
    EXPECT_EQ(AuxValueLength(std::byte{'A'}, dummy), 1U);
    EXPECT_EQ(AuxValueLength(std::byte{'c'}, dummy), 1U);
    EXPECT_EQ(AuxValueLength(std::byte{'C'}, dummy), 1U);
    EXPECT_EQ(AuxValueLength(std::byte{'s'}, dummy), 2U);
    EXPECT_EQ(AuxValueLength(std::byte{'S'}, dummy), 2U);
    EXPECT_EQ(AuxValueLength(std::byte{'i'}, dummy), 4U);
    EXPECT_EQ(AuxValueLength(std::byte{'I'}, dummy), 4U);
    EXPECT_EQ(AuxValueLength(std::byte{'f'}, dummy), 4U);
    EXPECT_EQ(AuxValueLength(std::byte{'d'}, dummy), 8U);
}

TEST(ZmiIndexAux, ZTerminatedString)
{
    using detail::AuxValueLength;
    std::vector<std::byte> data;
    AppendString(data, "abc", true);
    EXPECT_EQ(AuxValueLength(std::byte{'Z'}, data), 4U);
}

TEST(ZmiIndexAux, ZUnterminated)
{
    using detail::AuxValueLength;
    std::vector<std::byte> data;
    AppendString(data, "abc", false);
    EXPECT_FALSE(AuxValueLength(std::byte{'Z'}, data).has_value());
}

TEST(ZmiIndexAux, BArray)
{
    using detail::AuxValueLength;
    std::vector<std::byte> data;
    data.push_back(static_cast<std::byte>('S'));  // subtype uint16
    AppendU32LE(data, 3U);                        // count
    AppendU16LE(data, 10U);
    AppendU16LE(data, 20U);
    AppendU16LE(data, 30U);
    EXPECT_EQ(AuxValueLength(std::byte{'B'}, data), 1U + 4U + 3U * 2U);
}

TEST(ZmiIndexAux, UnknownTypeRejected)
{
    using detail::AuxValueLength;
    std::array<std::byte, 4> dummy{};
    EXPECT_FALSE(AuxValueLength(std::byte{'q'}, dummy).has_value());
}

// --- FindZTag ---

TEST(ZmiIndexAux, FindRgPresent)
{
    std::vector<std::byte> aux;
    AppendXiInt(aux, 'x', 'y', 7);
    AppendRgZ(aux, "deadbeef");
    AppendXiInt(aux, 'a', 'b', 9);
    const auto value{detail::FindZTag(aux, RG_TAG)};
    ASSERT_TRUE(value.has_value());
    EXPECT_EQ(*value, "deadbeef");
}

TEST(ZmiIndexAux, FindRgAbsent)
{
    std::vector<std::byte> aux;
    AppendXiInt(aux, 'x', 'y', 7);
    EXPECT_FALSE(detail::FindZTag(aux, RG_TAG).has_value());
}

TEST(ZmiIndexAux, FindRgMalformedNoNul)
{
    std::vector<std::byte> aux;
    aux.push_back(static_cast<std::byte>('R'));
    aux.push_back(static_cast<std::byte>('G'));
    aux.push_back(static_cast<std::byte>('Z'));
    AppendString(aux, "noterm", false);
    EXPECT_FALSE(detail::FindZTag(aux, RG_TAG).has_value());
}

// --- ParseRecordIdentity ---

TEST(ZmiIndexRecord, ParseIdentityWithRg)
{
    std::vector<std::byte> aux;
    AppendRgZ(aux, "deadbeef");
    const std::vector<std::byte> body{BuildRecordBody("movie/123/ccs", aux)};
    const auto [rgId, zmw]{detail::ParseRecordIdentity(body)};
    EXPECT_EQ(zmw, 123);
    EXPECT_EQ(rgId, static_cast<std::int32_t>(0xDEADBEEFU));
}

TEST(ZmiIndexRecord, ParseIdentityNoRg)
{
    std::vector<std::byte> aux;
    AppendXiInt(aux, 'q', 's', 42);
    const std::vector<std::byte> body{BuildRecordBody("movie/77/0_50", aux)};
    const auto [rgId, zmw]{detail::ParseRecordIdentity(body)};
    EXPECT_EQ(zmw, 77);
    EXPECT_EQ(rgId, 0);
}

TEST(ZmiIndexRecord, ParseIdentityNameNoSlash)
{
    std::vector<std::byte> aux;
    AppendRgZ(aux, "01");
    const std::vector<std::byte> body{BuildRecordBody("noslash", aux)};
    const auto [rgId, zmw]{detail::ParseRecordIdentity(body)};
    EXPECT_EQ(zmw, 0);
    EXPECT_EQ(rgId, 1);
}

// --- RecordParser ---

TEST(ZmiIndexParser, EmptyHeaderNoRecords)
{
    RecordParser parser;
    const std::vector<std::byte> header{BuildEmptyBamHeader()};
    parser.Feed(/*fileOffset=*/0x1000U, header);
    EXPECT_FALSE(parser.NextRecord().has_value());
    parser.Finish();
    EXPECT_EQ(parser.RecordsEmitted(), 0U);
}

TEST(ZmiIndexParser, SingleRecordSingleBlock)
{
    std::vector<std::byte> aux;
    AppendRgZ(aux, "00000005");
    const std::vector<std::byte> body{BuildRecordBody("m/42/ccs", aux)};
    const std::vector<std::byte> framed{WithBlockSize(body)};

    std::vector<std::byte> blockBytes{BuildEmptyBamHeader()};
    blockBytes.insert(blockBytes.end(), framed.begin(), framed.end());

    RecordParser parser;
    parser.Feed(/*fileOffset=*/0x10000U, blockBytes);
    const auto entry{parser.NextRecord()};
    ASSERT_TRUE(entry.has_value());
    EXPECT_EQ(entry->rgId, 5);
    EXPECT_EQ(entry->zmw, 42);
    EXPECT_EQ(entry->virtualOffset.BlockOffset(), 0x10000U);
    EXPECT_EQ(entry->virtualOffset.WithinBlockOffset(),
              static_cast<std::uint16_t>(BuildEmptyBamHeader().size()));

    EXPECT_FALSE(parser.NextRecord().has_value());
    parser.Finish();
    EXPECT_EQ(parser.RecordsEmitted(), 1U);
}

TEST(ZmiIndexParser, RecordSpansTwoBlocks)
{
    // Header lives entirely in block 0, record straddles into block 1.
    std::vector<std::byte> aux;
    AppendRgZ(aux, "0a");
    const std::vector<std::byte> body{BuildRecordBody("m/9/ccs", aux)};
    const std::vector<std::byte> framed{WithBlockSize(body)};

    // Build block 0: header + first few bytes of the framed record.
    std::vector<std::byte> block0{BuildEmptyBamHeader()};
    const std::size_t splitAt{6U};
    ASSERT_LT(splitAt, framed.size());
    block0.insert(block0.end(), framed.begin(), framed.begin() + splitAt);

    // Build block 1: remaining bytes of the record.
    std::vector<std::byte> block1{framed.begin() + splitAt, framed.end()};

    constexpr std::uint64_t fileOffset0{0x2000U};
    constexpr std::uint64_t fileOffset1{0x4000U};
    const std::size_t recordStartInBlock0{BuildEmptyBamHeader().size()};

    RecordParser parser;
    parser.Feed(fileOffset0, block0);
    EXPECT_FALSE(parser.NextRecord().has_value());
    parser.Feed(fileOffset1, block1);
    const auto entry{parser.NextRecord()};
    ASSERT_TRUE(entry.has_value());
    EXPECT_EQ(entry->zmw, 9);
    EXPECT_EQ(entry->rgId, 0x0A);
    EXPECT_EQ(entry->virtualOffset.BlockOffset(), fileOffset0);
    EXPECT_EQ(entry->virtualOffset.WithinBlockOffset(),
              static_cast<std::uint16_t>(recordStartInBlock0));

    EXPECT_FALSE(parser.NextRecord().has_value());
    parser.Finish();
    EXPECT_EQ(parser.RecordsEmitted(), 1U);
}

TEST(ZmiIndexParser, RecordStartsAtBlockBoundary)
{
    // Block 0 = header only. Record entirely in block 1.
    std::vector<std::byte> aux;
    AppendRgZ(aux, "ff");
    const std::vector<std::byte> body{BuildRecordBody("m/3/ccs", aux)};
    const std::vector<std::byte> framed{WithBlockSize(body)};

    const std::vector<std::byte> block0{BuildEmptyBamHeader()};
    const std::vector<std::byte> block1{framed};

    constexpr std::uint64_t fileOffset0{0x100U};
    constexpr std::uint64_t fileOffset1{0x200U};

    RecordParser parser;
    parser.Feed(fileOffset0, block0);
    EXPECT_FALSE(parser.NextRecord().has_value());
    parser.Feed(fileOffset1, block1);
    const auto entry{parser.NextRecord()};
    ASSERT_TRUE(entry.has_value());
    EXPECT_EQ(entry->virtualOffset.BlockOffset(), fileOffset1);
    EXPECT_EQ(entry->virtualOffset.WithinBlockOffset(), 0U);
    parser.Finish();
}

TEST(ZmiIndexParser, MultipleRecordsSameBlock)
{
    std::vector<std::byte> aux1;
    AppendRgZ(aux1, "01");
    std::vector<std::byte> aux2;
    AppendRgZ(aux2, "02");
    const std::vector<std::byte> framed1{WithBlockSize(BuildRecordBody("m/100/r1", aux1))};
    const std::vector<std::byte> framed2{WithBlockSize(BuildRecordBody("m/200/r2", aux2))};

    std::vector<std::byte> block{BuildEmptyBamHeader()};
    block.insert(block.end(), framed1.begin(), framed1.end());
    block.insert(block.end(), framed2.begin(), framed2.end());

    constexpr std::uint64_t fileOffset{0x5000U};
    const std::size_t headerSize{BuildEmptyBamHeader().size()};

    RecordParser parser;
    parser.Feed(fileOffset, block);

    const auto e1{parser.NextRecord()};
    ASSERT_TRUE(e1.has_value());
    EXPECT_EQ(e1->zmw, 100);
    EXPECT_EQ(e1->rgId, 1);
    EXPECT_EQ(e1->virtualOffset.BlockOffset(), fileOffset);
    EXPECT_EQ(e1->virtualOffset.WithinBlockOffset(), static_cast<std::uint16_t>(headerSize));

    const auto e2{parser.NextRecord()};
    ASSERT_TRUE(e2.has_value());
    EXPECT_EQ(e2->zmw, 200);
    EXPECT_EQ(e2->rgId, 2);
    EXPECT_EQ(e2->virtualOffset.BlockOffset(), fileOffset);
    EXPECT_EQ(e2->virtualOffset.WithinBlockOffset(),
              static_cast<std::uint16_t>(headerSize + framed1.size()));

    EXPECT_FALSE(parser.NextRecord().has_value());
    parser.Finish();
    EXPECT_EQ(parser.RecordsEmitted(), 2U);
}

TEST(ZmiIndexParser, RejectsBadMagic)
{
    RecordParser parser;
    std::vector<std::byte> badHeader{std::byte{'F'}, std::byte{'O'}, std::byte{'O'}, std::byte{0}};
    parser.Feed(0U, badHeader);
    EXPECT_THROW((void)parser.NextRecord(), std::runtime_error);
}

TEST(ZmiIndexParser, FinishThrowsOnPartialRecord)
{
    std::vector<std::byte> aux;
    AppendRgZ(aux, "01");
    const std::vector<std::byte> framed{WithBlockSize(BuildRecordBody("m/1/ccs", aux))};

    std::vector<std::byte> partial{BuildEmptyBamHeader()};
    partial.insert(partial.end(), framed.begin(), framed.begin() + 5);  // truncated

    RecordParser parser;
    parser.Feed(0U, partial);
    EXPECT_FALSE(parser.NextRecord().has_value());
    EXPECT_THROW(parser.Finish(), std::runtime_error);
}

}  // namespace ZmiIndex
}  // namespace Samoa
}  // namespace PacBio
