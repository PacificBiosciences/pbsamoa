#include <pbsamoa/core/RawRecord.hpp>

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <gtest/gtest.h>

#include <bit>
#include <span>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

namespace {

void WriteU16LE(std::byte* p, std::uint16_t value)
{
    p[0] = static_cast<std::byte>(value & 0xFFU);
    p[1] = static_cast<std::byte>((value >> 8U) & 0xFFU);
}

void WriteU32LE(std::byte* p, std::uint32_t value)
{
    p[0] = static_cast<std::byte>(value & 0xFFU);
    p[1] = static_cast<std::byte>((value >> 8U) & 0xFFU);
    p[2] = static_cast<std::byte>((value >> 16U) & 0xFFU);
    p[3] = static_cast<std::byte>((value >> 24U) & 0xFFU);
}

void WriteI32LE(std::byte* p, std::int32_t value)
{
    WriteU32LE(p, std::bit_cast<std::uint32_t>(value));
}

std::vector<std::byte> MakeTestRecordBytes()
{
    BamRecord rec;
    rec.Name("r001")
        .Flag(99)
        .RefId(0)
        .Pos(6)
        .MapQ(30)
        .Cigar(*ParseCigar("8M2I4M1D3M"))
        .NextRefId(0)
        .NextPos(36)
        .Tlen(39)
        .Sequence("TTAGATAAAGGATACTG")
        .Qualities({30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30});

    TagMap tags;
    tags.Set(TagKey{'N', 'M'}, TagValue{std::int64_t{1}});
    rec.Tags(std::move(tags));

    return rec.SerializeToBam();
}

std::vector<std::byte> MakeLittleEndianRawRecordBytes()
{
    // 32 fixed + name(4) + cigar(8) + packedSeq(1) + qual(1)
    std::vector<std::byte> data(46, std::byte{0});
    const std::int32_t refId{std::bit_cast<std::int32_t>(0x01020304U)};
    const std::int32_t pos{std::bit_cast<std::int32_t>(0x11223344U)};
    const std::int32_t nextRefId{std::bit_cast<std::int32_t>(0xA1B2C3D4U)};
    const std::int32_t nextPos{std::bit_cast<std::int32_t>(0x55667788U)};
    const std::int32_t tlen{std::bit_cast<std::int32_t>(0x90ABCDEFU)};

    WriteI32LE(std::data(data) + 0, refId);
    WriteI32LE(std::data(data) + 4, pos);
    data[8] = std::byte{4};                    // l_read_name ("rLE\0")
    data[9] = static_cast<std::byte>(0x7E);    // mapq
    WriteU16LE(std::data(data) + 10, 0x1234);  // bin
    WriteU16LE(std::data(data) + 12, 2);       // n_cigar_op
    WriteU16LE(std::data(data) + 14, 0xABCD);  // flag
    WriteU32LE(std::data(data) + 16, 1);       // l_seq
    WriteI32LE(std::data(data) + 20, nextRefId);
    WriteI32LE(std::data(data) + 24, nextPos);
    WriteI32LE(std::data(data) + 28, tlen);

    // name "rLE\0"
    data[32] = std::byte{'r'};
    data[33] = std::byte{'L'};
    data[34] = std::byte{'E'};
    data[35] = std::byte{0};

    // CIGAR raw words (little-endian uint32): 0x01020320 (len=0x010203, op=M),
    // 0x04050671 (len=0x0405067, op=I)
    WriteU32LE(std::data(data) + 36, 0x01020320U);
    WriteU32LE(std::data(data) + 40, 0x04050671U);

    // packed seq + qual
    data[44] = std::byte{0x10};             // "A"
    data[45] = static_cast<std::byte>(40);  // quality
    return data;
}

}  // namespace

TEST(RawRecord, ConstructFromSpan)
{
    const std::vector<std::byte> data{MakeTestRecordBytes()};
    const RawRecord view{std::span<const std::byte>{data}};
    EXPECT_FALSE(view.RawData().empty());
}

TEST(RawRecord, OwnsData)
{
    std::vector<std::byte> data{MakeTestRecordBytes()};
    const RawRecord view{std::span<const std::byte>{data}};

    // Overwrite original data — view must still work because it owns a copy
    std::ranges::fill(data, std::byte{0xFF});

    EXPECT_EQ(view.RefId(), 0);
    EXPECT_EQ(view.Pos(), 6);
    EXPECT_EQ(view.Name(), "r001");
}

TEST(RawRecord, FieldAccessors)
{
    const std::vector<std::byte> data{MakeTestRecordBytes()};
    const RawRecord view{std::span<const std::byte>{data}};

    EXPECT_EQ(view.RefId(), 0);
    EXPECT_EQ(view.Pos(), 6);
    EXPECT_EQ(view.MapQ(), 30u);
    EXPECT_EQ(view.Flag(), 99u);
    EXPECT_EQ(view.NextRefId(), 0);
    EXPECT_EQ(view.NextPos(), 36);
    EXPECT_EQ(view.Tlen(), 39);
    EXPECT_EQ(view.Name(), "r001");
    EXPECT_EQ(view.SeqLength(), 17);
    EXPECT_EQ(view.CigarOpCount(), 5u);
}

TEST(RawRecord, ToOwnedProducesEquivalentRecord)
{
    const std::vector<std::byte> data{MakeTestRecordBytes()};
    const RawRecord view{std::span<const std::byte>{data}};

    const BamRecord owned{view.ToOwned()};

    EXPECT_EQ(owned.Name(), "r001");
    EXPECT_EQ(owned.Flag(), 99u);
    EXPECT_EQ(owned.Pos(), 6);
    EXPECT_EQ(owned.Tlen(), 39);
    EXPECT_EQ(owned.Sequence(), "TTAGATAAAGGATACTG");
}

TEST(RawRecord, MoveSemantics)
{
    const std::vector<std::byte> data{MakeTestRecordBytes()};
    RawRecord view1{std::span<const std::byte>{data}};

    const RawRecord view2{std::move(view1)};
    EXPECT_EQ(view2.Pos(), 6);
    EXPECT_EQ(view2.Name(), "r001");
}

TEST(RawRecord, TooSmallForFixedFieldsThrows)
{
    const std::vector<std::byte> data(16, std::byte{0});
    EXPECT_THROW(RawRecord(std::span<const std::byte>{data}), std::invalid_argument);
}

TEST(RawRecord, EmptyDataThrows)
{
    const std::vector<std::byte> data;
    EXPECT_THROW(RawRecord(std::span<const std::byte>{data}), std::invalid_argument);
}

TEST(RawRecord, ZeroNameLengthThrows)
{
    // Create 32 bytes (minimum fixed fields) with l_read_name = 0 at offset 8
    std::vector<std::byte> data(32, std::byte{0});
    EXPECT_THROW(RawRecord(std::span<const std::byte>{data}), std::invalid_argument);
}

TEST(RawRecord, TruncatedVariableLengthFieldsThrows)
{
    // Create valid fixed fields but claim name length > remaining buffer
    std::vector<std::byte> data(33, std::byte{0});
    data[8] = std::byte{10};  // l_read_name = 10, but only 1 byte after fixed fields
    EXPECT_THROW(RawRecord(std::span<const std::byte>{data}), std::invalid_argument);
}

TEST(RawRecord, ExactMinimumSizeIsValid)
{
    // 32 fixed bytes + 1 byte name (l_read_name=1) + no CIGAR + no seq + no qual
    std::vector<std::byte> data(33, std::byte{0});
    data[8] = std::byte{1};  // l_read_name = 1 (just the NUL)
    EXPECT_NO_THROW(RawRecord(std::span<const std::byte>{data}));
}

TEST(RawRecord, FixedFieldsAndCigarDecodeFromLittleEndianBytes)
{
    const std::vector<std::byte> data{MakeLittleEndianRawRecordBytes()};
    const RawRecord view{std::span<const std::byte>{data}};

    EXPECT_EQ(view.RefId(), std::bit_cast<std::int32_t>(0x01020304U));
    EXPECT_EQ(view.Pos(), std::bit_cast<std::int32_t>(0x11223344U));
    EXPECT_EQ(view.MapQ(), 0x7EU);
    EXPECT_EQ(view.Bin(), 0x1234U);
    EXPECT_EQ(view.CigarOpCount(), 2U);
    EXPECT_EQ(view.Flag(), 0xABCDU);
    EXPECT_EQ(view.SeqLength(), 1U);
    EXPECT_EQ(view.NextRefId(), std::bit_cast<std::int32_t>(0xA1B2C3D4U));
    EXPECT_EQ(view.NextPos(), std::bit_cast<std::int32_t>(0x55667788U));
    EXPECT_EQ(view.Tlen(), std::bit_cast<std::int32_t>(0x90ABCDEFU));
    EXPECT_EQ(view.Name(), "rLE");

    const CigarView cigar{view.CigarOps()};
    ASSERT_EQ(std::size(cigar), 2U);
    EXPECT_EQ(cigar[0].RawValue(), 0x01020320U);
    EXPECT_EQ(cigar[1].RawValue(), 0x04050671U);
}

TEST(RawRecord, NonNullTerminatedNameThrows)
{
    std::vector<std::byte> data{MakeLittleEndianRawRecordBytes()};
    data[35] = std::byte{'X'};  // overwrite read-name NUL terminator
    EXPECT_THROW(RawRecord(std::span<const std::byte>{data}), std::invalid_argument);
}

}  // namespace Samoa
}  // namespace PacBio
