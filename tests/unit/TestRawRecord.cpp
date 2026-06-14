#include <pbsamoa/core/RawRecord.hpp>

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <gtest/gtest.h>

#include <array>
#include <bit>
#include <span>
#include <stdexcept>
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

// Build a BAM record carrying the kSmN long-CIGAR placeholder (n_cigar_op=2, cigar[0] =
// soft-clip of the whole read) per SAMv1 §4.2.2. \p withCgTag adds the real CIGAR in a
// CG:B,I tag; \p validPlaceholder=false makes cigar[0] a non-clip op so the placeholder
// guard must not fire. Seq is the 4-base "ACGT"; the CG CIGAR (1M1I2M) has query length 4.
std::vector<std::byte> MakeLongCigarPlaceholderBytes(bool withCgTag, bool validPlaceholder = true)
{
    const std::uint32_t seqLen{4};
    std::vector<std::byte> data(32, std::byte{0});
    WriteI32LE(std::data(data) + 0, 0);   // refId (>= 0)
    WriteI32LE(std::data(data) + 4, 0);   // pos (>= 0)
    data[8] = std::byte{3};               // l_read_name "cg\0"
    WriteU16LE(std::data(data) + 12, 2);  // n_cigar_op = 2 (placeholder)
    WriteU32LE(std::data(data) + 16, seqLen);
    WriteI32LE(std::data(data) + 20, -1);  // next_refID
    WriteI32LE(std::data(data) + 24, -1);  // next_pos

    const auto appendU32{[&data](std::uint32_t value) {
        std::array<std::byte, 4> buf{};
        WriteU32LE(std::data(buf), value);
        data.insert(std::end(data), std::begin(buf), std::end(buf));
    }};

    data.push_back(std::byte{'c'});
    data.push_back(std::byte{'g'});
    data.push_back(std::byte{0});

    appendU32(validPlaceholder ? ((seqLen << 4U) | 4U) : ((seqLen << 4U) | 0U));  // 4S or 4M
    appendU32((4U << 4U) | 3U);                                                   // 4N filler

    data.push_back(std::byte{(1 << 4) | 2});  // packed seq "AC"
    data.push_back(std::byte{(4 << 4) | 8});  // packed seq "GT"
    for (int i{0}; i < 4; ++i) {
        data.push_back(std::byte{0xFF});  // qual unavailable
    }

    if (withCgTag) {
        data.push_back(std::byte{'C'});
        data.push_back(std::byte{'G'});
        data.push_back(std::byte{'B'});
        data.push_back(std::byte{'I'});
        appendU32(3U);               // element count
        appendU32((1U << 4U) | 0U);  // 1M
        appendU32((1U << 4U) | 1U);  // 1I
        appendU32((2U << 4U) | 0U);  // 2M
    }
    return data;
}

// Like MakeLongCigarPlaceholderBytes(withCgTag=true) but lets the caller choose
// the B-array subtype byte so we can exercise 'i' (signed int32) alongside 'I'.
// The payload uint32 CIGAR ops are bit-identical for both subtypes.
std::vector<std::byte> MakeLongCigarPlaceholderBytesWithSubtype(char subtype)
{
    std::vector<std::byte> data{MakeLongCigarPlaceholderBytes(/*withCgTag=*/true)};
    // The subtype byte is the 4th byte of the CG tag: C G B <subtype>.
    for (std::size_t i{0}; (i + 3) < std::size(data); ++i) {
        if ((data[i] == std::byte{'C'}) && (data[i + 1] == std::byte{'G'}) &&
            (data[i + 2] == std::byte{'B'})) {
            data[i + 3] = static_cast<std::byte>(subtype);
            return data;
        }
    }
    throw std::logic_error{"CG tag not found in fixture"};
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

TEST(RawRecord, LongCigarPlaceholderExpandsFromCgTag)
{
    const std::vector<std::byte> data{MakeLongCigarPlaceholderBytes(/*withCgTag=*/true)};
    const RawRecord view{std::span<const std::byte>{data}};

    // The 2-op kSmN placeholder must be replaced by the 3-op CIGAR from the CG tag.
    const CigarView cigar{view.CigarOps()};
    ASSERT_EQ(std::size(cigar), 3U);
    EXPECT_EQ(cigar[0].RawValue(), (1U << 4U) | 0U);  // 1M
    EXPECT_EQ(cigar[1].RawValue(), (1U << 4U) | 1U);  // 1I
    EXPECT_EQ(cigar[2].RawValue(), (2U << 4U) | 0U);  // 2M

    // The now-redundant CG tag must be dropped from the decoded tags / owned record.
    EXPECT_FALSE(view.ParseTags().Contains(TagKey{'C', 'G'}));
    const BamRecord owned{view.ToOwned()};
    EXPECT_EQ(std::size(owned.Cigar()), 3U);
    EXPECT_FALSE(owned.Tags().Contains(TagKey{'C', 'G'}));
}

TEST(RawRecord, LongCigarPlaceholderWithoutCgTagIsNotExpanded)
{
    const std::vector<std::byte> data{MakeLongCigarPlaceholderBytes(/*withCgTag=*/false)};
    const RawRecord view{std::span<const std::byte>{data}};

    // No CG tag present — the placeholder must be left as the literal 2-op CIGAR.
    EXPECT_EQ(std::size(view.CigarOps()), 2U);
}

TEST(RawRecord, NonPlaceholderCigarWithCgTagIsNotExpanded)
{
    // cigar[0] is not a full-read soft clip, so the CG tag must be left untouched.
    const std::vector<std::byte> data{
        MakeLongCigarPlaceholderBytes(/*withCgTag=*/true, /*validPlaceholder=*/false)};
    const RawRecord view{std::span<const std::byte>{data}};

    EXPECT_EQ(std::size(view.CigarOps()), 2U);
    EXPECT_TRUE(view.ParseTags().Contains(TagKey{'C', 'G'}));
}

TEST(RawRecord, LongCigarCgTagSubtypeLowercaseIExpandsIdentically)
{
    // htslib bam_tag2cigar (sam.c:703) accepts CG:B,I *and* CG:B,i because the
    // payload uint32 encoding is the same for both subtypes. pbsamoa must match.
    const RawRecord viewUpper{
        std::span<const std::byte>{MakeLongCigarPlaceholderBytesWithSubtype('I')}};
    const RawRecord viewLower{
        std::span<const std::byte>{MakeLongCigarPlaceholderBytesWithSubtype('i')}};

    const CigarView cigarUpper{viewUpper.CigarOps()};
    const CigarView cigarLower{viewLower.CigarOps()};

    // Both must expand to the same 3-op CIGAR.
    ASSERT_EQ(std::size(cigarUpper), 3U);
    ASSERT_EQ(std::size(cigarLower), 3U);
    for (std::size_t idx{0}; idx < 3U; ++idx) {
        EXPECT_EQ(cigarUpper[idx].RawValue(), cigarLower[idx].RawValue()) << "op " << idx;
    }
}

TEST(RawRecord, MappedCigarQueryLengthMustMatchSeq)
{
    // Mapped record, l_seq=1, CIGAR "2M" (query length 2): htslib rejects the mismatch.
    std::vector<std::byte> data(32, std::byte{0});
    data[8] = std::byte{1};               // l_read_name = 1 (just the NUL)
    WriteU16LE(std::data(data) + 12, 1);  // n_cigar_op = 1
    WriteU16LE(std::data(data) + 14, 0);  // flag = 0 (mapped)
    WriteU32LE(std::data(data) + 16, 1);  // l_seq = 1
    data.push_back(std::byte{0});         // read name NUL
    std::array<std::byte, 4> cigar{};
    WriteU32LE(std::data(cigar), (2U << 4U) | 0U);  // 2M
    data.insert(std::end(data), std::begin(cigar), std::end(cigar));
    data.push_back(std::byte{0x10});  // packed seq, 1 base
    data.push_back(std::byte{0xFF});  // qual

    EXPECT_THROW(RawRecord(std::span<const std::byte>{data}), std::invalid_argument);
}

TEST(RawRecord, NonNullTerminatedNameThrows)
{
    std::vector<std::byte> data{MakeLittleEndianRawRecordBytes()};
    data[35] = std::byte{'X'};  // overwrite read-name NUL terminator
    EXPECT_THROW(RawRecord(std::span<const std::byte>{data}), std::invalid_argument);
}

// IsSecondary (0x100) and IsSupplementary (0x800) must be distinguishable —
// IsPrimary alone collapses both, but primary-relabel logic needs them apart.
TEST(RawRecord, IsSecondaryAndIsSupplementaryReflectFlag)
{
    BamRecord secondary;
    secondary.Name("r").Flag(0x100).Sequence("A").Qualities({30});
    const std::vector<std::byte> secData{secondary.SerializeToBam()};
    const RawRecord secView{std::span<const std::byte>{secData}};
    EXPECT_TRUE(secView.IsSecondary());
    EXPECT_FALSE(secView.IsSupplementary());
    EXPECT_FALSE(secView.IsPrimary());

    BamRecord supplementary;
    supplementary.Name("r").Flag(0x800).Sequence("A").Qualities({30});
    const std::vector<std::byte> supData{supplementary.SerializeToBam()};
    const RawRecord supView{std::span<const std::byte>{supData}};
    EXPECT_TRUE(supView.IsSupplementary());
    EXPECT_FALSE(supView.IsSecondary());
    EXPECT_FALSE(supView.IsPrimary());
}

}  // namespace Samoa
}  // namespace PacBio
