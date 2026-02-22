#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/CigarClipping.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/TagClipping.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <span>
#include <ranges>

#include <gtest/gtest.h>

#include <string>
#include <vector>

#include <cstdint>
#include <cstring>

namespace PacBio {
namespace Samoa {

namespace {
BamRecord MakeSimpleRecord()
{
    BamRecord rec;
    rec.Name("read1")
        .Flag(0)  // forward strand, mapped
        .RefId(0)
        .Pos(100)
        .MapQ(60)
        .Cigar(ParseCigar("10M"))
        .Sequence("ACGTACGTAC")
        .Qualities(std::vector<std::uint8_t>{30, 31, 32, 33, 34, 35, 36, 37, 38, 39});

    TagMap tags;
    tags.Set(TagKey{'q', 's'}, std::int64_t{100});
    tags.Set(TagKey{'q', 'e'}, std::int64_t{110});
    rec.Tags(std::move(tags));
    return rec;
}

BamRecord MakeRecordWithTags()
{
    auto rec{MakeSimpleRecord()};
    TagMap tags{rec.Tags()};

    // Add kinetics tag (per-base, forward orientation for fwd strand)
    TagArray ipArr{'C'};
    for (std::uint8_t i{0}; i < 10; ++i) {
        ipArr.AppendUInt8(i * 10);
    }
    tags.Set(TagKey{'i', 'p'}, std::move(ipArr));

    rec.Tags(std::move(tags));
    return rec;
}
}  // namespace

TEST(BamRecordClipping, ClipToQuery_Simple_Forward)
{
    auto rec{MakeSimpleRecord()};
    rec.Clip(ClipType::CLIP_TO_QUERY, 102, 108);

    EXPECT_EQ(rec.Sequence(), "GTACGT");
    ASSERT_EQ(std::size(rec.Qualities()), 6U);
    EXPECT_EQ(rec.Qualities()[0], 32);
    EXPECT_EQ(rec.Qualities()[5], 37);
    EXPECT_EQ(CigarToString(rec.Cigar()), "6M");
    EXPECT_EQ(rec.Pos(), 102);

    const auto* qs{rec.Tags().Get(TagKey{'q', 's'})};
    ASSERT_NE(qs, nullptr);
    EXPECT_EQ(std::get<std::int64_t>(*qs), 102);

    const auto* qe{rec.Tags().Get(TagKey{'q', 'e'})};
    ASSERT_NE(qe, nullptr);
    EXPECT_EQ(std::get<std::int64_t>(*qe), 108);
}

TEST(BamRecordClipping, ClipToReference_Simple_Forward)
{
    auto rec{MakeSimpleRecord()};
    rec.Clip(ClipType::CLIP_TO_REFERENCE, 102, 108);

    EXPECT_EQ(rec.Sequence(), "GTACGT");
    EXPECT_EQ(CigarToString(rec.Cigar()), "6M");
    EXPECT_EQ(rec.Pos(), 102);
}

TEST(BamRecordClipping, Clipped_ReturnsCopy)
{
    const auto original{MakeSimpleRecord()};
    const auto clipped{original.Clipped(ClipType::CLIP_TO_QUERY, 102, 108)};

    // Original unchanged
    EXPECT_EQ(original.Sequence(), "ACGTACGTAC");
    EXPECT_EQ(std::size(original.Qualities()), 10U);

    // Clipped is correct
    EXPECT_EQ(clipped.Sequence(), "GTACGT");
    EXPECT_EQ(std::size(clipped.Qualities()), 6U);
}

TEST(BamRecordClipping, ClipToReference_Unmapped_NoOp)
{
    BamRecord rec;
    rec.Name("unmapped")
        .Flag(0x4)  // unmapped
        .Sequence("ACGT")
        .Qualities(std::vector<std::uint8_t>{30, 31, 32, 33});

    rec.Clip(ClipType::CLIP_TO_REFERENCE, 0, 2);

    EXPECT_EQ(rec.Sequence(), "ACGT");
    EXPECT_EQ(std::size(rec.Qualities()), 4U);
}

TEST(BamRecordClipping, ClipToQuery_WithKineticsTags)
{
    auto rec{MakeRecordWithTags()};
    rec.Clip(ClipType::CLIP_TO_QUERY, 102, 108);

    const auto* ip{rec.Tags().Get(TagKey{'i', 'p'})};
    ASSERT_NE(ip, nullptr);
    const auto& arr{std::get<TagArray>(*ip)};
    EXPECT_EQ(arr.Count(), 6U);
    EXPECT_EQ(static_cast<std::uint8_t>(arr.Data()[0]), 20);  // index 2
    EXPECT_EQ(static_cast<std::uint8_t>(arr.Data()[5]), 70);  // index 7
}

TEST(BamRecordClipping, ClipToQuery_NoClipNeeded)
{
    auto rec{MakeSimpleRecord()};
    rec.Clip(ClipType::CLIP_TO_QUERY, 100, 110);

    EXPECT_EQ(rec.Sequence(), "ACGTACGTAC");
    EXPECT_EQ(CigarToString(rec.Cigar()), "10M");
    EXPECT_EQ(rec.Pos(), 100);
}

TEST(BamRecordClipping, ClipToQuery_WithMixedCigar)
{
    BamRecord rec;
    rec.Name("read2")
        .Flag(0)
        .RefId(0)
        .Pos(100)
        .MapQ(60)
        .Cigar(ParseCigar("3M2I5M"))  // query length = 10
        .Sequence("ACGTACGTAC")
        .Qualities(std::vector<std::uint8_t>{30, 31, 32, 33, 34, 35, 36, 37, 38, 39});

    TagMap tags;
    tags.Set(TagKey{'q', 's'}, std::int64_t{0});
    tags.Set(TagKey{'q', 'e'}, std::int64_t{10});
    rec.Tags(std::move(tags));

    rec.Clip(ClipType::CLIP_TO_QUERY, 4, 10);

    EXPECT_EQ(std::size(rec.Sequence()), 6U);
    EXPECT_EQ(CigarToString(rec.Cigar()), "1I5M");
    EXPECT_EQ(rec.Pos(), 103);
}

TEST(BamRecordClipping, ClipToReference_WithDeletion)
{
    BamRecord rec;
    rec.Name("read3")
        .Flag(0)
        .RefId(0)
        .Pos(50)
        .MapQ(60)
        .Cigar(ParseCigar("3M2D5M"))  // ref span = 10, query length = 8
        .Sequence("ABCDEFGH")
        .Qualities(std::vector<std::uint8_t>{10, 20, 30, 40, 50, 60, 70, 80});

    TagMap tags;
    tags.Set(TagKey{'q', 's'}, std::int64_t{0});
    tags.Set(TagKey{'q', 'e'}, std::int64_t{8});
    rec.Tags(std::move(tags));

    rec.Clip(ClipType::CLIP_TO_REFERENCE, 54, 58);

    EXPECT_EQ(CigarToString(rec.Cigar()), "1D3M");
    EXPECT_EQ(rec.Pos(), 54);
    EXPECT_EQ(std::size(rec.Sequence()), 3U);
    EXPECT_EQ(rec.Sequence(), "DEF");
}

TEST(BamRecordClipping, ClipToQuery_CCSRecord_NoQsQe)
{
    // CCS records have no qs/qe tags. QueryStart defaults to 0, QueryEnd to seqLen.
    BamRecord rec;
    rec.Name("movie/42/ccs")
        .Flag(0)
        .RefId(0)
        .Pos(100)
        .MapQ(60)
        .Cigar(ParseCigar("10M"))
        .Sequence("ACGTACGTAC")
        .Qualities(std::vector<std::uint8_t>{30, 31, 32, 33, 34, 35, 36, 37, 38, 39});

    rec.Clip(ClipType::CLIP_TO_QUERY, 2, 8);

    EXPECT_EQ(rec.Sequence(), "GTACGT");
    EXPECT_EQ(CigarToString(rec.Cigar()), "6M");
    EXPECT_EQ(rec.Pos(), 102);
}

TEST(BamRecordClipping, MutableTags_ReturnsReference)
{
    auto rec{MakeSimpleRecord()};
    TagMap& tags{rec.MutableTags()};
    tags.Set(TagKey{'X', 'Y'}, std::int64_t{42});

    const auto* val{rec.Tags().Get(TagKey{'X', 'Y'})};
    ASSERT_NE(val, nullptr);
    EXPECT_EQ(std::get<std::int64_t>(*val), 42);
}

// =====================================================================
// Helper: extract uint8 values from a TagArray for comparison
// =====================================================================
namespace {
std::vector<std::uint8_t> ExtractUInt8(const TagArray& arr)
{
    std::vector<std::uint8_t> result;
    result.reserve(arr.Count());
    for (std::uint32_t i{0}; i < arr.Count(); ++i) {
        result.push_back(static_cast<std::uint8_t>(arr.Data()[i]));
    }
    return result;
}

std::vector<std::uint16_t> ExtractUInt16(const TagArray& arr)
{
    std::vector<std::uint16_t> result;
    result.reserve(arr.Count());
    const auto data{arr.Data()};
    for (std::uint32_t i{0}; i < arr.Count(); ++i) {
        std::uint16_t v{0};
        std::memcpy(&v, &data[i * 2], 2);
        result.push_back(v);
    }
    return result;
}

// Build a mapped record with per-base tags for clipping pipeline tests.
// seq: 10 bases "AACCGTTAGC", quals: {10,10,20,20,30,40,40,10,30,20}
// qs=500, qe=510, ip/pw = same as quals (uint8), dq = string same as seq
BamRecord MakePipelineRecord(std::string_view cigarStr, bool isReverse, std::int32_t pos = 100)
{
    const std::string seq{"AACCGTTAGC"};
    const std::vector<std::uint8_t> quals{10, 10, 20, 20, 30, 40, 40, 10, 30, 20};

    BamRecord rec;
    rec.Name("read/42/0_10")
        .Flag(isReverse ? std::uint16_t{0x10} : std::uint16_t{0})
        .RefId(0)
        .Pos(pos)
        .MapQ(80)
        .Cigar(ParseCigar(cigarStr))
        .Sequence(seq)
        .Qualities(quals);

    TagMap tags;
    tags.Set(TagKey{'q', 's'}, std::int64_t{500});
    tags.Set(TagKey{'q', 'e'}, std::int64_t{510});

    // Per-base uint8 tags (SubstringClipStrategy)
    TagArray ipArr{'C'};
    TagArray pwArr{'C'};
    for (const std::uint8_t v : quals) {
        ipArr.AppendUInt8(v);
        pwArr.AppendUInt8(v);
    }
    tags.Set(TagKey{'i', 'p'}, std::move(ipArr));
    tags.Set(TagKey{'p', 'w'}, std::move(pwArr));

    // String per-base tag
    tags.Set(TagKey{'d', 'q'}, std::string{seq});

    rec.Tags(std::move(tags));
    return rec;
}
}  // namespace

// =====================================================================
// 1. Forward-strand query clipping with multiple CIGAR patterns
// =====================================================================

TEST(BamRecordClipping, ClipToQuery_Fwd_10Eq)
{
    // 10= forward, clip query [502,509) -> remove 2 front, 1 back
    auto rec{MakePipelineRecord("10=", false)};
    rec.Clip(ClipType::CLIP_TO_QUERY, 502, 509);

    EXPECT_EQ(rec.Sequence(), "CCGTTAG");
    EXPECT_EQ(CigarToString(rec.Cigar()), "7=");
    EXPECT_EQ(rec.Pos(), 102);  // advanced by 2
    EXPECT_EQ(rec.ReferenceEnd(), 109);

    // ip tag should be clipped [2..9) = {20,20,30,40,40,10,30}
    const auto* ip{rec.Tags().Get(TagKey{'i', 'p'})};
    ASSERT_NE(ip, nullptr);
    const auto ipVals{ExtractUInt8(std::get<TagArray>(*ip))};
    EXPECT_EQ(ipVals, (std::vector<std::uint8_t>{20, 20, 30, 40, 40, 10, 30}));

    // dq string tag
    const auto* dq{rec.Tags().Get(TagKey{'d', 'q'})};
    ASSERT_NE(dq, nullptr);
    EXPECT_EQ(std::get<std::string>(*dq), "CCGTTAG");
}

TEST(BamRecordClipping, ClipToQuery_Fwd_5Eq3D5Eq)
{
    // 5=3D5= forward, clip query [502,509)
    // Query length = 10, ref span = 5+3+5 = 13
    auto rec{MakePipelineRecord("5=3D5=", false)};
    rec.Clip(ClipType::CLIP_TO_QUERY, 502, 509);

    EXPECT_EQ(rec.Sequence(), "CCGTTAG");
    EXPECT_EQ(CigarToString(rec.Cigar()), "3=3D4=");
    EXPECT_EQ(rec.Pos(), 102);
    EXPECT_EQ(rec.ReferenceEnd(), 112);  // 102 + 3 + 3 + 4

    const auto* ip{rec.Tags().Get(TagKey{'i', 'p'})};
    ASSERT_NE(ip, nullptr);
    EXPECT_EQ(std::get<TagArray>(*ip).Count(), 7U);
}

TEST(BamRecordClipping, ClipToQuery_Fwd_4Eq1D2I2D4Eq)
{
    // 4=1D2I2D4= forward, clip query [502,509)
    // Query length = 4+2+4 = 10, ref span = 4+1+2+4 = 11
    auto rec{MakePipelineRecord("4=1D2I2D4=", false)};
    rec.Clip(ClipType::CLIP_TO_QUERY, 502, 509);

    EXPECT_EQ(rec.Sequence(), "CCGTTAG");
    EXPECT_EQ(CigarToString(rec.Cigar()), "2=1D2I2D3=");
    EXPECT_EQ(rec.Pos(), 102);
    // ref span of clipped: 2 + 1 + 2 + 3 = 8
    EXPECT_EQ(rec.ReferenceEnd(), 110);

    const auto* ip{rec.Tags().Get(TagKey{'i', 'p'})};
    ASSERT_NE(ip, nullptr);
    EXPECT_EQ(std::get<TagArray>(*ip).Count(), 7U);
    const auto ipVals{ExtractUInt8(std::get<TagArray>(*ip))};
    EXPECT_EQ(ipVals, (std::vector<std::uint8_t>{20, 20, 30, 40, 40, 10, 30}));
}

// =====================================================================
// 2. Reverse-strand query clipping
// =====================================================================

TEST(BamRecordClipping, ClipToQuery_Rev_10Eq)
{
    // Reverse strand: SEQ in BAM = revcomp of original query.
    // For reverse-strand, query clipping [502,509) removes 2 from original
    // front, 1 from original back. But in BAM SEQ order (reversed), that
    // means remove 1 from BAM front, 2 from BAM back.
    // So clipOffset=1, clipLength=7 into the BAM SEQ.
    auto rec{MakePipelineRecord("10=", true)};
    rec.Clip(ClipType::CLIP_TO_QUERY, 502, 509);

    // SEQ was "AACCGTTAGC" (10 bases), after clip: offset 1, length 7
    EXPECT_EQ(rec.Sequence(), "ACCGTTA");
    EXPECT_EQ(CigarToString(rec.Cigar()), "7=");
    EXPECT_EQ(rec.Pos(), 101);  // 100 + 1 (trimmed 1 from front of BAM)
    EXPECT_EQ(rec.ReferenceEnd(), 108);

    // ip tag (SubstringClipStrategy): clips same as SEQ [1..8)
    const auto* ip{rec.Tags().Get(TagKey{'i', 'p'})};
    ASSERT_NE(ip, nullptr);
    const auto ipVals{ExtractUInt8(std::get<TagArray>(*ip))};
    EXPECT_EQ(ipVals, (std::vector<std::uint8_t>{10, 20, 20, 30, 40, 40, 10}));
}

TEST(BamRecordClipping, ClipToQuery_Rev_5Eq3D5Eq)
{
    // 5=3D5= reverse, clip query [502,509)
    auto rec{MakePipelineRecord("5=3D5=", true)};
    rec.Clip(ClipType::CLIP_TO_QUERY, 502, 509);

    EXPECT_EQ(rec.Sequence(), "ACCGTTA");
    EXPECT_EQ(CigarToString(rec.Cigar()), "4=3D3=");
    EXPECT_EQ(rec.Pos(), 101);
    // ref span: 4 + 3 + 3 = 10, so refEnd = 111
    EXPECT_EQ(rec.ReferenceEnd(), 111);
}

TEST(BamRecordClipping, ClipToQuery_Rev_4Eq1D2I2D4Eq)
{
    // 4=1D2I2D4= reverse, clip query [502,509)
    // queryLen=10, refSpan=11. qs=500, qe=510.
    // frontRemove=2, backRemove=1, swap for reverse -> front=1, back=2
    // Walk front=1 through 4=: partial, 3= left, refAdvance=1, queryOffset=1
    // Walk back=2 through ...4=: shrink to 2=
    // Result: 3=1D2I2D2=, clipOffset=1, clipLength=7, newPos=101
    auto rec{MakePipelineRecord("4=1D2I2D4=", true)};
    rec.Clip(ClipType::CLIP_TO_QUERY, 502, 509);

    EXPECT_EQ(rec.Sequence(), "ACCGTTA");
    EXPECT_EQ(CigarToString(rec.Cigar()), "3=1D2I2D2=");
    EXPECT_EQ(rec.Pos(), 101);
    EXPECT_EQ(rec.ReferenceEnd(), 109);  // 101 + 3 + 1 + 2 + 2 = 109

    // ip: SubstringClip [1..8) = {10,20,20,30,40,40,10}
    const auto* ip{rec.Tags().Get(TagKey{'i', 'p'})};
    ASSERT_NE(ip, nullptr);
    const auto ipVals{ExtractUInt8(std::get<TagArray>(*ip))};
    EXPECT_EQ(ipVals, (std::vector<std::uint8_t>{10, 20, 20, 30, 40, 40, 10}));

    // dq: SubstringClip [1..8) = "ACCGTTA"
    const auto* dq{rec.Tags().Get(TagKey{'d', 'q'})};
    ASSERT_NE(dq, nullptr);
    EXPECT_EQ(std::get<std::string>(*dq), "ACCGTTA");
}

// =====================================================================
// 3. Forward-strand reference clipping
// =====================================================================

TEST(BamRecordClipping, ClipToReference_Fwd_10Eq)
{
    // 10= forward, ref clip [102,107) -> remove 2 ref bases front, 3 back
    auto rec{MakePipelineRecord("10=", false)};
    rec.Clip(ClipType::CLIP_TO_REFERENCE, 102, 107);

    EXPECT_EQ(rec.Sequence(), "CCGTT");
    EXPECT_EQ(CigarToString(rec.Cigar()), "5=");
    EXPECT_EQ(rec.Pos(), 102);
    EXPECT_EQ(rec.ReferenceEnd(), 107);

    // ip clipped to [2..7) = {20,20,30,40,40}
    const auto* ip{rec.Tags().Get(TagKey{'i', 'p'})};
    ASSERT_NE(ip, nullptr);
    const auto ipVals{ExtractUInt8(std::get<TagArray>(*ip))};
    EXPECT_EQ(ipVals, (std::vector<std::uint8_t>{20, 20, 30, 40, 40}));

    // qs/qe updated
    const auto* qs{rec.Tags().Get(TagKey{'q', 's'})};
    ASSERT_NE(qs, nullptr);
    EXPECT_EQ(std::get<std::int64_t>(*qs), 502);
    const auto* qe{rec.Tags().Get(TagKey{'q', 'e'})};
    ASSERT_NE(qe, nullptr);
    EXPECT_EQ(std::get<std::int64_t>(*qe), 507);
}

TEST(BamRecordClipping, ClipToReference_Fwd_5Eq3D5Eq)
{
    // 5=3D5= forward, ref span = 13, ref clip [102,107)
    // Front: remove 2 ref bases from 5=, leaving 3=
    // Back: need to remove 6 ref bases from end. 3D5= = 8 ref bases at end.
    //   Remove from back of 5= first: 5 < 6, fully consumed. Remove 1 from 3D: leaves 2D.
    //   But 3= remains. So clipped CIGAR = 3=2D
    auto rec{MakePipelineRecord("5=3D5=", false)};
    rec.Clip(ClipType::CLIP_TO_REFERENCE, 102, 107);

    EXPECT_EQ(rec.Sequence(), "CCG");
    EXPECT_EQ(CigarToString(rec.Cigar()), "3=2D");
    EXPECT_EQ(rec.Pos(), 102);
    EXPECT_EQ(rec.ReferenceEnd(), 107);  // 102 + 3 + 2
}

TEST(BamRecordClipping, ClipToReference_Fwd_4Eq1D2I2D4Eq)
{
    // 4=1D2I2D4= forward, ref span = 4+1+2+4 = 11, ref clip [102,107)
    // Front: remove 2 ref from 4= -> 2=
    // Back: remove 4 ref from end (refEnd = 111, need to remove 111-107 = 4)
    //   Trailing 4= has 4 ref, fully consumed. Removal stops.
    //   Remaining: 2=1D2I2D (ref span = 2+1+0+2 = 5, matches 107-102)
    auto rec{MakePipelineRecord("4=1D2I2D4=", false)};
    rec.Clip(ClipType::CLIP_TO_REFERENCE, 102, 107);

    EXPECT_EQ(CigarToString(rec.Cigar()), "2=1D2I2D");
    EXPECT_EQ(rec.Pos(), 102);
    EXPECT_EQ(rec.ReferenceEnd(), 107);  // 102 + 2 + 1 + 2 = 107
    EXPECT_EQ(rec.Sequence(), "CCGT");   // 2= + 2I = 4 query bases
}

// =====================================================================
// 4. Reverse-strand reference clipping
// =====================================================================

TEST(BamRecordClipping, ClipToReference_Rev_10Eq)
{
    // 10= reverse, ref clip [102,107)
    // In BAM, CIGAR/POS are in reference coords (same for both strands).
    // Remove 2 ref from front -> 2 BAM-seq bases from front
    // Remove 3 ref from back -> 3 BAM-seq bases from back
    // But for reverse strand, front/back query removal is swapped.
    // So clipOffset = seqLen - queryRemovedFront - queryRemovedBack... actually:
    // ClipCigarToReference swaps front/back for reverse.
    auto rec{MakePipelineRecord("10=", true)};
    rec.Clip(ClipType::CLIP_TO_REFERENCE, 102, 107);

    EXPECT_EQ(CigarToString(rec.Cigar()), "5=");
    EXPECT_EQ(rec.Pos(), 102);
    EXPECT_EQ(rec.ReferenceEnd(), 107);
    // For reverse: queryRemovedFront=2, queryRemovedBack=3 (from ref),
    // then swap -> clipOffset=3 (from BAM-SEQ start), clipLength=5
    EXPECT_EQ(rec.Sequence(), "CGTTA");  // BAM SEQ [3..8)
    EXPECT_EQ(std::size(rec.Qualities()), 5U);
}

TEST(BamRecordClipping, ClipToReference_Rev_5Eq3D5Eq)
{
    // 5=3D5= reverse, ref clip [102,107)
    auto rec{MakePipelineRecord("5=3D5=", true)};
    rec.Clip(ClipType::CLIP_TO_REFERENCE, 102, 107);

    EXPECT_EQ(CigarToString(rec.Cigar()), "3=2D");
    EXPECT_EQ(rec.Pos(), 102);
    // queryRemovedFront(before swap) = 2, queryRemovedBack(before swap) = 5
    // After swap for reverse: clipOffset=5, clipLength=10-2-5=3
    // BAM SEQ = "AACCGTTAGC"[5..8) = "TTA"
    EXPECT_EQ(rec.Sequence(), "TTA");
    EXPECT_EQ(std::size(rec.Qualities()), 3U);
}

TEST(BamRecordClipping, ClipToReference_Rev_4Eq1D2I2D4Eq)
{
    // 4=1D2I2D4= reverse, ref clip [102,107)
    // refSpan=11, refEnd=111. frontRefRemove=2, backRefRemove=4.
    // Walk front=2: 4= partial -> 2=, queryRemovedFront=2
    // Walk back=4: 4= (4 ref, 4 query) fully consumed, queryRemovedBack=4
    // Remaining CIGAR: 2=1D2I2D
    // Reverse swap: clipOffset=4, clipLength=10-2-4=4
    // BAM SEQ "AACCGTTAGC"[4..8) = "GTTA"
    auto rec{MakePipelineRecord("4=1D2I2D4=", true)};
    rec.Clip(ClipType::CLIP_TO_REFERENCE, 102, 107);

    EXPECT_EQ(CigarToString(rec.Cigar()), "2=1D2I2D");
    EXPECT_EQ(rec.Pos(), 102);
    EXPECT_EQ(rec.ReferenceEnd(), 107);  // 102 + 2 + 1 + 2 = 107
    EXPECT_EQ(rec.Sequence(), "GTTA");   // 2= + 2I = 4 query bases
    EXPECT_EQ(std::size(rec.Qualities()), 4U);

    // ip: SubstringClip [4..8) = {30,40,40,10}
    const auto* ip{rec.Tags().Get(TagKey{'i', 'p'})};
    ASSERT_NE(ip, nullptr);
    const auto ipVals{ExtractUInt8(std::get<TagArray>(*ip))};
    EXPECT_EQ(ipVals, (std::vector<std::uint8_t>{30, 40, 40, 10}));
}

// =====================================================================
// 5. Soft-clip handling with reference clipping
// =====================================================================

TEST(BamRecordClipping, ClipToReference_SoftClips_Forward)
{
    // 2S10=3S: 15 query bases, 10 aligned ref bases, pos=100, refEnd=110
    // ref clip [102,107) -> remove 2 ref from front, 3 from back of aligned
    BamRecord rec;
    const std::string seq{"TTAACCGTTAGCAAA"};
    std::vector<std::uint8_t> quals(15);
    for (std::uint8_t i{0}; i < 15; ++i) {
        quals[i] = i + 10;
    }
    rec.Name("softclip_test")
        .Flag(0)
        .RefId(0)
        .Pos(100)
        .MapQ(80)
        .Cigar(ParseCigar("2S10=3S"))
        .Sequence(seq)
        .Qualities(std::move(quals));

    TagMap tags;
    tags.Set(TagKey{'q', 's'}, std::int64_t{500});
    tags.Set(TagKey{'q', 'e'}, std::int64_t{515});

    TagArray ipArr{'C'};
    for (std::uint8_t i{0}; i < 15; ++i) {
        ipArr.AppendUInt8(i * 5);
    }
    tags.Set(TagKey{'i', 'p'}, std::move(ipArr));
    rec.Tags(std::move(tags));

    rec.Clip(ClipType::CLIP_TO_REFERENCE, 102, 107);

    // Ref clip removes soft clips entirely, plus trims aligned portion.
    // Front: 2S (2 query, 0 ref) + 2= (2 ref, 2 query) -> 4 query removed front
    // Back: 3S (3 query, 0 ref) + 3= (3 ref, 3 query) -> 6 query removed back
    // Remaining: 5= (5 query bases), clipOffset=4, clipLength=5
    EXPECT_EQ(CigarToString(rec.Cigar()), "5=");
    EXPECT_EQ(rec.Pos(), 102);
    EXPECT_EQ(rec.ReferenceEnd(), 107);
    EXPECT_EQ(rec.Sequence(), "CCGTT");  // seq[4..9)

    // ip tag clipped to [4..9) = {20, 25, 30, 35, 40}
    const auto* ip{rec.Tags().Get(TagKey{'i', 'p'})};
    ASSERT_NE(ip, nullptr);
    const auto ipVals{ExtractUInt8(std::get<TagArray>(*ip))};
    EXPECT_EQ(ipVals, (std::vector<std::uint8_t>{20, 25, 30, 35, 40}));
}

// =====================================================================
// 6. CCS records (no qs/qe tags) with per-base tags
// =====================================================================

TEST(BamRecordClipping, ClipToQuery_CCS_WithTags)
{
    // CCS: no qs/qe, queryStart=0, queryEnd=seqLen
    BamRecord rec;
    rec.Name("movie/42/ccs")
        .Flag(0)
        .RefId(0)
        .Pos(100)
        .MapQ(60)
        .Cigar(ParseCigar("10="))
        .Sequence("AACCGTTAGC")
        .Qualities(std::vector<std::uint8_t>{10, 10, 20, 20, 30, 40, 40, 10, 30, 20});

    TagMap tags;
    TagArray ipArr{'C'};
    for (std::uint8_t i{0}; i < 10; ++i) {
        ipArr.AppendUInt8(i * 10);
    }
    tags.Set(TagKey{'i', 'p'}, std::move(ipArr));
    rec.Tags(std::move(tags));

    rec.Clip(ClipType::CLIP_TO_QUERY, 2, 7);

    EXPECT_EQ(rec.Sequence(), "CCGTT");
    EXPECT_EQ(CigarToString(rec.Cigar()), "5=");
    EXPECT_EQ(rec.Pos(), 102);

    // ip clipped to [2..7) = {20,30,40,50,60}
    const auto* ip{rec.Tags().Get(TagKey{'i', 'p'})};
    ASSERT_NE(ip, nullptr);
    const auto ipVals{ExtractUInt8(std::get<TagArray>(*ip))};
    EXPECT_EQ(ipVals, (std::vector<std::uint8_t>{20, 30, 40, 50, 60}));

    // Verify no qs/qe was created
    EXPECT_EQ(rec.Tags().Get(TagKey{'q', 's'}), nullptr);
    EXPECT_EQ(rec.Tags().Get(TagKey{'q', 'e'}), nullptr);
}

TEST(BamRecordClipping, ClipToReference_CCS_WithTags)
{
    BamRecord rec;
    rec.Name("movie/42/ccs")
        .Flag(0)
        .RefId(0)
        .Pos(100)
        .MapQ(60)
        .Cigar(ParseCigar("5=3D5="))
        .Sequence("AACCGTTAGC")
        .Qualities(std::vector<std::uint8_t>{10, 10, 20, 20, 30, 40, 40, 10, 30, 20});

    TagMap tags;
    TagArray ipArr{'C'};
    for (std::uint8_t i{0}; i < 10; ++i) {
        ipArr.AppendUInt8(i * 10);
    }
    tags.Set(TagKey{'i', 'p'}, std::move(ipArr));
    rec.Tags(std::move(tags));

    rec.Clip(ClipType::CLIP_TO_REFERENCE, 102, 107);

    // Front: remove 2 ref from 5= -> 3=
    // Back: refEnd=113, need to remove 6. 5= (5 ref) consumed, 1 from 3D -> 2D
    // CIGAR: 3=2D, clipOffset=2, clipLength=3
    EXPECT_EQ(CigarToString(rec.Cigar()), "3=2D");
    EXPECT_EQ(rec.Pos(), 102);
    EXPECT_EQ(rec.Sequence(), "CCG");

    // No qs/qe should be created
    EXPECT_EQ(rec.Tags().Get(TagKey{'q', 's'}), nullptr);
}

// =====================================================================
// 7. Flanking inserts
// =====================================================================

TEST(BamRecordClipping, FlankingInserts_Forward_KeepVsExcise)
{
    // CIGAR: 3=6I10=6I1= (26 query bases, ref span = 3+10+1 = 14)
    // ref clip [103,113) should leave 6I10=6I with keep, 10= with excise.
    const std::string seq(26, 'A');
    std::vector<std::uint8_t> quals(26, 30);
    for (std::size_t i{0}; i < 26; ++i) {
        quals[i] = static_cast<std::uint8_t>(i);
    }

    BamRecord rec;
    rec.Name("flanking_fwd")
        .Flag(0)
        .RefId(0)
        .Pos(100)
        .MapQ(80)
        .Cigar(ParseCigar("3=6I10=6I1="))
        .Sequence(seq)
        .Qualities(quals);

    TagMap tags;
    tags.Set(TagKey{'q', 's'}, std::int64_t{500});
    tags.Set(TagKey{'q', 'e'}, std::int64_t{526});
    rec.Tags(std::move(tags));

    // With inserts kept
    {
        auto kept{rec};
        kept.Clip(ClipType::CLIP_TO_REFERENCE, 103, 113, false);
        EXPECT_EQ(CigarToString(kept.Cigar()), "6I10=6I");
        EXPECT_EQ(std::size(kept.Sequence()), 22U);
        EXPECT_EQ(kept.Pos(), 103);
        EXPECT_EQ(kept.ReferenceEnd(), 113);
    }

    // With inserts excised
    {
        auto excised{rec};
        excised.Clip(ClipType::CLIP_TO_REFERENCE, 103, 113, true);
        EXPECT_EQ(CigarToString(excised.Cigar()), "10=");
        EXPECT_EQ(std::size(excised.Sequence()), 10U);
        EXPECT_EQ(excised.Pos(), 103);
        EXPECT_EQ(excised.ReferenceEnd(), 113);
    }
}

TEST(BamRecordClipping, FlankingInserts_Reverse_KeepVsExcise)
{
    // Same CIGAR reverse strand
    const std::string seq(26, 'A');
    const std::vector<std::uint8_t> quals(26, 30);

    BamRecord rec;
    rec.Name("flanking_rev")
        .Flag(0x10)  // reverse
        .RefId(0)
        .Pos(100)
        .MapQ(80)
        .Cigar(ParseCigar("3=6I10=6I1="))
        .Sequence(seq)
        .Qualities(quals);

    TagMap tags;
    tags.Set(TagKey{'q', 's'}, std::int64_t{500});
    tags.Set(TagKey{'q', 'e'}, std::int64_t{526});
    rec.Tags(std::move(tags));

    // With inserts kept
    {
        auto kept{rec};
        kept.Clip(ClipType::CLIP_TO_REFERENCE, 103, 113, false);
        EXPECT_EQ(CigarToString(kept.Cigar()), "6I10=6I");
        EXPECT_EQ(std::size(kept.Sequence()), 22U);
        EXPECT_EQ(kept.Pos(), 103);
    }

    // With inserts excised
    {
        auto excised{rec};
        excised.Clip(ClipType::CLIP_TO_REFERENCE, 103, 113, true);
        EXPECT_EQ(CigarToString(excised.Cigar()), "10=");
        EXPECT_EQ(std::size(excised.Sequence()), 10U);
        EXPECT_EQ(excised.Pos(), 103);
    }
}

TEST(BamRecordClipping, FlankingInserts_QueryClip_DoesNotExcise)
{
    // exciseFlankingInserts should be ignored for query clipping (pbbam behavior)
    const std::string seq(15, 'A');
    const std::vector<std::uint8_t> quals(15, 30);

    BamRecord rec;
    rec.Name("flanking_query")
        .Flag(0)
        .RefId(0)
        .Pos(100)
        .MapQ(80)
        .Cigar(ParseCigar("4I5=6I"))
        .Sequence(seq)
        .Qualities(quals);

    TagMap tags;
    tags.Set(TagKey{'q', 's'}, std::int64_t{500});
    tags.Set(TagKey{'q', 'e'}, std::int64_t{515});
    rec.Tags(std::move(tags));

    rec.Clip(ClipType::CLIP_TO_QUERY, 502, 512, true);

    // Query clipping is purely offset-based; flanking insert excision has no effect.
    EXPECT_EQ(CigarToString(rec.Cigar()), "2I5=3I");
    EXPECT_EQ(std::size(rec.Sequence()), 10U);
    EXPECT_EQ(rec.Pos(), 100);
    EXPECT_EQ(rec.ReferenceEnd(), 105);
}

// =====================================================================
// 8. CCS kinetics: fi/fp (forward) and ri/rp (reverse)
// =====================================================================

TEST(BamRecordClipping, CCSKinetics_FiRiFpRp)
{
    // fi/fp use SubstringClipStrategy, ri/rp use ReverseSubstringClipStrategy.
    // 10-base CCS, clip [2,7) -> 5 bases retained.
    BamRecord rec;
    rec.Name("movie/42/ccs")
        .Flag(0)
        .RefId(0)
        .Pos(100)
        .MapQ(60)
        .Cigar(ParseCigar("10="))
        .Sequence("AACCGTTAGC")
        .Qualities(std::vector<std::uint8_t>{10, 10, 20, 20, 30, 40, 40, 10, 30, 20});

    TagMap tags;

    // fi: {0,10,20,30,40,50,60,70,80,90} uint16
    TagArray fiArr{'S'};
    TagArray fpArr{'S'};
    TagArray riArr{'S'};
    TagArray rpArr{'S'};
    for (std::uint16_t i{0}; i < 10; ++i) {
        fiArr.AppendUInt16(i * 10);
        fpArr.AppendUInt16(i * 10 + 2);
        riArr.AppendUInt16(i * 10 + 4);
        rpArr.AppendUInt16(i * 10 + 6);
    }
    tags.Set(TagKey{'f', 'i'}, std::move(fiArr));
    tags.Set(TagKey{'f', 'p'}, std::move(fpArr));
    tags.Set(TagKey{'r', 'i'}, std::move(riArr));
    tags.Set(TagKey{'r', 'p'}, std::move(rpArr));
    rec.Tags(std::move(tags));

    rec.Clip(ClipType::CLIP_TO_QUERY, 2, 7);

    EXPECT_EQ(rec.Sequence(), "CCGTT");

    // fi/fp: SubstringClip [2..7) = indices 2,3,4,5,6
    const auto* fi{rec.Tags().Get(TagKey{'f', 'i'})};
    ASSERT_NE(fi, nullptr);
    const auto fiVals{ExtractUInt16(std::get<TagArray>(*fi))};
    EXPECT_EQ(fiVals, (std::vector<std::uint16_t>{20, 30, 40, 50, 60}));

    const auto* fp{rec.Tags().Get(TagKey{'f', 'p'})};
    ASSERT_NE(fp, nullptr);
    const auto fpVals{ExtractUInt16(std::get<TagArray>(*fp))};
    EXPECT_EQ(fpVals, (std::vector<std::uint16_t>{22, 32, 42, 52, 62}));

    // ri/rp: ReverseSubstringClip. reverseOffset = 10 - (2+5) = 3, so [3..8)
    const auto* ri{rec.Tags().Get(TagKey{'r', 'i'})};
    ASSERT_NE(ri, nullptr);
    const auto riVals{ExtractUInt16(std::get<TagArray>(*ri))};
    EXPECT_EQ(riVals, (std::vector<std::uint16_t>{34, 44, 54, 64, 74}));

    const auto* rp{rec.Tags().Get(TagKey{'r', 'p'})};
    ASSERT_NE(rp, nullptr);
    const auto rpVals{ExtractUInt16(std::get<TagArray>(*rp))};
    EXPECT_EQ(rpVals, (std::vector<std::uint16_t>{36, 46, 56, 66, 76}));
}

TEST(BamRecordClipping, CCSKinetics_SingleBase)
{
    BamRecord rec;
    rec.Name("movie/42/ccs")
        .Flag(0)
        .RefId(0)
        .Pos(100)
        .MapQ(60)
        .Cigar(ParseCigar("10="))
        .Sequence("AACCGTTAGC")
        .Qualities(std::vector<std::uint8_t>{10, 10, 20, 20, 30, 40, 40, 10, 30, 20});

    TagMap tags;
    TagArray fiArr{'S'};
    TagArray riArr{'S'};
    for (std::uint16_t i{0}; i < 10; ++i) {
        fiArr.AppendUInt16(i * 10);
        riArr.AppendUInt16(i * 10 + 4);
    }
    tags.Set(TagKey{'f', 'i'}, std::move(fiArr));
    tags.Set(TagKey{'r', 'i'}, std::move(riArr));
    rec.Tags(std::move(tags));

    rec.Clip(ClipType::CLIP_TO_QUERY, 6, 7);

    EXPECT_EQ(rec.Sequence(), "T");

    // fi: SubstringClip [6..7) = {60}
    const auto fiVals{ExtractUInt16(std::get<TagArray>(*rec.Tags().Get(TagKey{'f', 'i'})))};
    EXPECT_EQ(fiVals, (std::vector<std::uint16_t>{60}));

    // ri: ReverseSubstringClip. reverseOffset = 10 - (6+1) = 3, [3..4) = {34}
    const auto riVals{ExtractUInt16(std::get<TagArray>(*rec.Tags().Get(TagKey{'r', 'i'})))};
    EXPECT_EQ(riVals, (std::vector<std::uint16_t>{34}));
}

TEST(BamRecordClipping, CCSKinetics_NoClip)
{
    BamRecord rec;
    rec.Name("movie/42/ccs")
        .Flag(0)
        .RefId(0)
        .Pos(100)
        .MapQ(60)
        .Cigar(ParseCigar("10="))
        .Sequence("AACCGTTAGC")
        .Qualities(std::vector<std::uint8_t>{10, 10, 20, 20, 30, 40, 40, 10, 30, 20});

    TagMap tags;
    TagArray fiArr{'S'};
    TagArray riArr{'S'};
    for (std::uint16_t i{0}; i < 10; ++i) {
        fiArr.AppendUInt16(i * 10);
        riArr.AppendUInt16(i * 10 + 4);
    }
    tags.Set(TagKey{'f', 'i'}, std::move(fiArr));
    tags.Set(TagKey{'r', 'i'}, std::move(riArr));
    rec.Tags(std::move(tags));

    rec.Clip(ClipType::CLIP_TO_QUERY, 0, 10);

    EXPECT_EQ(rec.Sequence(), "AACCGTTAGC");
    const auto fiVals{ExtractUInt16(std::get<TagArray>(*rec.Tags().Get(TagKey{'f', 'i'})))};
    EXPECT_EQ(fiVals, (std::vector<std::uint16_t>{0, 10, 20, 30, 40, 50, 60, 70, 80, 90}));
    const auto riVals{ExtractUInt16(std::get<TagArray>(*rec.Tags().Get(TagKey{'r', 'i'})))};
    EXPECT_EQ(riVals, (std::vector<std::uint16_t>{4, 14, 24, 34, 44, 54, 64, 74, 84, 94}));
}

// =====================================================================
// 9. Basemods (MM/ML) through full BamRecord::Clip pipeline
// =====================================================================

TEST(BamRecordClipping, Basemods_MMML_FullPipeline)
{
    // Sequence: "ACTCCACGACTCGTCACACTCACGTCTCA" (29 bases)
    //  A C T C C A C G A C T  C  G  T  C  A  C  A  C  T  C  A  C  G  T  C  T  C  A
    //  0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28
    // C positions: 1, 3, 4, 6, 9, 11, 14, 16, 18, 20, 22, 25, 27
    //
    // MM: "C+m,3,1,4;" -> skip 3 Cs, mod at 4th C (C@6); skip 1 C, mod at next (C@11);
    //   skip 4 Cs (14,16,18,20), mod at C@22.
    // ML: {18, 128, 234} (3 mod sites: C@6, C@11, C@22)
    //
    // Clip [3,24) -> clipOffset=3, clipLength=21. Remaining: "CCACGACTCGTCACACTCACG"
    // C's in clipped seq (offset in clipped): 0,1,3,6,8,11,13,15,17,19
    // Mods from C@6 (now at clipped offset 3), C@11 (offset 8), C@22 (offset 19)
    // New MM: skip 2 Cs -> mod, skip 1 C -> mod, skip 4 Cs -> mod = "C+m,2,1,4;"
    // ML: {18, 128, 234} unchanged (all 3 mods retained)

    const std::string seq{"ACTCCACGACTCGTCACACTCACGTCTCA"};
    const std::vector<std::uint8_t> quals(29, 30);

    BamRecord rec;
    rec.Name("movie/42/ccs")
        .Flag(0)
        .RefId(0)
        .Pos(100)
        .MapQ(60)
        .Cigar(ParseCigar("29="))
        .Sequence(seq)
        .Qualities(quals);

    TagMap tags;
    tags.Set(TagKey{'M', 'M'}, std::string{"C+m,3,1,4;"});

    TagArray mlArr{'C'};
    mlArr.AppendUInt8(18);
    mlArr.AppendUInt8(128);
    mlArr.AppendUInt8(234);
    tags.Set(TagKey{'M', 'L'}, std::move(mlArr));
    rec.Tags(std::move(tags));

    rec.Clip(ClipType::CLIP_TO_QUERY, 3, 24);

    EXPECT_EQ(rec.Sequence(), "CCACGACTCGTCACACTCACG");
    EXPECT_EQ(std::size(rec.Sequence()), 21U);

    // Check MM tag
    const auto* mm{rec.Tags().Get(TagKey{'M', 'M'})};
    ASSERT_NE(mm, nullptr);
    const auto& mmStr{std::get<std::string>(*mm)};
    EXPECT_EQ(mmStr, "C+m,2,1,4;");

    // Check ML tag — all 3 mods retained
    const auto* ml{rec.Tags().Get(TagKey{'M', 'L'})};
    ASSERT_NE(ml, nullptr);
    const auto mlVals{ExtractUInt8(std::get<TagArray>(*ml))};
    EXPECT_EQ(mlVals, (std::vector<std::uint8_t>{18, 128, 234}));
}

TEST(BamRecordClipping, Basemods_MMML_LostAllMods)
{
    // Clip a region with no modification sites -> empty MM skips, empty ML
    const std::string seq{"ACTCCACGACTCGTCACACTCACGTCTCA"};
    const std::vector<std::uint8_t> quals(29, 30);

    BamRecord rec;
    rec.Name("movie/42/ccs")
        .Flag(0)
        .RefId(0)
        .Pos(100)
        .MapQ(60)
        .Cigar(ParseCigar("29="))
        .Sequence(seq)
        .Qualities(quals);

    TagMap tags;
    // MM: "C+m,3,1,4;" -> mods at C@6, C@11, C@22
    tags.Set(TagKey{'M', 'M'}, std::string{"C+m,3,1,4;"});
    TagArray mlArr{'C'};
    mlArr.AppendUInt8(18);
    mlArr.AppendUInt8(128);
    mlArr.AppendUInt8(234);
    tags.Set(TagKey{'M', 'L'}, std::move(mlArr));
    rec.Tags(std::move(tags));

    // Clip [1,5) — "CTCC" — C positions in original [1,5): C@1, C@3, C@4
    // None of these are mod sites (mods at C@6, C@11, C@22)
    rec.Clip(ClipType::CLIP_TO_QUERY, 1, 5);

    EXPECT_EQ(rec.Sequence(), "CTCC");

    const auto* mm{rec.Tags().Get(TagKey{'M', 'M'})};
    ASSERT_NE(mm, nullptr);
    EXPECT_EQ(std::get<std::string>(*mm), "C+m;");

    const auto* ml{rec.Tags().Get(TagKey{'M', 'L'})};
    ASSERT_NE(ml, nullptr);
    EXPECT_EQ(std::get<TagArray>(*ml).Count(), 0U);
}

// =====================================================================
// 10. Pileup (sa RLE tag) through full BamRecord::Clip pipeline
// =====================================================================

TEST(BamRecordClipping, Pileup_SaTag_FullPipeline)
{
    // Sequence: 29 bases, CIGAR: 29=
    const std::string seq{"ACTCCACGACTCGTCACACTCACGTCTCA"};
    const std::vector<std::uint8_t> quals(29, 30);

    BamRecord rec;
    rec.Name("movie/42/ccs")
        .Flag(0)
        .RefId(0)
        .Pos(100)
        .MapQ(60)
        .Cigar(ParseCigar("29="))
        .Sequence(seq)
        .Qualities(quals);

    TagMap tags;
    // sa: RLE pairs {runLen, cov}: {4, 3, 15, 4, 10, 2}
    // 4 bases at cov 3, 15 bases at cov 4, 10 bases at cov 2 -> total 29 bases (matches seqLen)
    TagArray saArr{'C'};
    saArr.AppendUInt8(4);   // run 1: length
    saArr.AppendUInt8(3);   // run 1: coverage
    saArr.AppendUInt8(15);  // run 2: length
    saArr.AppendUInt8(4);   // run 2: coverage
    saArr.AppendUInt8(10);  // run 3: length
    saArr.AppendUInt8(2);   // run 3: coverage
    tags.Set(TagKey{'s', 'a'}, std::move(saArr));

    // sm/sx: per-base arrays (29 elements)
    TagArray smArr{'C'};
    TagArray sxArr{'C'};
    for (std::uint8_t i{0}; i < 29; ++i) {
        smArr.AppendUInt8(i);
        sxArr.AppendUInt8(28 - i);
    }
    tags.Set(TagKey{'s', 'm'}, std::move(smArr));
    tags.Set(TagKey{'s', 'x'}, std::move(sxArr));

    rec.Tags(std::move(tags));

    // Clip [3,24) -> offset=3, length=21
    rec.Clip(ClipType::CLIP_TO_QUERY, 3, 24);

    EXPECT_EQ(rec.Sequence(), "CCACGACTCGTCACACTCACG");

    // sa: clip [3..24) from {4,3, 15,4, 10,2}, seqLen=29
    //   prefix = 3 bases. First run: 4@3. 3 < 4, split -> 1@3.
    //   suffix = 29-24 = 5 bases. Last run: 10@2. 5 < 10, split -> 5@2.
    //   Result: {1,3, 15,4, 5,2}
    const auto* sa{rec.Tags().Get(TagKey{'s', 'a'})};
    ASSERT_NE(sa, nullptr);
    const auto saVals{ExtractUInt8(std::get<TagArray>(*sa))};
    EXPECT_EQ(saVals, (std::vector<std::uint8_t>{1, 3, 15, 4, 5, 2}));

    // sm: clipped [3..24) = {3,4,...,23}
    const auto* sm{rec.Tags().Get(TagKey{'s', 'm'})};
    ASSERT_NE(sm, nullptr);
    const auto smVals{ExtractUInt8(std::get<TagArray>(*sm))};
    ASSERT_EQ(std::size(smVals), 21U);
    EXPECT_EQ(smVals[0], 3);
    EXPECT_EQ(smVals[20], 23);

    // sx: clipped [3..24) -> sx[3]=25, sx[23]=5
    const auto* sx{rec.Tags().Get(TagKey{'s', 'x'})};
    ASSERT_NE(sx, nullptr);
    const auto sxVals{ExtractUInt8(std::get<TagArray>(*sx))};
    ASSERT_EQ(std::size(sxVals), 21U);
    EXPECT_EQ(sxVals[0], 25);
    EXPECT_EQ(sxVals[20], 5);
}

TEST(BamRecordClipping, Pileup_SaTag_AlignedBoundary)
{
    // Clip exactly at run boundaries
    // sa: {5,3, 5,4, 5,2} -> 15 bases
    const std::string seq(15, 'A');
    const std::vector<std::uint8_t> quals(15, 30);

    BamRecord rec;
    rec.Name("pileup_boundary")
        .Flag(0)
        .RefId(0)
        .Pos(100)
        .MapQ(60)
        .Cigar(ParseCigar("15="))
        .Sequence(seq)
        .Qualities(quals);

    TagMap tags;
    TagArray saArr{'C'};
    saArr.AppendUInt8(5);
    saArr.AppendUInt8(3);
    saArr.AppendUInt8(5);
    saArr.AppendUInt8(4);
    saArr.AppendUInt8(5);
    saArr.AppendUInt8(2);
    tags.Set(TagKey{'s', 'a'}, std::move(saArr));
    rec.Tags(std::move(tags));

    // Clip [5,10) — exactly the middle run
    rec.Clip(ClipType::CLIP_TO_QUERY, 5, 10);

    const auto* sa{rec.Tags().Get(TagKey{'s', 'a'})};
    ASSERT_NE(sa, nullptr);
    const auto saVals{ExtractUInt8(std::get<TagArray>(*sa))};
    EXPECT_EQ(saVals, (std::vector<std::uint8_t>{5, 4}));
}

// =====================================================================
// 11. Serialization round-trip: clip -> serialize -> deserialize -> verify
// =====================================================================

namespace {
BamRecord RoundTrip(const BamRecord& rec)
{
    const std::vector<std::byte> bytes{rec.SerializeToBam()};
    const RawRecord raw{std::span<const std::byte>{bytes}};
    return raw.ToOwned();
}
}  // namespace

TEST(BamRecordClipping, RoundTrip_SimpleForwardClip)
{
    // Create record with tags, clip to query, serialize, deserialize, verify all fields
    auto rec{MakeRecordWithTags()};
    rec.Clip(ClipType::CLIP_TO_QUERY, 102, 108);

    const auto rt{RoundTrip(rec)};

    // Core fields
    EXPECT_EQ(rt.Name(), rec.Name());
    EXPECT_EQ(rt.Flag(), rec.Flag());
    EXPECT_EQ(rt.RefId(), rec.RefId());
    EXPECT_EQ(rt.Pos(), rec.Pos());
    EXPECT_EQ(rt.MapQ(), rec.MapQ());
    EXPECT_EQ(CigarToString(rt.Cigar()), CigarToString(rec.Cigar()));
    EXPECT_EQ(rt.Sequence(), rec.Sequence());
    EXPECT_TRUE(std::ranges::equal(rt.Qualities(), rec.Qualities()));

    // Tags: qs/qe int64
    const auto* rtQs{rt.Tags().Get(TagKey{'q', 's'})};
    const auto* origQs{rec.Tags().Get(TagKey{'q', 's'})};
    ASSERT_NE(rtQs, nullptr);
    ASSERT_NE(origQs, nullptr);
    EXPECT_EQ(std::get<std::int64_t>(*rtQs), std::get<std::int64_t>(*origQs));

    const auto* rtQe{rt.Tags().Get(TagKey{'q', 'e'})};
    const auto* origQe{rec.Tags().Get(TagKey{'q', 'e'})};
    ASSERT_NE(rtQe, nullptr);
    ASSERT_NE(origQe, nullptr);
    EXPECT_EQ(std::get<std::int64_t>(*rtQe), std::get<std::int64_t>(*origQe));

    // Tags: ip array
    const auto* rtIp{rt.Tags().Get(TagKey{'i', 'p'})};
    const auto* origIp{rec.Tags().Get(TagKey{'i', 'p'})};
    ASSERT_NE(rtIp, nullptr);
    ASSERT_NE(origIp, nullptr);
    EXPECT_EQ(std::get<TagArray>(*rtIp), std::get<TagArray>(*origIp));
}

TEST(BamRecordClipping, RoundTrip_ReferenceClipWithDeletion)
{
    // Record with 3M2D5M CIGAR, clip to reference, round-trip, verify.
    // Use valid nucleotides (BAM 4-bit encoding only supports =ACMGRSVTWYHKDBN).
    BamRecord rec;
    rec.Name("read3")
        .Flag(0)
        .RefId(0)
        .Pos(50)
        .MapQ(60)
        .Cigar(ParseCigar("3M2D5M"))
        .Sequence("ACGTAGCT")
        .Qualities(std::vector<std::uint8_t>{10, 20, 30, 40, 50, 60, 70, 80});

    TagMap tags;
    tags.Set(TagKey{'q', 's'}, std::int64_t{0});
    tags.Set(TagKey{'q', 'e'}, std::int64_t{8});
    tags.Set(TagKey{'R', 'G'}, std::string{"test_group"});
    rec.Tags(std::move(tags));

    // ref clip [54,58): front removes 2 ref from 3M -> 1M left, then 2D, back removes 2 from 5M
    // Clipped CIGAR: 1M2D3M (1 ref from tail of 3M + 2D + 3 ref from 5M = 6 ref span = 58-52? no)
    // Actually: refSpan = 3+2+5 = 10, pos=50, refEnd=60.
    // Front: remove 54-50 = 4 ref from front. Walk 3M(3 ref, 3 query consumed), then 2D(2 ref,
    //   0 query). 3+2=5 > 4. Back up: after 3M we have 3 consumed, need 1 more from 2D -> 1D.
    //   Remaining from front walk: 1D + 5M, queryRemovedFront=3.
    // Back: remove 60-58 = 2 ref from back. Walk 5M -> 2 consumed from back, 3M left.
    //   Remaining: 1D + 3M, queryRemovedBack=2.
    // clipOffset=3, clipLength=8-3-2=3, newPos=54.
    rec.Clip(ClipType::CLIP_TO_REFERENCE, 54, 58);

    EXPECT_EQ(CigarToString(rec.Cigar()), "1D3M");
    EXPECT_EQ(rec.Pos(), 54);
    EXPECT_EQ(rec.Sequence(), "TAG");

    const auto rt{RoundTrip(rec)};

    // Core fields
    EXPECT_EQ(rt.Name(), rec.Name());
    EXPECT_EQ(rt.Flag(), rec.Flag());
    EXPECT_EQ(rt.RefId(), rec.RefId());
    EXPECT_EQ(rt.Pos(), rec.Pos());
    EXPECT_EQ(rt.MapQ(), rec.MapQ());
    EXPECT_EQ(CigarToString(rt.Cigar()), "1D3M");
    EXPECT_EQ(rt.Sequence(), "TAG");
    EXPECT_TRUE(std::ranges::equal(rt.Qualities(), rec.Qualities()));

    // String tag survives round-trip
    const auto* rtRg{rt.Tags().Get(TagKey{'R', 'G'})};
    ASSERT_NE(rtRg, nullptr);
    EXPECT_EQ(std::get<std::string>(*rtRg), "test_group");

    // Int64 tags survive
    const auto* rtQs{rt.Tags().Get(TagKey{'q', 's'})};
    const auto* origQs{rec.Tags().Get(TagKey{'q', 's'})};
    ASSERT_NE(rtQs, nullptr);
    ASSERT_NE(origQs, nullptr);
    EXPECT_EQ(std::get<std::int64_t>(*rtQs), std::get<std::int64_t>(*origQs));
}

TEST(BamRecordClipping, RoundTrip_CCSKinetics)
{
    // Record with fi/fp/ri/rp, clip, round-trip, verify tag arrays survive
    BamRecord rec;
    rec.Name("movie/42/ccs")
        .Flag(0)
        .RefId(0)
        .Pos(100)
        .MapQ(60)
        .Cigar(ParseCigar("10="))
        .Sequence("AACCGTTAGC")
        .Qualities(std::vector<std::uint8_t>{10, 10, 20, 20, 30, 40, 40, 10, 30, 20});

    TagMap tags;
    TagArray fiArr{'S'};
    TagArray fpArr{'S'};
    TagArray riArr{'S'};
    TagArray rpArr{'S'};
    for (std::uint16_t i{0}; i < 10; ++i) {
        fiArr.AppendUInt16(i * 10);
        fpArr.AppendUInt16(i * 10 + 2);
        riArr.AppendUInt16(i * 10 + 4);
        rpArr.AppendUInt16(i * 10 + 6);
    }
    tags.Set(TagKey{'f', 'i'}, std::move(fiArr));
    tags.Set(TagKey{'f', 'p'}, std::move(fpArr));
    tags.Set(TagKey{'r', 'i'}, std::move(riArr));
    tags.Set(TagKey{'r', 'p'}, std::move(rpArr));
    rec.Tags(std::move(tags));

    rec.Clip(ClipType::CLIP_TO_QUERY, 2, 7);

    EXPECT_EQ(rec.Sequence(), "CCGTT");

    const auto rt{RoundTrip(rec)};

    // Core fields
    EXPECT_EQ(rt.Name(), rec.Name());
    EXPECT_EQ(rt.Pos(), rec.Pos());
    EXPECT_EQ(CigarToString(rt.Cigar()), "5=");
    EXPECT_EQ(rt.Sequence(), "CCGTT");
    EXPECT_TRUE(std::ranges::equal(rt.Qualities(), rec.Qualities()));

    // fi: SubstringClip [2..7) = {20, 30, 40, 50, 60}
    const auto* rtFi{rt.Tags().Get(TagKey{'f', 'i'})};
    ASSERT_NE(rtFi, nullptr);
    const auto rtFiVals{ExtractUInt16(std::get<TagArray>(*rtFi))};
    EXPECT_EQ(rtFiVals, (std::vector<std::uint16_t>{20, 30, 40, 50, 60}));

    // fp: {22, 32, 42, 52, 62}
    const auto* rtFp{rt.Tags().Get(TagKey{'f', 'p'})};
    ASSERT_NE(rtFp, nullptr);
    const auto rtFpVals{ExtractUInt16(std::get<TagArray>(*rtFp))};
    EXPECT_EQ(rtFpVals, (std::vector<std::uint16_t>{22, 32, 42, 52, 62}));

    // ri: ReverseSubstringClip. reverseOffset = 10 - (2+5) = 3, [3..8) = {34, 44, 54, 64, 74}
    const auto* rtRi{rt.Tags().Get(TagKey{'r', 'i'})};
    ASSERT_NE(rtRi, nullptr);
    const auto rtRiVals{ExtractUInt16(std::get<TagArray>(*rtRi))};
    EXPECT_EQ(rtRiVals, (std::vector<std::uint16_t>{34, 44, 54, 64, 74}));

    // rp: {36, 46, 56, 66, 76}
    const auto* rtRp{rt.Tags().Get(TagKey{'r', 'p'})};
    ASSERT_NE(rtRp, nullptr);
    const auto rtRpVals{ExtractUInt16(std::get<TagArray>(*rtRp))};
    EXPECT_EQ(rtRpVals, (std::vector<std::uint16_t>{36, 46, 56, 66, 76}));

    // Verify all arrays match the original clipped record exactly
    EXPECT_EQ(std::get<TagArray>(*rtFi), std::get<TagArray>(*rec.Tags().Get(TagKey{'f', 'i'})));
    EXPECT_EQ(std::get<TagArray>(*rtFp), std::get<TagArray>(*rec.Tags().Get(TagKey{'f', 'p'})));
    EXPECT_EQ(std::get<TagArray>(*rtRi), std::get<TagArray>(*rec.Tags().Get(TagKey{'r', 'i'})));
    EXPECT_EQ(std::get<TagArray>(*rtRp), std::get<TagArray>(*rec.Tags().Get(TagKey{'r', 'p'})));
}

TEST(BamRecordClipping, RoundTrip_BasemodsMMML)
{
    // Record with MM/ML, clip, round-trip, verify MM string and ML array survive
    const std::string seq{"ACTCCACGACTCGTCACACTCACGTCTCA"};
    const std::vector<std::uint8_t> quals(29, 30);

    BamRecord rec;
    rec.Name("movie/42/ccs")
        .Flag(0)
        .RefId(0)
        .Pos(100)
        .MapQ(60)
        .Cigar(ParseCigar("29="))
        .Sequence(seq)
        .Qualities(quals);

    TagMap tags;
    tags.Set(TagKey{'M', 'M'}, std::string{"C+m,3,1,4;"});
    TagArray mlArr{'C'};
    mlArr.AppendUInt8(18);
    mlArr.AppendUInt8(128);
    mlArr.AppendUInt8(234);
    tags.Set(TagKey{'M', 'L'}, std::move(mlArr));
    rec.Tags(std::move(tags));

    // Clip [3,24) -> retains all 3 mod sites (C@6, C@11, C@22)
    rec.Clip(ClipType::CLIP_TO_QUERY, 3, 24);

    EXPECT_EQ(rec.Sequence(), "CCACGACTCGTCACACTCACG");

    const auto rt{RoundTrip(rec)};

    // Core fields
    EXPECT_EQ(rt.Name(), rec.Name());
    EXPECT_EQ(rt.Pos(), rec.Pos());
    EXPECT_EQ(CigarToString(rt.Cigar()), "21=");
    EXPECT_EQ(rt.Sequence(), "CCACGACTCGTCACACTCACG");
    EXPECT_TRUE(std::ranges::equal(rt.Qualities(), rec.Qualities()));

    // MM string tag survives round-trip
    const auto* rtMm{rt.Tags().Get(TagKey{'M', 'M'})};
    ASSERT_NE(rtMm, nullptr);
    EXPECT_EQ(std::get<std::string>(*rtMm), "C+m,2,1,4;");

    // ML array tag survives round-trip — all 3 mods retained
    const auto* rtMl{rt.Tags().Get(TagKey{'M', 'L'})};
    ASSERT_NE(rtMl, nullptr);
    const auto rtMlVals{ExtractUInt8(std::get<TagArray>(*rtMl))};
    EXPECT_EQ(rtMlVals, (std::vector<std::uint8_t>{18, 128, 234}));

    // Verify ML array matches original exactly
    EXPECT_EQ(std::get<TagArray>(*rtMl), std::get<TagArray>(*rec.Tags().Get(TagKey{'M', 'L'})));
}

}  // namespace Samoa
}  // namespace PacBio
