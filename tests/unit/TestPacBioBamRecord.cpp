#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <string>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {
namespace {

BamRecord MakePacBioRecord(std::string name)
{
    BamRecord record;
    record.Name(std::move(name));
    return record;
}

BamRecord MakeFullRecord()
{
    BamRecord record;
    record.Name("movie1/12345/100_2000");
    record.MutableTags().Set(RG_TAG, std::string{"abcd1234"});
    record.MutableTags().Set(TagKey{'q', 's'}, std::int64_t{100});
    record.MutableTags().Set(TagKey{'q', 'e'}, std::int64_t{2000});

    // SNR: 4 floats
    TagArray snrArray{'f'};
    snrArray.AppendFloat(3.0f);
    snrArray.AppendFloat(4.0f);
    snrArray.AppendFloat(2.5f);
    snrArray.AppendFloat(5.0f);
    record.MutableTags().Set(TagKey{'s', 'n'}, std::move(snrArray));

    // Read accuracy
    record.MutableTags().Set(TagKey{'r', 'q'}, 0.95f);

    // Context flags
    record.MutableTags().Set(TagKey{'c', 'x'}, std::int64_t{3});

    // Wall times
    record.MutableTags().Set(TagKey{'w', 's'}, std::int64_t{1000});
    record.MutableTags().Set(TagKey{'w', 'e'}, std::int64_t{2000});

    return record;
}

}  // namespace

TEST(PacBioBamRecord, MovieName)
{
    const auto record{MakePacBioRecord("movie1/12345/100_2000")};
    EXPECT_EQ(record.MovieName(), "movie1");
}

TEST(PacBioBamRecord, HoleNumber)
{
    const auto record{MakePacBioRecord("movie1/12345/100_2000")};
    EXPECT_EQ(record.HoleNumber(), 12345);
}

TEST(PacBioBamRecord, QueryStartFromTag)
{
    const auto record{MakeFullRecord()};
    EXPECT_EQ(record.QueryStart(), 100);
}

TEST(PacBioBamRecord, QueryEndFromTag)
{
    const auto record{MakeFullRecord()};
    EXPECT_EQ(record.QueryEnd(), 2000);
}

TEST(PacBioBamRecord, QueryStartFromName)
{
    const auto record{MakePacBioRecord("movie1/12345/100_2000")};
    EXPECT_EQ(record.QueryStart(), 100);
}

TEST(PacBioBamRecord, QueryEndFromName)
{
    const auto record{MakePacBioRecord("movie1/12345/100_2000")};
    EXPECT_EQ(record.QueryEnd(), 2000);
}

TEST(PacBioBamRecord, ReadGroupId)
{
    const auto record{MakeFullRecord()};
    EXPECT_EQ(record.ReadGroupId(), "abcd1234");
}

TEST(PacBioBamRecord, ReadGroupIdMissingThrows)
{
    const auto record{MakePacBioRecord("movie1/12345/100_2000")};
    EXPECT_THROW(record.ReadGroupId(), std::runtime_error);
}

TEST(PacBioBamRecord, LocalContextFlagsPresent)
{
    const auto record{MakeFullRecord()};
    const auto flags{record.LocalContextFlags()};
    ASSERT_TRUE(flags);
    EXPECT_EQ(static_cast<int>(*flags), 3);
}

TEST(PacBioBamRecord, LocalContextFlagsAbsent)
{
    const auto record{MakePacBioRecord("movie1/12345/100_2000")};
    EXPECT_FALSE(record.LocalContextFlags());
}

TEST(PacBioBamRecord, LocalContextFlagsSetter)
{
    auto record{MakeFullRecord()};
    record.LocalContextFlags(static_cast<Data::LocalContextFlags>(7));
    const auto flags{record.LocalContextFlags()};
    ASSERT_TRUE(flags);
    EXPECT_EQ(static_cast<int>(*flags), 7);
}

TEST(PacBioBamRecord, SignalToNoise)
{
    const auto record{MakeFullRecord()};
    const auto snr{record.SignalToNoise()};
    EXPECT_FLOAT_EQ(snr[0], 3.0f);
    EXPECT_FLOAT_EQ(snr[1], 4.0f);
    EXPECT_FLOAT_EQ(snr[2], 2.5f);
    EXPECT_FLOAT_EQ(snr[3], 5.0f);
}

TEST(PacBioBamRecord, ReadAccuracy)
{
    const auto record{MakeFullRecord()};
    EXPECT_FLOAT_EQ(static_cast<float>(record.ReadAccuracy()), 0.95f);
}

TEST(PacBioBamRecord, PulseWidthAbsent)
{
    const auto record{MakePacBioRecord("movie1/12345/100_2000")};
    EXPECT_FALSE(record.PulseWidth());
}

TEST(PacBioBamRecord, IPDAbsent)
{
    const auto record{MakePacBioRecord("movie1/12345/100_2000")};
    EXPECT_FALSE(record.IPD());
}

TEST(PacBioBamRecord, WallStartPresent)
{
    const auto record{MakeFullRecord()};
    const auto ws{record.WallStart()};
    ASSERT_TRUE(ws);
    EXPECT_EQ(*ws, 1000);
}

TEST(PacBioBamRecord, WallEndPresent)
{
    const auto record{MakeFullRecord()};
    const auto we{record.WallEnd()};
    ASSERT_TRUE(we);
    EXPECT_EQ(*we, 2000);
}

TEST(PacBioBamRecord, WallStartAbsent)
{
    const auto record{MakePacBioRecord("movie1/12345/100_2000")};
    EXPECT_FALSE(record.WallStart());
}

TEST(PacBioBamRecord, WallEndAbsent)
{
    const auto record{MakePacBioRecord("movie1/12345/100_2000")};
    EXPECT_FALSE(record.WallEnd());
}

// Full-length CCS/HiFi/IsoSeq read names end in a non-interval segment ("ccs")
// and carry no qs/qe tag. QueryStart/QueryEnd must not throw — they return the
// whole-read interval [0, seqLen), matching pbbam semantics.
TEST(PacBioBamRecord, QueryStartCcsNameNoIntervalReturnsZero)
{
    auto record{MakePacBioRecord("movie1/12345/ccs")};
    record.Sequence("ACGTACGT");
    EXPECT_EQ(record.QueryStart(), 0);
}

TEST(PacBioBamRecord, QueryEndCcsNameNoIntervalReturnsSeqLength)
{
    auto record{MakePacBioRecord("movie1/12345/ccs")};
    record.Sequence("ACGTACGT");
    EXPECT_EQ(record.QueryEnd(), 8);
}

TEST(PacBioBamRecord, IsSecondaryReflectsFlag)
{
    BamRecord record;
    record.Flag(0x100);
    EXPECT_TRUE(record.IsSecondary());
    EXPECT_FALSE(record.IsSupplementary());
    EXPECT_FALSE(record.IsPrimary());
}

TEST(PacBioBamRecord, IsSupplementaryReflectsFlag)
{
    BamRecord record;
    record.Flag(0x800);
    EXPECT_TRUE(record.IsSupplementary());
    EXPECT_FALSE(record.IsSecondary());
    EXPECT_FALSE(record.IsPrimary());
}

// Forward strand: AlignedStart adds the leading soft-clip to QueryStart;
// AlignedEnd subtracts the trailing soft-clip from QueryEnd.
TEST(PacBioBamRecord, AlignedStartEndForwardStrandSoftClips)
{
    BamRecord record;
    record.Name("movie/1/0_20");  // QueryStart=0, QueryEnd=20 from name
    record.Sequence(std::string(20, 'A'));
    record.Flag(0x0);  // mapped, forward
    record.Cigar(*ParseCigar("3S10M7S"));
    EXPECT_EQ(record.AlignedStart(), 3);
    EXPECT_EQ(record.AlignedEnd(), 13);
}

// Reverse strand: the CIGAR is reference-oriented, so soft-clips are consumed
// from the opposite query end (polymerase coordinates).
TEST(PacBioBamRecord, AlignedStartEndReverseStrandClipsFromOppositeEnd)
{
    BamRecord record;
    record.Name("movie/1/0_20");
    record.Sequence(std::string(20, 'A'));
    record.Flag(0x10);  // mapped, reverse
    record.Cigar(*ParseCigar("3S10M7S"));
    EXPECT_EQ(record.AlignedStart(), 7);
    EXPECT_EQ(record.AlignedEnd(), 17);
}

// Map on a forward alignment sets alignment fields and leaves SEQ/QUAL intact.
TEST(PacBioBamRecord, MapForwardSetsFieldsAndKeepsSequence)
{
    BamRecord record;
    record.Name("movie/1/0_8");
    record.Flag(0x4);  // unmapped
    record.Sequence("ACGTACGT");
    record.Qualities({1, 2, 3, 4, 5, 6, 7, 8});
    record.Map(2, 100, false, *ParseCigar("8M"), 60);
    EXPECT_TRUE(record.IsMapped());
    EXPECT_FALSE(record.IsReverseStrand());
    EXPECT_EQ(record.RefId(), 2);
    EXPECT_EQ(record.Pos(), 100);
    EXPECT_EQ(record.MapQ(), 60u);
    EXPECT_EQ(record.Sequence(), "ACGTACGT");
    EXPECT_EQ(CigarToString(record.Cigar()), "8M");
}

// Map on a reverse alignment reverse-complements SEQ and reverses QUAL, matching
// the BAM convention (SEQ/QUAL stored in alignment orientation).
TEST(PacBioBamRecord, MapReverseReverseComplementsSeqAndReversesQual)
{
    BamRecord record;
    record.Name("movie/1/0_8");
    record.Flag(0x4);
    record.Sequence("AAACGGTT");
    record.Qualities({1, 2, 3, 4, 5, 6, 7, 8});
    record.Map(0, 50, true, *ParseCigar("8M"), 30);
    EXPECT_TRUE(record.IsMapped());
    EXPECT_TRUE(record.IsReverseStrand());
    EXPECT_EQ(record.Sequence(), "AACCGTTT");
    const std::vector<std::uint8_t> expectedQual{8, 7, 6, 5, 4, 3, 2, 1};
    EXPECT_TRUE(std::ranges::equal(record.Qualities(), expectedQual));
}

}  // namespace Samoa
}  // namespace PacBio
