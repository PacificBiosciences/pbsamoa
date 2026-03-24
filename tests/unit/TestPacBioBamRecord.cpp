#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <gtest/gtest.h>

#include <string>

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

TEST(PacBioBamRecord, FullName)
{
    const auto record{MakePacBioRecord("movie1/12345/100_2000")};
    EXPECT_EQ(record.FullName(), "movie1/12345/100_2000");
}

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
    ASSERT_TRUE(flags.has_value());
    EXPECT_EQ(static_cast<int>(*flags), 3);
}

TEST(PacBioBamRecord, LocalContextFlagsAbsent)
{
    const auto record{MakePacBioRecord("movie1/12345/100_2000")};
    EXPECT_FALSE(record.LocalContextFlags().has_value());
}

TEST(PacBioBamRecord, LocalContextFlagsSetter)
{
    auto record{MakeFullRecord()};
    record.LocalContextFlags(static_cast<Data::LocalContextFlags>(7));
    const auto flags{record.LocalContextFlags()};
    ASSERT_TRUE(flags.has_value());
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
    EXPECT_FALSE(record.PulseWidth().has_value());
}

TEST(PacBioBamRecord, IPDAbsent)
{
    const auto record{MakePacBioRecord("movie1/12345/100_2000")};
    EXPECT_FALSE(record.IPD().has_value());
}

TEST(PacBioBamRecord, WallStartPresent)
{
    const auto record{MakeFullRecord()};
    const auto ws{record.WallStart()};
    ASSERT_TRUE(ws.has_value());
    EXPECT_EQ(*ws, 1000);
}

TEST(PacBioBamRecord, WallEndPresent)
{
    const auto record{MakeFullRecord()};
    const auto we{record.WallEnd()};
    ASSERT_TRUE(we.has_value());
    EXPECT_EQ(*we, 2000);
}

TEST(PacBioBamRecord, WallStartAbsent)
{
    const auto record{MakePacBioRecord("movie1/12345/100_2000")};
    EXPECT_FALSE(record.WallStart().has_value());
}

TEST(PacBioBamRecord, WallEndAbsent)
{
    const auto record{MakePacBioRecord("movie1/12345/100_2000")};
    EXPECT_FALSE(record.WallEnd().has_value());
}

}  // namespace Samoa
}  // namespace PacBio
