#include <pbsamoa/core/SamHeader.hpp>

#include <gtest/gtest.h>

#include <string>

namespace PacBio {
namespace Samoa {

TEST(MakeReadGroupId, BasicMovieAndReadType)
{
    const std::string id{MakeReadGroupId("movie1", "SUBREAD")};
    EXPECT_EQ(std::size(id), 8U);
    // Same inputs always produce the same ID.
    EXPECT_EQ(id, MakeReadGroupId("movie1", "SUBREAD"));
}

TEST(MakeReadGroupId, DifferentInputsProduceDifferentIds)
{
    EXPECT_NE(MakeReadGroupId("movie1", "SUBREAD"), MakeReadGroupId("movie2", "SUBREAD"));
    EXPECT_NE(MakeReadGroupId("movie1", "SUBREAD"), MakeReadGroupId("movie1", "CCS"));
}

TEST(MakeReadGroupId, StrandAffectsId)
{
    const std::string unstrand{MakeReadGroupId("movie1", "SUBREAD")};
    const std::string fwd{MakeReadGroupId("movie1", "SUBREAD", Data::Strand::FORWARD)};
    const std::string rev{MakeReadGroupId("movie1", "SUBREAD", Data::Strand::REVERSE)};
    EXPECT_NE(unstrand, fwd);
    EXPECT_NE(unstrand, rev);
    EXPECT_NE(fwd, rev);
}

TEST(MakeReadGroupId, UnmappedStrandSameAsNoStrand)
{
    EXPECT_EQ(MakeReadGroupId("movie1", "CCS"),
              MakeReadGroupId("movie1", "CCS", Data::Strand::UNMAPPED));
}

TEST(MakeReadGroupId, WithBarcode)
{
    const std::string id{MakeReadGroupId("movie1", "SUBREAD", "0--0")};
    const std::string base{MakeReadGroupId("movie1", "SUBREAD")};
    EXPECT_EQ(id, base + "/0--0");
}

TEST(ReadGroupBaseId, StripsBarcodes)
{
    EXPECT_EQ(ReadGroupBaseId("abcd1234/0--0"), "abcd1234");
    EXPECT_EQ(ReadGroupBaseId("abcd1234"), "abcd1234");
}

// --- ReadGroup DS field accessors ---

TEST(ReadGroupDs, DsFieldAccess)
{
    ReadGroup rg{"test"};
    rg.SetTag("DS", "READTYPE=SUBREAD;BINDINGKIT=101;SEQUENCINGKIT=202");

    ASSERT_NE(rg.ReadType(), nullptr);
    EXPECT_EQ(*rg.ReadType(), "SUBREAD");
    ASSERT_NE(rg.BindingKit(), nullptr);
    EXPECT_EQ(*rg.BindingKit(), "101");
    ASSERT_NE(rg.SequencingKit(), nullptr);
    EXPECT_EQ(*rg.SequencingKit(), "202");
    EXPECT_EQ(rg.BasecallerVersion(), nullptr);
}

TEST(ReadGroupDs, DsFieldGenericAccess)
{
    ReadGroup rg{"test"};
    rg.SetTag("DS", "READTYPE=SUBREAD;CUSTOM=value");

    ASSERT_NE(rg.DsField("CUSTOM"), nullptr);
    EXPECT_EQ(*rg.DsField("CUSTOM"), "value");
    EXPECT_EQ(rg.DsField("MISSING"), nullptr);
}

TEST(ReadGroupDs, SetDsFieldUpdatesExisting)
{
    ReadGroup rg{"test"};
    rg.SetTag("DS", "READTYPE=SUBREAD;BINDINGKIT=101");

    rg.SetDsField("READTYPE", "CCS");
    ASSERT_NE(rg.ReadType(), nullptr);
    EXPECT_EQ(*rg.ReadType(), "CCS");
    // Other fields preserved
    ASSERT_NE(rg.BindingKit(), nullptr);
    EXPECT_EQ(*rg.BindingKit(), "101");
}

TEST(ReadGroupDs, SetDsFieldInsertsNew)
{
    ReadGroup rg{"test"};
    rg.SetTag("DS", "READTYPE=SUBREAD");

    rg.SetDsField("STRAND", "FORWARD");
    ASSERT_NE(rg.DsField("STRAND"), nullptr);
    EXPECT_EQ(*rg.DsField("STRAND"), "FORWARD");
}

TEST(ReadGroupDs, NoDsTagReturnsNullptr)
{
    const ReadGroup rg{"test"};
    EXPECT_EQ(rg.ReadType(), nullptr);
    EXPECT_EQ(rg.BindingKit(), nullptr);
    EXPECT_EQ(rg.DsField("ANYTHING"), nullptr);
}

TEST(ReadGroupDs, ParsedDsFieldsRoundtrip)
{
    ReadGroup rg{"test"};
    rg.SetTag("DS", "A=1;B=2;C=3");

    const auto fields{rg.ParsedDsFields()};
    ASSERT_EQ(std::size(fields), 3U);
    EXPECT_EQ(fields[0].first, "A");
    EXPECT_EQ(fields[0].second, "1");

    rg.SetDsFields(fields);
    const std::string* dsTag{rg.GetTag("DS")};
    ASSERT_NE(dsTag, nullptr);
    EXPECT_EQ(*dsTag, "A=1;B=2;C=3");
}

TEST(ReadGroupDs, ParsesWhitespaceAndTrailingDelimiter)
{
    ReadGroup rg{"test"};
    rg.SetTag("DS", " READTYPE = SUBREAD ; CUSTOM = value ;");

    const auto fields{rg.ParsedDsFields()};
    ASSERT_EQ(std::size(fields), 2U);
    EXPECT_EQ(fields[0].first, "READTYPE");
    EXPECT_EQ(fields[0].second, "SUBREAD");
    EXPECT_EQ(fields[1].first, "CUSTOM");
    EXPECT_EQ(fields[1].second, "value");
}

TEST(ReadGroupDs, MovieNameFromPuTag)
{
    ReadGroup rg{"test"};
    EXPECT_EQ(rg.MovieName(), nullptr);

    rg.SetTag("PU", "m84001_230101_000000");
    ASSERT_NE(rg.MovieName(), nullptr);
    EXPECT_EQ(*rg.MovieName(), "m84001_230101_000000");
}

TEST(ReadGroupDs, SetTagDsInvalidatesCache)
{
    ReadGroup rg{"test"};
    rg.SetTag("DS", "READTYPE=SUBREAD");
    ASSERT_NE(rg.ReadType(), nullptr);
    EXPECT_EQ(*rg.ReadType(), "SUBREAD");

    // Overwrite DS tag directly — cache should be invalidated
    rg.SetTag("DS", "READTYPE=CCS");
    ASSERT_NE(rg.ReadType(), nullptr);
    EXPECT_EQ(*rg.ReadType(), "CCS");
}

}  // namespace Samoa
}  // namespace PacBio
