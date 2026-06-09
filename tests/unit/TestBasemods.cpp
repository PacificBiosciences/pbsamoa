#include <pbsamoa/core/Basemods.hpp>

#include <gtest/gtest.h>

#include <string>
#include <vector>

namespace PacBio {
namespace Samoa {

TEST(Basemods, ChebiNumericModCodeParsed)
{
    // A numeric ChEBI mod code is a single code, not a base followed by a skip count
    // (htslib sam_mods.c). "C+76792,5,12" must keep the digits in the prefix.
    const std::vector<BasemodRecord> records{ParseBasemodString("C+76792,5,12;")};
    ASSERT_EQ(std::size(records), 1U);
    EXPECT_EQ(records[0].Prefix, "C+76792");
    ASSERT_EQ(std::size(records[0].Skips), 2U);
    EXPECT_EQ(records[0].Skips[0], 5);
    EXPECT_EQ(records[0].Skips[1], 12);
}

TEST(Basemods, ChebiModCodeRoundTrips)
{
    const std::vector<BasemodRecord> records{ParseBasemodString("C+76792,5,12;")};
    EXPECT_EQ(WriteBasemodString(records), "C+76792,5,12;");
}

TEST(Basemods, SingleLetterModCodeStillParses)
{
    const std::vector<BasemodRecord> records{ParseBasemodString("C+m?,1,3,0;A+a,2;")};
    ASSERT_EQ(std::size(records), 2U);
    EXPECT_EQ(records[0].Prefix, "C+m?");
    EXPECT_EQ(records[0].Skips, (std::vector<std::int32_t>{1, 3, 0}));
    EXPECT_EQ(records[1].Prefix, "A+a");
    EXPECT_EQ(records[1].Skips, (std::vector<std::int32_t>{2}));
}

TEST(Basemods, MultiModCodeParsedAndRoundTrips)
{
    // "C+mh" carries two modifications (m, h) sharing one delta list.
    const std::vector<BasemodRecord> records{ParseBasemodString("C+mh,5,12;")};
    ASSERT_EQ(std::size(records), 1U);
    EXPECT_EQ(records[0].Prefix, "C+mh");
    EXPECT_EQ(records[0].Skips, (std::vector<std::int32_t>{5, 12}));
    EXPECT_EQ(WriteBasemodString(records), "C+mh,5,12;");
}

TEST(Basemods, ModCodeCountReflectsStride)
{
    EXPECT_EQ(ModCodeCount("C+m"), 1U);
    EXPECT_EQ(ModCodeCount("C+m?"), 1U);
    EXPECT_EQ(ModCodeCount("C+mh"), 2U);
    EXPECT_EQ(ModCodeCount("C+mhf."), 3U);
    EXPECT_EQ(ModCodeCount("C+76792"), 1U);  // numeric ChEBI code = single modification
}

TEST(Basemods, NBaseModCountsAllBases)
{
    // The 'N' canonical base is a wildcard matching every base, so on a pure-ACGT read the
    // mod sites must still be located. Mods at the 1st/2nd/3rd base ("N+n,0,0,0").
    BasemodRecord rec;
    rec.Prefix = "N+n";
    rec.Skips = {0, 0, 0};

    // Clip away the first two bases; retain a 4-base window.
    const BasemodClipWindow window{ClipBasemodRecord(rec, "ACGTACGT", /*clipOffset=*/2,
                                                     /*clipLength=*/4)};
    EXPECT_EQ(window.FrontRemoved, 2U);  // sites at bases 0 and 1 fall before the window
    EXPECT_EQ(window.Retained, 1U);      // site at base 2 is retained
}

TEST(Basemods, ReverseStrandClipMirrorsWindowAndComplementsBase)
{
    // A reverse-strand read stores SEQ as the reverse-complement of the original 5'->3'
    // read, while MM coordinates stay in original orientation. Original read "CCCAAA" has
    // three modified C's (C+m,0,0,0); the stored SEQ is revcomp("CCCAAA") = "TTTGGG".
    // Keeping SEQ[3:6] ("GGG") keeps exactly the original "CCC", so all three calls survive.
    BasemodRecord rec;
    rec.Prefix = "C+m";
    rec.Skips = {0, 0, 0};

    // Walking the stored SEQ as-if-forward finds no literal 'C' in "GGG" and wrongly drops
    // every call — this is the pre-fix behaviour and why reverse handling is required.
    const BasemodClipWindow fwd{
        ClipBasemodRecord(rec, "TTTGGG", /*clipOffset=*/3, /*clipLength=*/3, /*reverse=*/false)};
    EXPECT_EQ(fwd.Retained, 0U);

    const BasemodClipWindow rev{
        ClipBasemodRecord(rec, "TTTGGG", /*clipOffset=*/3, /*clipLength=*/3, /*reverse=*/true)};
    EXPECT_EQ(rev.FrontRemoved, 0U);
    EXPECT_EQ(rev.Retained, 3U);
    EXPECT_EQ(rev.RetainedSkips, (std::vector<std::int32_t>{0, 0, 0}));
}

TEST(Basemods, ReverseStrandClipFrontRemovesTrailingOriginalCalls)
{
    // Same reverse read; keeping SEQ[0:3] ("TTT") keeps the original trailing "AAA", which
    // carries no C, so all three C+m calls fall before the (mirrored) window and are removed.
    BasemodRecord rec;
    rec.Prefix = "C+m";
    rec.Skips = {0, 0, 0};
    const BasemodClipWindow rev{
        ClipBasemodRecord(rec, "TTTGGG", /*clipOffset=*/0, /*clipLength=*/3, /*reverse=*/true)};
    EXPECT_EQ(rev.FrontRemoved, 3U);
    EXPECT_EQ(rev.Retained, 0U);
}

}  // namespace Samoa
}  // namespace PacBio
