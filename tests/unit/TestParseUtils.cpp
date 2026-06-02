#include "../../src/tools/ParseUtils.hpp"

#include <gtest/gtest.h>

namespace PacBio {
namespace Samoa {
namespace Tools {

TEST(ParseUtils, ParseRegionConvertsToZeroBasedCoordinates)
{
    const auto region{ParseRegion("chr1:5-9", "ref:start-end")};
    ASSERT_TRUE(region);
    EXPECT_EQ(region->RefName, "chr1");
    EXPECT_EQ(region->Beg, 4);
    EXPECT_EQ(region->End, 9);
}

TEST(ParseUtils, ParseRegionSplitsOnRightmostColonForColonRefNames)
{
    // Reference names may contain ':' (e.g. HLA alleles); the rightmost colon delimits.
    const auto region{ParseRegion("HLA-DRB1*12:17:100-200", "ref:start-end")};
    ASSERT_TRUE(region);
    EXPECT_EQ(region->RefName, "HLA-DRB1*12:17");
    EXPECT_EQ(region->Beg, 99);
    EXPECT_EQ(region->End, 200);
}

TEST(ParseUtils, ParseRegionBareNameIsWholeReference)
{
    // A bare reference name selects the whole reference (htslib accepts this).
    const auto region{ParseRegion("chr1", "ref:start-end")};
    ASSERT_TRUE(region);
    EXPECT_EQ(region->RefName, "chr1");
    EXPECT_EQ(region->Beg, 0);
    EXPECT_EQ(region->End, std::numeric_limits<std::int32_t>::max());
}

TEST(ParseUtils, ParseRegionOpenEndedForms)
{
    // "chr:100" and "chr:100-" run from the start coordinate to the end of the reference.
    for (const std::string_view text : {"chr1:100", "chr1:100-"}) {
        const auto region{ParseRegion(text, "ref:start-end")};
        ASSERT_TRUE(region) << text;
        EXPECT_EQ(region->Beg, 99);
        EXPECT_EQ(region->End, std::numeric_limits<std::int32_t>::max());
    }

    // "chr:-200" runs from the start of the reference to coordinate 200.
    const auto fromStart{ParseRegion("chr1:-200", "ref:start-end")};
    ASSERT_TRUE(fromStart);
    EXPECT_EQ(fromStart->Beg, 0);
    EXPECT_EQ(fromStart->End, 200);
}

TEST(ParseUtils, ParseRegionRejectsInvalidStartValue)
{
    const auto region{ParseRegion("chr1:foo-9", "ref:start-end")};
    ASSERT_FALSE(region);
    EXPECT_EQ(region.error(), "invalid region start: foo");
}

TEST(ParseUtils, ParseRegionRejectsInvalidEndValue)
{
    const auto region{ParseRegion("chr1:5-foo", "ref:start-end")};
    ASSERT_FALSE(region);
    EXPECT_EQ(region.error(), "invalid region end: foo");
}

TEST(ParseUtils, ParseRegionRejectsEndBeforeStart)
{
    const auto region{ParseRegion("chr1:9-5", "ref:start-end")};
    ASSERT_FALSE(region);
    EXPECT_EQ(region.error(), "region end must be >= start");
}

}  // namespace Tools
}  // namespace Samoa
}  // namespace PacBio
