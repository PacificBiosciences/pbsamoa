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

TEST(ParseUtils, ParseRegionRejectsInvalidFormat)
{
    const auto region{ParseRegion("chr1", "ref:start-end")};
    ASSERT_FALSE(region);
    EXPECT_EQ(region.error(), "invalid region format 'chr1' (expected ref:start-end)");
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

TEST(ParseUtils, ParseRegionRejectsInvalidCoordinateRange)
{
    const auto region{ParseRegion("chr1:0-9", "ref:start-end")};
    ASSERT_FALSE(region);
    EXPECT_EQ(region.error(), "region must satisfy start >= 1 and end >= start");
}

}  // namespace Tools
}  // namespace Samoa
}  // namespace PacBio
