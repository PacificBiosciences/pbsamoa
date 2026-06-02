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

}  // namespace Samoa
}  // namespace PacBio
