#include "../../src/PulseBitset.hpp"

#include <gtest/gtest.h>

#include <string>

namespace PacBio {
namespace Samoa {

TEST(PulseBitset, ConstructFromPulseCall)
{
    const PulseBitset bs{"AaCgTt"};
    EXPECT_EQ(bs.Size(), 6U);
    EXPECT_EQ(bs.Count(), 3U);
    EXPECT_TRUE(bs.IsSet(0));
    EXPECT_FALSE(bs.IsSet(1));
    EXPECT_TRUE(bs.IsSet(2));
    EXPECT_FALSE(bs.IsSet(3));
    EXPECT_TRUE(bs.IsSet(4));
    EXPECT_FALSE(bs.IsSet(5));
}

TEST(PulseBitset, FindFirstAndNext)
{
    const PulseBitset bs{"aACgTt"};  // bits: 0,1,1,0,1,0
    EXPECT_EQ(bs.FindFirst(), 1U);
    EXPECT_EQ(bs.FindNext(1), 2U);
    EXPECT_EQ(bs.FindNext(2), 4U);
    EXPECT_EQ(bs.FindNext(4), PulseBitset::NPOS);
}

TEST(PulseBitset, FindFirst_empty)
{
    const PulseBitset bs{"aaa"};
    EXPECT_EQ(bs.FindFirst(), PulseBitset::NPOS);
}

TEST(PulseBitset, LargerThan64Bits)
{
    std::string pulses(64, 'a');
    pulses += "ACGTAC";
    const PulseBitset bs{pulses};
    EXPECT_EQ(bs.Size(), 70U);
    EXPECT_EQ(bs.Count(), 6U);
    EXPECT_EQ(bs.FindFirst(), 64U);
    EXPECT_EQ(bs.FindNext(64), 65U);
}

TEST(PulseBitset, FindNthBase)
{
    const PulseBitset bs{"aACgTtGc"};  // basecalled at: 1(A), 2(C), 4(T), 6(G)
    EXPECT_EQ(bs.FindNthBase(0), 1U);
    EXPECT_EQ(bs.FindNthBase(1), 2U);
    EXPECT_EQ(bs.FindNthBase(2), 4U);
    EXPECT_EQ(bs.FindNthBase(3), 6U);
    EXPECT_EQ(bs.FindNthBase(4), PulseBitset::NPOS);
}

TEST(PulseBitset, Empty)
{
    const PulseBitset bs{""};
    EXPECT_EQ(bs.Size(), 0U);
    EXPECT_EQ(bs.Count(), 0U);
    EXPECT_EQ(bs.FindFirst(), PulseBitset::NPOS);
}

}  // namespace Samoa
}  // namespace PacBio
