#include "../../src/PulseBitset.hpp"

#include <gtest/gtest.h>

namespace PacBio {
namespace Samoa {

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
}

}  // namespace Samoa
}  // namespace PacBio
