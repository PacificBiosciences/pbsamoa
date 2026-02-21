#include <pbsamoa/core/Bgzf.hpp>

#include <gtest/gtest.h>

namespace PacBio {
namespace Samoa {

TEST(VirtualOffset, DefaultConstructsToZero)
{
    constexpr VirtualOffset vo;
    EXPECT_EQ(vo.BlockOffset(), 0u);
    EXPECT_EQ(vo.WithinBlockOffset(), 0u);
    EXPECT_EQ(vo.Value(), 0u);
}

TEST(VirtualOffset, ConstructFromRawValue)
{
    // Block offset 100, within-block offset 50
    // raw = (100 << 16) | 50 = 6553650
    constexpr VirtualOffset vo{6553650u};
    EXPECT_EQ(vo.BlockOffset(), 100u);
    EXPECT_EQ(vo.WithinBlockOffset(), 50u);
    EXPECT_EQ(vo.Value(), 6553650u);
}

TEST(VirtualOffset, ConstructFromComponents)
{
    constexpr VirtualOffset vo{100u, 50u};
    EXPECT_EQ(vo.BlockOffset(), 100u);
    EXPECT_EQ(vo.WithinBlockOffset(), 50u);
}

TEST(VirtualOffset, MaxWithinBlockOffset)
{
    // Within-block offset is 16 bits, max = 65535
    constexpr VirtualOffset vo{0u, 65535u};
    EXPECT_EQ(vo.WithinBlockOffset(), 65535u);
    EXPECT_EQ(vo.BlockOffset(), 0u);
}

TEST(VirtualOffset, LargeBlockOffset)
{
    // Block offset uses upper 48 bits
    const std::uint64_t largeOffset{0x0000FFFFFFFFu};
    const VirtualOffset vo{largeOffset, 0u};
    EXPECT_EQ(vo.BlockOffset(), largeOffset);
    EXPECT_EQ(vo.WithinBlockOffset(), 0u);
}

TEST(VirtualOffset, ComparisonOrdering)
{
    constexpr VirtualOffset a{100u, 0u};
    constexpr VirtualOffset b{100u, 50u};
    constexpr VirtualOffset c{200u, 0u};

    EXPECT_LT(a, b);
    EXPECT_LT(b, c);
    EXPECT_LT(a, c);
    EXPECT_EQ(a, a);
    EXPECT_NE(a, b);
}

TEST(VirtualOffset, Equality)
{
    constexpr VirtualOffset a{42u, 7u};
    constexpr VirtualOffset b{42u, 7u};
    EXPECT_EQ(a, b);
}

TEST(VirtualOffset, RoundTripComponentsToValue)
{
    for (const std::uint64_t block : {0u, 1u, 1000u, 0xFFFFFFu}) {
        for (const std::uint16_t within : {0, 1, 100, 65535}) {
            const VirtualOffset vo{block, within};
            const VirtualOffset vo2{vo.Value()};
            EXPECT_EQ(vo, vo2);
            EXPECT_EQ(vo.BlockOffset(), vo2.BlockOffset());
            EXPECT_EQ(vo.WithinBlockOffset(), vo2.WithinBlockOffset());
        }
    }
}

}  // namespace Samoa
}  // namespace PacBio
