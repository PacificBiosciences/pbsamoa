#include <pbsamoa/io/BamRawReader.hpp>

#include <gtest/gtest.h>

namespace PacBio {
namespace Samoa {

TEST(ByteLimit, ConstructFromValue)
{
    constexpr ByteLimit limit{64 * 1024 * 1024};
    EXPECT_EQ(limit.Value(), 64u * 1024 * 1024);
}

TEST(ByteLimit, DefaultIsMaxSize)
{
    constexpr ByteLimit limit;
    EXPECT_GT(limit.Value(), 0u);
}

TEST(ByteLimit, UserDefinedLiteral)
{
    constexpr ByteLimit limit{64U * 1024U * 1024U};
    EXPECT_EQ(limit.Value(), 64u * 1024 * 1024);
}

TEST(ErrorPolicy, ThrowPolicyThrows)
{
    const ThrowPolicy policy;
    EXPECT_THROW(policy.OnCorruptRecord("test error"), std::runtime_error);
}

TEST(ErrorPolicy, SkipPolicyDoesNotThrow)
{
    SkipPolicy policy;
    EXPECT_NO_THROW(policy.OnCorruptRecord("test error"));
    EXPECT_EQ(policy.SkippedCount(), 1u);
}

}  // namespace Samoa
}  // namespace PacBio
