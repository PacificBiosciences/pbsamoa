#include <pbsamoa/index/ZmwIndex.hpp>
#include <pbsamoa/index/ZmwWhitelist.hpp>
#include <pbsamoa/io/ZmiWriter.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <format>
#include <thread>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {

class ZmwWhitelistTest : public ::testing::Test
{
protected:
    std::filesystem::path zmiPath_;

    void SetUp() override
    {
        const std::size_t threadHash{std::hash<std::thread::id>{}(std::this_thread::get_id())};
        zmiPath_ = std::filesystem::temp_directory_path() /
                   std::format("pbsamoa_whitelist_test_{}.zmi", threadHash);

        // Write a ZMI with known data:
        //   rgId=100, zmw=10, offset=1000
        //   rgId=100, zmw=10, offset=1100  (second subread)
        //   rgId=100, zmw=20, offset=2000
        //   rgId=100, zmw=30, offset=3000
        //   rgId=200, zmw=10, offset=4000  (same zmw, different RG)
        //   rgId=200, zmw=40, offset=5000
        ZmiWriter writer{zmiPath_};
        writer.AddRecord(100, 10, 1000);
        writer.AddRecord(100, 10, 1100);
        writer.AddRecord(100, 20, 2000);
        writer.AddRecord(100, 30, 3000);
        writer.AddRecord(200, 10, 4000);
        writer.AddRecord(200, 40, 5000);
    }

    void TearDown() override { std::filesystem::remove(zmiPath_); }
};

TEST_F(ZmwWhitelistTest, ResolveByHoleNumber)
{
    const ZmwIndex index{ZmwIndex::FromZmi(zmiPath_)};
    const ZmwWhitelist whitelist{std::vector<std::int32_t>{10, 30}};

    const auto offsets{whitelist.Resolve(index)};

    // zmw=10: offsets 1000, 1100, 4000 (all RGs)
    // zmw=30: offset 3000
    // sorted in file order: 1000, 1100, 3000, 4000
    ASSERT_EQ(std::size(offsets), 4u);
    EXPECT_EQ(offsets[0], 1000);
    EXPECT_EQ(offsets[1], 1100);
    EXPECT_EQ(offsets[2], 3000);
    EXPECT_EQ(offsets[3], 4000);
}

TEST_F(ZmwWhitelistTest, ResolveByIdentity)
{
    const ZmwIndex index{ZmwIndex::FromZmi(zmiPath_)};
    const ZmwWhitelist whitelist{
        std::vector<ZmwIdentity>{ZmwIdentity{100, 10}, ZmwIdentity{200, 40}}};

    const auto offsets{whitelist.Resolve(index)};

    // (100,10): 1000, 1100
    // (200,40): 5000
    // sorted: 1000, 1100, 5000
    ASSERT_EQ(std::size(offsets), 3u);
    EXPECT_EQ(offsets[0], 1000);
    EXPECT_EQ(offsets[1], 1100);
    EXPECT_EQ(offsets[2], 5000);
}

TEST_F(ZmwWhitelistTest, ResolveEmptyWhitelist)
{
    const ZmwIndex index{ZmwIndex::FromZmi(zmiPath_)};
    const ZmwWhitelist whitelist{std::vector<std::int32_t>{}};

    const auto offsets{whitelist.Resolve(index)};
    EXPECT_TRUE(std::empty(offsets));
}

TEST_F(ZmwWhitelistTest, ResolveNonexistentZmw)
{
    const ZmwIndex index{ZmwIndex::FromZmi(zmiPath_)};
    const ZmwWhitelist whitelist{std::vector<std::int32_t>{9999}};

    const auto offsets{whitelist.Resolve(index)};
    EXPECT_TRUE(std::empty(offsets));
}

TEST_F(ZmwWhitelistTest, ResolveDeduplicated)
{
    const ZmwIndex index{ZmwIndex::FromZmi(zmiPath_)};
    // Duplicate zmw in input — offsets should not be duplicated
    const ZmwWhitelist whitelist{std::vector<std::int32_t>{10, 10}};

    const auto offsets{whitelist.Resolve(index)};

    // zmw=10 offsets (all RGs): 1000, 1100, 4000 — each appears once
    ASSERT_EQ(std::size(offsets), 3u);
    EXPECT_TRUE(std::ranges::is_sorted(offsets));
}

}  // namespace Samoa
}  // namespace PacBio
