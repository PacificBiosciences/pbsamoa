#include <pbsamoa/core/Bgzf.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <format>
#include <mutex>
#include <optional>
#include <span>
#include <thread>
#include <vector>

namespace PacBio {
namespace Samoa {

class BgzfWriterTest : public ::testing::Test
{
protected:
    std::filesystem::path tmpPath_{};

    void SetUp() override
    {
        const std::size_t threadHash{std::hash<std::thread::id>{}(std::this_thread::get_id())};
        tmpPath_ = std::filesystem::temp_directory_path() /
                   std::format("pbsamoa_test_{}.bgzf", threadHash);
    }

    void TearDown() override { std::filesystem::remove(tmpPath_); }
};

TEST(BgzfWriterConfigTest, DefaultValues)
{
    const BgzfWriterConfig config{};
    EXPECT_EQ(config.CompressionLevel, 6);
    EXPECT_EQ(config.BgzfWorkers, 4U);
    EXPECT_EQ(config.InputQueueCapacity, 256U);
    EXPECT_EQ(config.BlocksPerBatch, 32);
}

TEST(BgzfWriteMetricsTest, DefaultValues)
{
    const BgzfWriteMetrics metrics{};
    EXPECT_EQ(metrics.CallerStalls, 0U);
    EXPECT_EQ(metrics.PackerStalls, 0U);
    EXPECT_EQ(metrics.CompressNs, 0U);
    EXPECT_EQ(metrics.BytesCompressed, 0U);
    EXPECT_EQ(metrics.BlocksWritten, 0U);
    EXPECT_EQ(metrics.IoWriteNs, 0U);
    EXPECT_EQ(metrics.WriterStalls, 0U);
    EXPECT_EQ(metrics.CallbackNs, 0U);
}

TEST_F(BgzfWriterTest, WriteAndReadBack)
{
    const std::string original{"Hello, BGZF round-trip!"};

    {
        BgzfWriter writer{tmpPath_};
        writer.Write(std::as_bytes(std::span{std::data(original), std::size(original)}));
    }

    BgzfReader reader{tmpPath_};
    EXPECT_TRUE(reader.HasEofMarker());

    std::vector<std::byte> buffer{};
    buffer.resize(65536U);
    const std::optional<std::size_t> bytesRead{reader.ReadBlock(buffer)};
    ASSERT_TRUE(bytesRead.has_value());
    ASSERT_EQ(*bytesRead, std::size(original));

    const std::string result{reinterpret_cast<const char*>(std::data(buffer)), *bytesRead};
    EXPECT_EQ(result, original);
}

TEST_F(BgzfWriterTest, LargeDataSpansMultipleBlocks)
{
    std::vector<std::byte> input{};
    input.resize(200000U);
    for (std::size_t i{0}; i < std::size(input); ++i) {
        input[i] = static_cast<std::byte>(i & 0xFFU);
    }

    {
        BgzfWriter writer{tmpPath_};
        writer.Write(input);
    }

    BgzfReader reader{tmpPath_};
    std::vector<std::byte> output{};
    std::vector<std::byte> buffer{};
    buffer.resize(65536U);
    while (true) {
        const std::optional<std::size_t> bytesRead{reader.ReadBlock(buffer)};
        if ((!bytesRead.has_value()) || (*bytesRead == 0U)) {
            break;
        }
        output.insert(std::ranges::end(output), std::ranges::begin(buffer),
                      std::ranges::begin(buffer) + *bytesRead);
    }

    ASSERT_EQ(std::size(output), std::size(input));
    EXPECT_TRUE(std::ranges::equal(output, input));
}

TEST_F(BgzfWriterTest, EmptyFileHasEofMarker)
{
    {
        const BgzfWriter writer{tmpPath_};
    }

    BgzfReader reader{tmpPath_};
    EXPECT_TRUE(reader.HasEofMarker());

    std::vector<std::byte> buffer{};
    buffer.resize(65536U);
    const std::optional<std::size_t> bytesRead{reader.ReadBlock(buffer)};
    EXPECT_TRUE((!bytesRead.has_value()) || (*bytesRead == 0U));
}

TEST_F(BgzfWriterTest, MultipleSmallWrites)
{
    const std::string part1{"Hello, "};
    const std::string part2{"world!"};

    {
        BgzfWriter writer{tmpPath_};
        writer.Write(std::as_bytes(std::span{std::data(part1), std::size(part1)}));
        writer.Write(std::as_bytes(std::span{std::data(part2), std::size(part2)}));
    }

    BgzfReader reader{tmpPath_};
    std::vector<std::byte> buffer{};
    buffer.resize(65536U);
    const std::optional<std::size_t> bytesRead{reader.ReadBlock(buffer)};
    ASSERT_TRUE(bytesRead.has_value());

    const std::string result{reinterpret_cast<const char*>(std::data(buffer)), *bytesRead};
    EXPECT_EQ(result, "Hello, world!");
}

TEST_F(BgzfWriterTest, ConstructWithConfig)
{
    const BgzfWriterConfig config{
        .CompressionLevel = 4,
        .BgzfWorkers = 2,
    };
    BgzfWriter writer{tmpPath_, config};
    writer.Close();
    EXPECT_TRUE(std::filesystem::exists(tmpPath_));
}

TEST_F(BgzfWriterTest, WriteWithCallback)
{
    std::vector<std::pair<std::int64_t, std::vector<std::byte>>> received;
    std::mutex mu;

    auto callback = [&](std::int64_t offset, std::span<const std::byte> data) {
        const std::lock_guard lock{mu};
        received.emplace_back(offset, std::vector<std::byte>{std::begin(data), std::end(data)});
    };

    BgzfWriter writer{tmpPath_, BgzfWriterConfig{.BgzfWorkers = 2}};
    writer.SetCallback(std::move(callback));

    for (int i{0}; i < 3; ++i) {
        std::vector<std::byte> record(100);
        std::ranges::fill(record, static_cast<std::byte>(i));
        const PendingCallback cb{
            .rawData = record,
            .active = true,
        };
        writer.Write(record, cb);
    }
    writer.Close();

    const std::lock_guard lock{mu};
    ASSERT_EQ(std::size(received), 3U);
    for (std::size_t i{1}; i < std::size(received); ++i) {
        EXPECT_GE(received[i].first, received[i - 1].first);
    }
}

TEST_F(BgzfWriterTest, MetricsPopulatedAfterWrite)
{
    BgzfWriter writer{tmpPath_, BgzfWriterConfig{.BgzfWorkers = 2}};
    std::vector<std::byte> chunk(32000U);
    for (int i{0}; i < 20; ++i) {
        writer.Write(chunk);
    }
    writer.Close();

    const BgzfWriteMetrics metrics{writer.GetMetrics()};
    EXPECT_GT(metrics.BlocksWritten, 0U);
    EXPECT_GT(metrics.BytesCompressed, 0U);
    EXPECT_GT(metrics.CompressNs, 0U);
    EXPECT_GT(metrics.IoWriteNs, 0U);
}

TEST_F(BgzfWriterTest, OversizedRecordIsSplitAcrossBlocksAndCallbackFiresOnce)
{
    std::vector<std::byte> input(200000U);
    for (std::size_t i{0}; i < std::size(input); ++i) {
        input[i] = static_cast<std::byte>(i & 0xFFU);
    }

    std::vector<std::vector<std::byte>> callbackData;
    std::vector<std::int64_t> callbackOffsets;
    std::mutex mu;

    BgzfWriter writer{tmpPath_, BgzfWriterConfig{
                                    .BgzfWorkers = 2,
                                    .BlocksPerBatch = 1,
                                }};
    writer.SetCallback([&](std::int64_t offset, std::span<const std::byte> rawData) {
        const std::lock_guard lock{mu};
        callbackOffsets.push_back(offset);
        callbackData.emplace_back(std::ranges::begin(rawData), std::ranges::end(rawData));
    });

    const PendingCallback cb{
        .rawData = input,
        .active = true,
    };
    writer.Write(input, cb);
    writer.Close();

    {
        const std::lock_guard lock{mu};
        ASSERT_EQ(std::size(callbackOffsets), 1U);
        ASSERT_EQ(std::size(callbackData), 1U);
        EXPECT_TRUE(std::ranges::equal(callbackData.front(), input));
    }

    BgzfReader reader{tmpPath_};
    std::vector<std::byte> output{};
    std::vector<std::byte> buffer(65536U);
    while (true) {
        const std::optional<std::size_t> bytesRead{reader.ReadBlock(buffer)};
        ASSERT_TRUE(bytesRead.has_value());
        if (*bytesRead == 0U) {
            break;
        }
        output.insert(std::ranges::end(output), std::ranges::begin(buffer),
                      std::ranges::begin(buffer) + static_cast<std::ptrdiff_t>(*bytesRead));
    }
    EXPECT_TRUE(std::ranges::equal(output, input));

    const BgzfWriteMetrics metrics{writer.GetMetrics()};
    EXPECT_GE(metrics.BlocksWritten, 2U);
}

}  // namespace Samoa
}  // namespace PacBio
