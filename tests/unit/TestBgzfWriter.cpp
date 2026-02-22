#include <pbsamoa/core/Bgzf.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <format>
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

TEST_F(BgzfWriterTest, TellReturnsVirtualOffset)
{
    BgzfWriter writer{tmpPath_};

    // Before any writes, Tell should be (0, 0)
    const VirtualOffset pos0{writer.Tell()};
    EXPECT_EQ(pos0.BlockOffset(), 0u);
    EXPECT_EQ(pos0.WithinBlockOffset(), 0u);

    // Write some data (smaller than one block)
    const std::string data{"Hello"};
    writer.Write(std::as_bytes(std::span{std::data(data), std::size(data)}));

    // Within-block offset should advance, block offset still 0
    const VirtualOffset pos1{writer.Tell()};
    EXPECT_EQ(pos1.BlockOffset(), 0u);
    EXPECT_EQ(pos1.WithinBlockOffset(), 5u);

    // Write enough to force a block flush (MAX_UNCOMPRESSED_SIZE = 0xFF00 = 65280)
    std::vector<std::byte> filler(65280 - 5);
    writer.Write(filler);

    // After flush: block offset should be nonzero, within-block offset 0
    const VirtualOffset pos2{writer.Tell()};
    EXPECT_GT(pos2.BlockOffset(), 0u);
    EXPECT_EQ(pos2.WithinBlockOffset(), 0u);
}

}  // namespace Samoa
}  // namespace PacBio
