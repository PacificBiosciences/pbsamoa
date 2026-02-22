#include "TestData.hpp"

#include <pbsamoa/core/Bgzf.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <optional>
#include <span>
#include <vector>
#include <cstddef>

namespace PacBio {
namespace Samoa {

TEST(BgzfReader, OpenValidBam)
{
    const std::filesystem::path path{tests::DataDir / "header_only.bam"};
    const BgzfReader reader{path};
    EXPECT_TRUE(reader.HasEofMarker());
}

TEST(BgzfReader, ReadFirstBlock)
{
    const std::filesystem::path path{tests::DataDir / "header_only.bam"};
    BgzfReader reader{path};

    std::vector<std::byte> buffer{};
    buffer.resize(65536U);
    const std::optional<std::size_t> bytesRead{reader.ReadBlock(buffer)};
    ASSERT_TRUE(bytesRead.has_value());
    EXPECT_GT(*bytesRead, 0U);
}

TEST(BgzfReader, ReadAllBlocks)
{
    const std::filesystem::path path{tests::DataDir / "spec_example.bam"};
    BgzfReader reader{path};

    std::size_t totalDecompressed{0};
    std::size_t blockCount{0};
    std::vector<std::byte> buffer{};
    buffer.resize(65536U);

    while (true) {
        const std::optional<std::size_t> bytesRead{reader.ReadBlock(buffer)};
        if ((!bytesRead.has_value()) || (*bytesRead == 0U)) {
            break;
        }
        totalDecompressed += *bytesRead;
        ++blockCount;
    }

    EXPECT_GT(blockCount, 0U);
    EXPECT_GT(totalDecompressed, 0U);
}

TEST(BgzfReader, EofMarkerDetected)
{
    const std::filesystem::path path{tests::DataDir / "spec_example.bam"};
    const BgzfReader reader{path};
    EXPECT_TRUE(reader.HasEofMarker());
}

TEST(BgzfReader, SeekToVirtualOffset)
{
    const std::filesystem::path path{tests::DataDir / "spec_example.bam"};
    BgzfReader reader{path};

    std::vector<std::byte> firstRead{};
    firstRead.resize(65536U);
    const std::optional<std::size_t> firstReadSize{reader.ReadBlock(firstRead)};
    ASSERT_TRUE(firstReadSize.has_value());

    reader.Seek(VirtualOffset{0U, 0U});

    std::vector<std::byte> secondRead{};
    secondRead.resize(65536U);
    const std::optional<std::size_t> secondReadSize{reader.ReadBlock(secondRead)};
    ASSERT_TRUE(secondReadSize.has_value());
    ASSERT_EQ(*firstReadSize, *secondReadSize);

    EXPECT_TRUE(std::ranges::equal(std::span<const std::byte>{firstRead}.first(*firstReadSize),
                                   std::span<const std::byte>{secondRead}.first(*secondReadSize)));
}

TEST(BgzfReader, ThrowOnNonexistentFile)
{
    EXPECT_THROW(BgzfReader{"/nonexistent/file.bam"}, std::runtime_error);
}

}  // namespace Samoa
}  // namespace PacBio
