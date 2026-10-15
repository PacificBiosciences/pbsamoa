#include "TestData.hpp"
#include "TestTempDir.hpp"

#include <pbsamoa/core/Bgzf.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <optional>
#include <span>
#include <vector>

#include <cstdint>

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
    ASSERT_TRUE(bytesRead);
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
        if (!bytesRead || (*bytesRead == 0U)) {
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
    ASSERT_TRUE(firstReadSize);

    reader.Seek(VirtualOffset{0U, 0U});

    std::vector<std::byte> secondRead{};
    secondRead.resize(65536U);
    const std::optional<std::size_t> secondReadSize{reader.ReadBlock(secondRead)};
    ASSERT_TRUE(secondReadSize);
    ASSERT_EQ(*firstReadSize, *secondReadSize);

    EXPECT_TRUE(std::ranges::equal(std::span<const std::byte>{firstRead}.first(*firstReadSize),
                                   std::span<const std::byte>{secondRead}.first(*secondReadSize)));
}

TEST(BgzfReader, TellInBlockOnlyModeGivesTheBlockReadLast)
{
    // Five blocks: four full ones and one partly filled one.
    std::vector<std::byte> text(4 * BGZF_MAX_UNCOMPRESSED_BLOCK + 100);
    for (std::size_t i = 0; i < std::size(text); ++i) {
        text[i] = std::byte(i % 251);
    }
    BgzfBlockEncoder encoder{6};
    encoder.Write(text);
    encoder.Finish();
    const std::vector<std::byte> blocks{encoder.TakeBlocks()};
    std::vector<std::uint64_t> offsets;
    for (std::size_t at = 0; at < std::size(blocks);) {
        offsets.push_back(at);
        const std::optional<BgzfBlockInfo> info{
            ParseBgzfBlockHeader(std::span<const std::byte>{blocks}.subspan(at))};
        ASSERT_TRUE(info);
        at += info->blockSize;
    }
    ASSERT_EQ(std::size(offsets), 5U);

    const tests::TempDirGuard dir{"bgzf_reader_tell"};
    const std::filesystem::path path{dir.File("blocks.gz")};
    {
        std::ofstream out{path, std::ios::binary};
        out.write(reinterpret_cast<const char*>(std::data(blocks)), std::ssize(blocks));
        out.write(reinterpret_cast<const char*>(std::data(BGZF_EOF_MARKER)),
                  std::ssize(BGZF_EOF_MARKER));
    }

    BgzfReader reader{path};
    EXPECT_EQ(reader.Tell(), (VirtualOffset{0U, 0U}));
    std::vector<std::byte> buffer(65536U);
    for (const std::uint64_t offset : offsets) {
        const std::optional<std::size_t> bytesRead{reader.ReadBlock(buffer)};
        ASSERT_TRUE(bytesRead);
        ASSERT_GT(*bytesRead, 0U);
        EXPECT_EQ(reader.Tell(), (VirtualOffset{offset, 0U}));
    }
    const std::optional<std::size_t> atEnd{reader.ReadBlock(buffer)};
    ASSERT_TRUE(atEnd);
    EXPECT_EQ(*atEnd, 0U);

    reader.Seek(VirtualOffset{offsets[2], 0U});
    EXPECT_EQ(reader.Tell(), (VirtualOffset{offsets[2], 0U}));
    ASSERT_TRUE(reader.ReadBlock(buffer));
    EXPECT_EQ(reader.Tell(), (VirtualOffset{offsets[2], 0U}));
    ASSERT_TRUE(reader.ReadBlock(buffer));
    EXPECT_EQ(reader.Tell(), (VirtualOffset{offsets[3], 0U}));
}

TEST(BgzfReader, ThrowOnNonexistentFile)
{
    EXPECT_THROW(BgzfReader{"/nonexistent/file.bam"}, std::runtime_error);
}

}  // namespace Samoa
}  // namespace PacBio
