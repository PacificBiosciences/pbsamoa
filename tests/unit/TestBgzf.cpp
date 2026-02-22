#include <pbsamoa/core/Bgzf.hpp>

#include <libdeflate.h>

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

namespace {

constexpr std::array<std::uint8_t, 28> EOF_BYTES{
    0x1F, 0x8B, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xFF, 0x06, 0x00, 0x42, 0x43,
    0x02, 0x00, 0x1B, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
};

struct LibdeflateCompressorDeleter
{
    void operator()(libdeflate_compressor* compressor) const
    {
        libdeflate_free_compressor(compressor);
    }
};

using CompressorPtr = std::unique_ptr<libdeflate_compressor, LibdeflateCompressorDeleter>;

std::vector<std::uint8_t> MakeBgzfBlock(std::span<const std::byte> input)
{
    const CompressorPtr compressor{libdeflate_alloc_compressor(6)};
    if (!compressor) {
        throw std::runtime_error{"Failed to allocate libdeflate compressor"};
    }

    std::vector<std::uint8_t> cdata{};
    cdata.resize(libdeflate_deflate_compress_bound(compressor.get(), std::size(input)));
    const std::size_t compressedSize{libdeflate_deflate_compress(
        compressor.get(), std::data(input), std::size(input), std::data(cdata), std::size(cdata))};
    cdata.resize(compressedSize);

    const std::uint32_t crc{libdeflate_crc32(0, std::data(input), std::size(input))};
    const std::uint32_t isize = std::size(input);

    const std::uint16_t xlen{6U};
    const std::uint32_t blockSize = 10U + 2U + xlen + compressedSize + 8U;
    const std::uint16_t bsize = blockSize - 1U;

    std::vector<std::uint8_t> block{};
    block.reserve(blockSize);

    block.insert(std::ranges::end(block), {31U, 139U, 8U, 4U});
    block.insert(std::ranges::end(block), {0U, 0U, 0U, 0U});
    block.insert(std::ranges::end(block), {0U, 0U});
    block.push_back(static_cast<std::uint8_t>(xlen & 0xFFU));
    block.push_back(static_cast<std::uint8_t>(xlen >> 8U));
    block.insert(std::ranges::end(block), {66U, 67U});
    block.insert(std::ranges::end(block), {2U, 0U});
    block.push_back(static_cast<std::uint8_t>(bsize & 0xFFU));
    block.push_back(static_cast<std::uint8_t>(bsize >> 8U));
    block.insert(std::ranges::end(block), std::ranges::begin(cdata), std::ranges::end(cdata));

    for (std::int32_t i{0}; i < 4; ++i) {
        block.push_back(static_cast<std::uint8_t>((crc >> (i * 8)) & 0xFFU));
    }
    for (std::int32_t i{0}; i < 4; ++i) {
        block.push_back(static_cast<std::uint8_t>((isize >> (i * 8)) & 0xFFU));
    }

    return block;
}

}  // namespace

TEST(BgzfBlockHeader, ParseEofMarker)
{
    const std::optional<BgzfBlockInfo> result{
        ParseBgzfBlockHeader(std::as_bytes(std::span{EOF_BYTES}))};
    ASSERT_TRUE(result.has_value());

    const BgzfBlockInfo& info{*result};
    EXPECT_EQ(info.blockSize, 28U);
    EXPECT_EQ(info.compressedDataOffset, 18U);
    EXPECT_EQ(info.compressedDataSize, 2U);
}

TEST(BgzfBlockHeader, RejectBadMagic)
{
    std::array<std::uint8_t, 28> bad{EOF_BYTES};
    bad[0] = 0x00;
    const std::optional<BgzfBlockInfo> result{ParseBgzfBlockHeader(std::as_bytes(std::span{bad}))};
    EXPECT_FALSE(result.has_value());
}

TEST(BgzfBlockHeader, RejectMissingExtraFlag)
{
    std::array<std::uint8_t, 28> bad{EOF_BYTES};
    bad[3] = 0x00;
    const std::optional<BgzfBlockInfo> result{ParseBgzfBlockHeader(std::as_bytes(std::span{bad}))};
    EXPECT_FALSE(result.has_value());
}

TEST(BgzfBlockHeader, RejectTruncatedInput)
{
    const std::optional<BgzfBlockInfo> result{
        ParseBgzfBlockHeader(std::as_bytes(std::span{EOF_BYTES}).subspan(0U, 10U))};
    EXPECT_FALSE(result.has_value());
}

TEST(BgzfBlockHeader, EofMarkerConstantMatchesSpec)
{
    ASSERT_EQ(std::size(BGZF_EOF_MARKER), 28U);
    for (std::size_t i{0}; i < 28U; ++i) {
        EXPECT_EQ(static_cast<std::uint8_t>(BGZF_EOF_MARKER[i]), EOF_BYTES[i])
            << "mismatch at byte " << i;
    }
}

TEST(BgzfBlockHeader, CheckEofMarker)
{
    EXPECT_TRUE(IsBgzfEofMarker(std::as_bytes(std::span{EOF_BYTES})));

    std::array<std::uint8_t, 28> bad{EOF_BYTES};
    bad[27] = 0xFF;
    EXPECT_FALSE(IsBgzfEofMarker(std::as_bytes(std::span{bad})));
}

TEST(BgzfDecompress, RoundTripSmallData)
{
    const std::string original{"Hello, BGZF world!"};
    const std::span<const std::byte> inputSpan{
        std::as_bytes(std::span{std::data(original), std::size(original)})};
    const std::vector<std::uint8_t> block{MakeBgzfBlock(inputSpan)};

    const std::span<const std::byte> blockSpan{std::as_bytes(std::span{block})};
    const std::optional<BgzfBlockInfo> info{ParseBgzfBlockHeader(blockSpan)};
    ASSERT_TRUE(info.has_value());

    std::vector<std::byte> output{};
    output.resize(std::size(original));
    const std::optional<std::size_t> result{DecompressBgzfBlock(blockSpan, *info, output)};
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(*result, std::size(original));

    const std::string decoded{reinterpret_cast<const char*>(std::data(output)), std::size(output)};
    EXPECT_EQ(decoded, original);
}

TEST(BgzfDecompress, EmptyBlock)
{
    const std::span<const std::byte> blockSpan{std::as_bytes(std::span{EOF_BYTES})};
    const std::optional<BgzfBlockInfo> info{ParseBgzfBlockHeader(blockSpan)};
    ASSERT_TRUE(info.has_value());

    std::vector<std::byte> output{};
    output.resize(65536U);
    const std::optional<std::size_t> result{DecompressBgzfBlock(blockSpan, *info, output)};
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(*result, 0U);
}

TEST(BgzfDecompress, MaxBlockData)
{
    std::vector<std::byte> input{};
    input.resize(65536U);
    for (std::size_t i{0}; i < std::size(input); ++i) {
        input[i] = static_cast<std::byte>(i & 0xFFU);
    }
    const std::vector<std::uint8_t> block{MakeBgzfBlock(input)};

    const std::span<const std::byte> blockSpan{std::as_bytes(std::span{block})};
    const std::optional<BgzfBlockInfo> info{ParseBgzfBlockHeader(blockSpan)};
    ASSERT_TRUE(info.has_value());

    std::vector<std::byte> output{};
    output.resize(65536U);
    const std::optional<std::size_t> result{DecompressBgzfBlock(blockSpan, *info, output)};
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(*result, 65536U);
    EXPECT_TRUE(std::ranges::equal(output, input));
}

}  // namespace Samoa
}  // namespace PacBio
