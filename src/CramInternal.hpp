#ifndef PBSAMOA_CRAM_INTERNAL_HPP
#define PBSAMOA_CRAM_INTERNAL_HPP

#include <pbsamoa/core/Tags.hpp>

#include "LibdeflateUtils.hpp"

#include <compare>
#include <filesystem>
#include <span>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

struct TagTriple
{
    char Tag1{};
    char Tag2{};
    char Type{};

    auto operator<=>(const TagTriple&) const = default;
};

inline constexpr std::int32_t TagContentId(char tag1, char tag2, char type)
{
    return (static_cast<std::int32_t>(static_cast<std::uint8_t>(tag1)) << 16) |
           (static_cast<std::int32_t>(static_cast<std::uint8_t>(tag2)) << 8) |
           static_cast<std::int32_t>(static_cast<std::uint8_t>(type));
}

inline constexpr std::int32_t TagContentId(TagKey key, char type)
{
    return TagContentId(key.First(), key.Second(), type);
}

struct EncodedTagValuePayload
{
    char Type{};
    std::vector<std::byte> Payload;
};

EncodedTagValuePayload EncodeTagValueToBamPayload(const TagValue& value);
TagValue DecodeTagValueFromBamPayload(char type, std::span<const std::byte> payload);

struct CramGzipDecompressorContext
{
    LibdeflateDecompressorPtr Decompressor{libdeflate_alloc_decompressor()};
};

inline constexpr int DEFAULT_GZIP_COMPRESSION_LEVEL = 6;

struct CramGzipCompressorContext
{
    LibdeflateCompressorPtr Compressor{libdeflate_alloc_compressor(DEFAULT_GZIP_COMPRESSION_LEVEL)};
};

// CRAM record flags (CF data series)
inline constexpr std::int32_t CRAM_FLAG_QUALITY_AS_ARRAY = 0x1;
inline constexpr std::int32_t CRAM_FLAG_DETACHED = 0x2;
inline constexpr std::int32_t CRAM_FLAG_HAS_MATE_DOWNSTREAM = 0x4;
inline constexpr std::int32_t CRAM_FLAG_SEQUENCE_OMITTED = 0x8;

inline constexpr TagKey RG_TAG{'R', 'G'};

inline std::filesystem::path DefaultCraiPath(const std::filesystem::path& cramPath)
{
    auto p = cramPath;
    p += ".crai";
    return p;
}

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_CRAM_INTERNAL_HPP
