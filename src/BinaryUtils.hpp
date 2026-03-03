#ifndef PBSAMOA_BINARYUTILS_HPP
#define PBSAMOA_BINARYUTILS_HPP

#include <algorithm>
#include <charconv>
#include <filesystem>
#include <format>
#include <fstream>
#include <span>
#include <stdexcept>
#include <string_view>
#include <vector>

#include <bit>
#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief Read a little-endian uint16 from raw bytes.
inline std::uint16_t ReadU16LE(const std::byte* p)
{
    const std::uint16_t low{std::to_integer<std::uint16_t>(p[0])};
    const std::uint16_t high{std::to_integer<std::uint16_t>(p[1])};
    return low | (high << 8U);
}

/// \brief Read a little-endian uint32 from raw bytes.
inline std::uint32_t ReadU32LE(const std::byte* p)
{
    return std::to_integer<std::uint32_t>(p[0]) | (std::to_integer<std::uint32_t>(p[1]) << 8U) |
           (std::to_integer<std::uint32_t>(p[2]) << 16U) |
           (std::to_integer<std::uint32_t>(p[3]) << 24U);
}

/// \brief Read a little-endian int32 from raw bytes.
inline std::int32_t ReadI32LE(const std::byte* p)
{
    return std::bit_cast<std::int32_t>(ReadU32LE(p));
}

/// \brief Read a little-endian uint64 from raw bytes.
inline std::uint64_t ReadU64LE(const std::byte* p)
{
    return std::to_integer<std::uint64_t>(p[0]) | (std::to_integer<std::uint64_t>(p[1]) << 8U) |
           (std::to_integer<std::uint64_t>(p[2]) << 16U) |
           (std::to_integer<std::uint64_t>(p[3]) << 24U) |
           (std::to_integer<std::uint64_t>(p[4]) << 32U) |
           (std::to_integer<std::uint64_t>(p[5]) << 40U) |
           (std::to_integer<std::uint64_t>(p[6]) << 48U) |
           (std::to_integer<std::uint64_t>(p[7]) << 56U);
}

/// \brief Read a little-endian int64 from raw bytes.
inline std::int64_t ReadI64LE(const std::byte* p)
{
    return std::bit_cast<std::int64_t>(ReadU64LE(p));
}

/// \brief Write a little-endian int32 to a byte vector.
inline void WriteI32LE(std::vector<std::byte>& out, std::int32_t value)
{
    const auto pos = std::size(out);
    out.resize(pos + sizeof(value));
    std::ranges::copy_n(reinterpret_cast<const std::byte*>(&value), sizeof(value),
                        out.data() + pos);
}

/// \brief Write a little-endian uint32 to a byte vector.
inline void WriteU32LE(std::vector<std::byte>& out, std::uint32_t value)
{
    const auto pos = std::size(out);
    out.resize(pos + sizeof(value));
    std::ranges::copy_n(reinterpret_cast<const std::byte*>(&value), sizeof(value),
                        out.data() + pos);
}

inline constexpr std::size_t MAX_DECOMPRESSED_BLOCK_SIZE{65536U};

/// \brief Compute total BAM header size from binary data.
///
/// Returns the total number of bytes from BAM magic through end of reference
/// dictionary, or 0 if insufficient data.
inline std::size_t ComputeHeaderSize(std::span<const std::byte> data)
{
    // Need at least: magic(4) + l_text(4) = 8
    if (std::size(data) < 8) {
        return 0;
    }

    const std::uint32_t lText{ReadU32LE(std::data(data) + 4)};
    // Need: magic(4) + l_text(4) + text(lText) + n_ref(4)
    const std::size_t minSize{8 + lText + 4};
    if (std::size(data) < minSize) {
        return 0;
    }

    const std::uint32_t nRef{ReadU32LE(std::data(data) + 8 + lText)};
    std::size_t offset{8 + lText + 4};

    for (std::uint32_t i{0}; i < nRef; ++i) {
        if (offset + 4 > std::size(data)) {
            return 0;
        }
        const std::uint32_t lName{ReadU32LE(std::data(data) + offset)};
        offset += 4;
        // name(lName) + l_ref(4)
        if (offset + lName + 4 > std::size(data)) {
            return 0;
        }
        offset += lName + 4;
    }

    return offset;
}

/// \brief Read entire file contents as a byte vector.
inline std::vector<std::byte> ReadAllBytes(const std::filesystem::path& path)
{
    std::ifstream in{path, std::ios::binary};
    if (!in.is_open()) {
        throw std::runtime_error{"cannot open file: " + path.string()};
    }

    in.seekg(0, std::ios::end);
    if (in.fail()) {
        throw std::runtime_error{"failed to read file: " + path.string()};
    }
    const std::streamoff fileSize = in.tellg();
    if (fileSize < 0) {
        throw std::runtime_error{"failed to read file: " + path.string()};
    }
    in.seekg(0, std::ios::beg);
    if (in.fail()) {
        throw std::runtime_error{"failed to read file: " + path.string()};
    }

    std::vector<std::byte> data(static_cast<std::size_t>(fileSize));
    if (fileSize > 0) {
        in.read(reinterpret_cast<char*>(data.data()), fileSize);
        if (in.gcount() != fileSize) {
            throw std::runtime_error{"failed to read file: " + path.string()};
        }
    }
    return data;
}

/// \brief Parse a string_view as an integer using std::from_chars (throws on failure).
///
/// \param[in] value      The string to parse.
/// \param[in] context    Human-readable context for error messages.
template <typename Integer>
Integer ParseInteger(std::string_view value, std::string_view context)
{
    Integer parsed{};
    const auto [ptr, ec] = std::from_chars(value.data(), value.data() + std::size(value), parsed);
    if (ec != std::errc{} || ptr != (value.data() + std::size(value))) {
        throw std::runtime_error{std::format("{}: not numeric: '{}'", context, value)};
    }
    return parsed;
}

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_BINARYUTILS_HPP
