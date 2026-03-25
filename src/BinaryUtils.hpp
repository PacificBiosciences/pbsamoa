#ifndef PBSAMOA_BINARYUTILS_HPP
#define PBSAMOA_BINARYUTILS_HPP

#include <pbsamoa/core/Endian.hpp>

#include <algorithm>
#include <charconv>
#include <concepts>
#include <filesystem>
#include <format>
#include <fstream>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief Write a little-endian integral value to raw bytes (unaligned-safe).
template <std::integral T>
inline void WriteLE(std::byte* dst, T value)
{
    std::ranges::copy_n(reinterpret_cast<const std::byte*>(&value), sizeof(value), dst);
}

template <std::integral T>
inline void AppendLE(std::vector<std::byte>& out, T value)
{
    const std::size_t appendOffset{std::size(out)};
    out.resize(appendOffset + sizeof(value));
    WriteLE(std::data(out) + appendOffset, value);
}

/// \brief Write a little-endian int32 to a byte vector.
inline void WriteI32LE(std::vector<std::byte>& out, std::int32_t value) { AppendLE(out, value); }

/// \brief Write a little-endian uint32 to a byte vector.
inline void WriteU32LE(std::vector<std::byte>& out, std::uint32_t value) { AppendLE(out, value); }

inline constexpr std::size_t MAX_DECOMPRESSED_BLOCK_SIZE{65536U};

[[noreturn]] inline void ThrowReadFileError(const std::string& pathText)
{
    throw std::runtime_error{"failed to read file: " + pathText};
}

/// \brief Compute total BAM header size from binary data.
///
/// Returns the total number of bytes from BAM magic through end of reference
/// dictionary, or 0 if insufficient data.
inline std::size_t ComputeHeaderSize(std::span<const std::byte> data)
{
    constexpr std::size_t MAGIC_SIZE{4};
    constexpr std::size_t LTEXT_SIZE{4};
    constexpr std::size_t NREF_SIZE{4};
    constexpr std::size_t REF_NAME_LENGTH_SIZE{4};
    constexpr std::size_t REF_LENGTH_SIZE{4};
    constexpr std::size_t PREFIX_SIZE{MAGIC_SIZE + LTEXT_SIZE};
    const std::byte* const dataBegin{std::data(data)};
    const std::size_t dataSize{std::size(data)};

    // Need at least: magic(4) + l_text(4)
    if (dataSize < PREFIX_SIZE) {
        return 0;
    }

    const std::uint32_t lText{ReadU32LE(dataBegin + MAGIC_SIZE)};
    // Need: magic(4) + l_text(4) + text(lText) + n_ref(4)
    const std::size_t minSize{PREFIX_SIZE + lText + NREF_SIZE};
    if (dataSize < minSize) {
        return 0;
    }

    const std::uint32_t nRef{ReadU32LE(dataBegin + PREFIX_SIZE + lText)};
    std::size_t offset{PREFIX_SIZE + lText + NREF_SIZE};

    for (std::uint32_t i{0}; i < nRef; ++i) {
        if (offset + REF_NAME_LENGTH_SIZE > dataSize) {
            return 0;
        }
        const std::uint32_t lName{ReadU32LE(dataBegin + offset)};
        offset += REF_NAME_LENGTH_SIZE;
        // name(lName) + l_ref(4)
        if (offset + lName + REF_LENGTH_SIZE > dataSize) {
            return 0;
        }
        offset += lName + REF_LENGTH_SIZE;
    }

    return offset;
}

/// \brief Read entire file contents as a byte vector.
inline std::vector<std::byte> ReadAllBytes(const std::filesystem::path& path)
{
    const std::string pathText{path.string()};
    std::ifstream in{path, std::ios::binary};
    if (!in.is_open()) {
        throw std::runtime_error{"cannot open file: " + pathText};
    }

    in.seekg(0, std::ios::end);
    if (in.fail()) {
        ThrowReadFileError(pathText);
    }
    const std::streamoff fileSize{in.tellg()};
    if (fileSize < 0) {
        ThrowReadFileError(pathText);
    }
    in.seekg(0, std::ios::beg);
    if (in.fail()) {
        ThrowReadFileError(pathText);
    }

    const std::streamsize readSize{fileSize};
    std::vector<std::byte> data(static_cast<std::size_t>(readSize));
    if (fileSize > 0) {
        in.read(reinterpret_cast<char*>(data.data()), readSize);
        if (in.gcount() != readSize) {
            ThrowReadFileError(pathText);
        }
    }
    return data;
}

/// \brief Parse a string_view as an integer using std::from_chars (throws on
/// failure).
///
/// \param[in] value      The string to parse.
/// \param[in] context    Human-readable context for error messages.
template <typename Integer>
Integer ParseInteger(std::string_view value, std::string_view context)
{
    Integer parsed{};
    const char* const begin{std::data(value)};
    const char* const end{begin + std::size(value)};
    const std::from_chars_result parseResult{std::from_chars(begin, end, parsed)};
    if ((parseResult.ec != std::errc{}) || (parseResult.ptr != end)) {
        throw std::runtime_error{std::format("{}: not numeric: '{}'", context, value)};
    }
    return parsed;
}

template <std::signed_integral Integer>
inline constexpr Integer OneBasedPositionOrZero(Integer zeroBasedPos)
{
    if (zeroBasedPos < 0) {
        return 0;
    }
    return zeroBasedPos + 1;
}

template <std::signed_integral Integer>
inline constexpr Integer NonEmptyAlignmentEnd(Integer pos, Integer end)
{
    if (end > pos) {
        return end;
    }
    return pos + 1;
}

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_BINARYUTILS_HPP
