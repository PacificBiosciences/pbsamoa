#ifndef PBSAMOA_TOOLS_CLIUTILS_HPP
#define PBSAMOA_TOOLS_CLIUTILS_HPP

#include "ParseUtils.hpp"

#include <pbsamoa/io/BamSort.hpp>

#include <array>
#include <format>
#include <limits>
#include <stdexcept>
#include <string_view>

#include <cstdint>

namespace PacBio {
namespace Samoa {
namespace Tools {

// --- Value interpreters for pbcopper CLIv2 option strings ---------------------
// CLIv2 does the scanning/typing; these turn an already-extracted option value
// into a domain type (or validate it).

/// Interpret an --order value.  Throws on an unrecognized value.
inline SortOrder ParseSortOrder(std::string_view value)
{
    if (value == "coordinate") {
        return SortOrder::COORDINATE;
    }
    if (value == "queryname") {
        return SortOrder::QUERY_NAME;
    }
    if (value == "tag") {
        return SortOrder::TAG;
    }
    throw std::runtime_error{std::format("invalid --order: {}", value)};
}

/// Interpret a --tag value (exactly 2 characters).  Throws otherwise.
inline std::array<char, 2> ParseSortTag(std::string_view value)
{
    if (std::size(value) != 2) {
        throw std::runtime_error{std::format("--tag must be exactly 2 characters: {}", value)};
    }
    return {value[0], value[1]};
}

/// Interpret a memory SIZE value with an optional K/M/G suffix (base 1024).
/// \p option names the flag for the error message.
inline ByteLimit ParseMemory(std::string_view text, std::string_view option = "--memory")
{
    if (text.empty()) {
        throw std::runtime_error{std::format("invalid {}: empty value", option)};
    }

    std::uint64_t multiplier{1};
    switch (text.back()) {
        case 'k':
        case 'K':
            multiplier = std::uint64_t{1024};
            break;
        case 'm':
        case 'M':
            multiplier = std::uint64_t{1024} * 1024;
            break;
        case 'g':
        case 'G':
            multiplier = std::uint64_t{1024} * 1024 * 1024;
            break;
        default:
            break;
    }
    const std::string_view digits{multiplier != 1 ? text.substr(0, std::size(text) - 1) : text};

    const std::uint64_t value{ParseIntegerOrThrow<std::uint64_t>(digits, option)};
    if (value > (std::numeric_limits<std::size_t>::max() / multiplier)) {
        throw std::runtime_error{std::format("{} is too large: {}", option, text)};
    }
    return ByteLimit{value * multiplier};
}

/// Validate a compression level in [1, 12].  Throws otherwise.
inline std::int32_t CheckCompressionLevel(std::int32_t level)
{
    if ((level < 1) || (level > 12)) {
        throw std::runtime_error{"--compression must be in [1, 12]"};
    }
    return level;
}

}  // namespace Tools
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_CLIUTILS_HPP
