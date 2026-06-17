#ifndef PBSAMOA_TOOLS_CLIUTILS_HPP
#define PBSAMOA_TOOLS_CLIUTILS_HPP

#include "ParseUtils.hpp"

#include <pbsamoa/io/BamSort.hpp>

#include <array>
#include <format>
#include <stdexcept>
#include <string>
#include <string_view>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {
namespace Tools {

inline bool HasFollowingArgument(int index, int argc) { return (index + 1) < argc; }

inline std::string_view NextArgumentView(char* const* argv, int& index)
{
    ++index;
    return argv[index];
}

inline const char* RequireOptionValue(int argc, char** argv, int& currentIndex,
                                      std::string_view option)
{
    if (!HasFollowingArgument(currentIndex, argc)) {
        throw std::runtime_error{std::format("missing value for {}", option)};
    }
    ++currentIndex;
    return argv[currentIndex];
}

template <typename IntType>
IntType ParseNextIntegerArgument(char* const* argv, int& index, std::string_view name)
{
    return ParseIntegerOrThrow<IntType>(NextArgumentView(argv, index), name);
}

template <typename IntType>
IntType ParseIntOption(int argc, char** argv, int& currentIndex, std::string_view option,
                       std::string_view parseName)
{
    return ParseIntegerOrThrow<IntType>(RequireOptionValue(argc, argv, currentIndex, option),
                                        parseName);
}

/// Build a "pbsamoa <tool> arg1 arg2 ..." command-line string for @PG CL tags.
inline std::string BuildCommandLine(std::string_view toolName, int argc, char** argv)
{
    std::string commandLine{"pbsamoa "};
    commandLine += toolName;
    for (int i{0}; i < argc; ++i) {
        commandLine += ' ';
        commandLine += argv[i];
    }
    return commandLine;
}

/// Parse --order coordinate|queryname|tag.  Throws on invalid value.
inline SortOrder ParseSortOrder(int argc, char** argv, int& i)
{
    const std::string_view value{RequireOptionValue(argc, argv, i, "--order")};
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

/// Parse --tag XX (exactly 2 chars).  Throws on invalid value.
inline std::array<char, 2> ParseSortTag(int argc, char** argv, int& i)
{
    const std::string_view value{RequireOptionValue(argc, argv, i, "--tag")};
    if (std::size(value) != 2) {
        throw std::runtime_error{std::format("--tag must be exactly 2 characters: {}", value)};
    }
    return {value[0], value[1]};
}

/// Parse --threads N (>= 0); 0 means auto.  Returns std::size_t.
inline std::size_t ParseThreadsOption(int argc, char** argv, int& i)
{
    const std::int32_t threads{ParseIntegerOrThrow<std::int32_t>(
        RequireOptionValue(argc, argv, i, "--threads"), "--threads")};
    if (threads < 0) {
        throw std::runtime_error{"--threads must be >= 0"};
    }
    return static_cast<std::size_t>(threads);
}

/// Parse --compression L ([1, 12]).  Returns the level as int.
inline std::int32_t ParseCompressionLevelOption(int argc, char** argv, int& i)
{
    const std::int32_t level{ParseIntegerOrThrow<std::int32_t>(
        RequireOptionValue(argc, argv, i, "--compression"), "--compression")};
    if ((level < 1) || (level > 12)) {
        throw std::runtime_error{"--compression must be in [1, 12]"};
    }
    return level;
}

}  // namespace Tools
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_CLIUTILS_HPP
