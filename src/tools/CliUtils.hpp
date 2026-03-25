#ifndef PBSAMOA_TOOLS_CLIUTILS_HPP
#define PBSAMOA_TOOLS_CLIUTILS_HPP

#include "ParseUtils.hpp"

#include <format>
#include <stdexcept>
#include <string_view>

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

}  // namespace Tools
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_CLIUTILS_HPP
