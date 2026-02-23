#ifndef PBSAMOA_TOOLS_PARSEUTILS_HPP
#define PBSAMOA_TOOLS_PARSEUTILS_HPP

#include <charconv>
#include <expected>
#include <format>
#include <string>
#include <string_view>
#include <type_traits>

namespace PacBio {
namespace Samoa {
namespace Tools {

template <typename IntType>
std::expected<IntType, std::string> ParseInteger(std::string_view text, std::string_view name)
{
    static_assert(std::is_integral_v<IntType>, "ParseInteger requires an integral target type");

    IntType value{};
    const char* const begin{text.data()};
    const char* const end{text.data() + text.size()};
    const auto [ptr, ec]{std::from_chars(begin, end, value)};

    if ((ec != std::errc{}) || (ptr != end)) {
        return std::unexpected{std::format("invalid {}: {}", name, text)};
    }
    return value;
}

}  // namespace Tools
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_PARSEUTILS_HPP
