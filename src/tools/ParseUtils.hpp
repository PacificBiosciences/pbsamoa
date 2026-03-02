#ifndef PBSAMOA_TOOLS_PARSEUTILS_HPP
#define PBSAMOA_TOOLS_PARSEUTILS_HPP

#include <algorithm>
#include <charconv>
#include <expected>
#include <format>
#include <string>
#include <string_view>
#include <thread>
#include <type_traits>

#include <cstddef>
#include <cstdint>

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

/// \brief Resolve a user-requested worker count to an actual thread count.
///
/// \param requested  User-requested count (< 0 means auto-detect)
/// \param explicitCap  Maximum when user provides explicit count
/// \param autoCap  Maximum when auto-detecting (default: 8)
inline std::size_t ResolveNumWorkers(std::int32_t requested, std::int32_t explicitCap,
                                     std::int32_t autoCap = 8)
{
    const std::int32_t hwThreads = std::ranges::max(std::thread::hardware_concurrency(), 1U);

    if (requested < 0) {
        return static_cast<std::size_t>(std::ranges::min(hwThreads, autoCap));
    }
    return static_cast<std::size_t>(std::ranges::min(requested, explicitCap));
}

}  // namespace Tools
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_PARSEUTILS_HPP
