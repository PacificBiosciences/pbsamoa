#ifndef PBSAMOA_SAMFIELDUTILS_HPP
#define PBSAMOA_SAMFIELDUTILS_HPP

#include "WriterUtils.hpp"

#include <string_view>

#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief Returns "*" if \p refId is unmapped (< 0), empty string_view otherwise.
///
/// Callers should output the header reference name when the return value is empty.
[[nodiscard]] constexpr std::string_view RnameSentinel(std::int32_t refId) noexcept
{
    if (refId < 0) {
        return "*";
    }
    return {};
}

/// \brief Returns the SAM RNEXT sentinel for the given ref/next-ref combination.
///
/// Returns "*" when \p nextRefId is unmapped (< 0), "=" when it matches \p refId,
/// and an empty string_view when the header reference name should be used.
[[nodiscard]] constexpr std::string_view RnextSentinel(std::int32_t refId,
                                                       std::int32_t nextRefId) noexcept
{
    if (nextRefId < 0) {
        return "*";
    }
    if (nextRefId == refId) {
        return "=";
    }
    return {};
}

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_SAMFIELDUTILS_HPP
