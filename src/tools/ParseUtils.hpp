#ifndef PBSAMOA_TOOLS_PARSEUTILS_HPP
#define PBSAMOA_TOOLS_PARSEUTILS_HPP

#include <algorithm>
#include <charconv>
#include <concepts>
#include <expected>
#include <format>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {
namespace Tools {

struct ParsedRegion
{
    std::string RefName;
    std::int32_t Beg{0};
    std::int32_t End{0};
};

template <std::integral IntType>
[[nodiscard]] std::expected<IntType, std::string> ParseInteger(std::string_view text,
                                                               std::string_view name)
{
    IntType value{};
    const char* const begin{text.data()};
    const char* const end{text.data() + text.size()};
    const auto [ptr, ec]{std::from_chars(begin, end, value)};

    if ((ec != std::errc{}) || (ptr != end)) {
        return std::unexpected{std::format("invalid {}: {}", name, text)};
    }
    return value;
}

struct ParsedRegionFields
{
    std::string_view RefName;
    std::string_view StartText;
    std::string_view EndText;
};

[[nodiscard]] inline std::expected<ParsedRegionFields, std::string> SplitRegionFields(
    std::string_view text, std::string_view /*expectedFormat*/)
{
    // Split on the rightmost colon so reference names that themselves contain ':'
    // (e.g. "HLA-DRB1*12:17") are handled, matching htslib hts_parse_region (hts.c).
    const std::size_t colonPos{text.rfind(':')};
    if (colonPos == std::string_view::npos) {
        // Bare reference name -> whole reference (htslib accepts this).
        return ParsedRegionFields{.RefName = text, .StartText = {}, .EndText = {}};
    }

    const std::string_view refName{text.substr(0, colonPos)};
    const std::string_view rest{text.substr(colonPos + 1)};
    const std::size_t dashPos{rest.find('-')};
    if (dashPos == std::string_view::npos) {
        // "chr:100" -> from start to the end of the reference.
        return ParsedRegionFields{.RefName = refName, .StartText = rest, .EndText = {}};
    }

    return ParsedRegionFields{
        .RefName = refName,
        .StartText = rest.substr(0, dashPos),
        .EndText = rest.substr(dashPos + 1),
    };
}

[[nodiscard]] inline std::expected<ParsedRegion, std::string> BuildParsedRegion(
    const ParsedRegionFields& fields)
{
    // Open-ended end runs to the end of the reference (htslib uses HTS_POS_MAX; here the
    // coordinate type is 32-bit, so INT32_MAX is the sentinel Query treats as "no upper bound").
    constexpr std::int32_t OPEN_END{std::numeric_limits<std::int32_t>::max()};

    // Start: empty or "0" -> whole-reference start; otherwise a 1-based coordinate.
    std::int32_t beg{0};
    if (!std::empty(fields.StartText)) {
        const auto parsedStart{ParseInteger<std::int32_t>(fields.StartText, "region start")};
        if (!parsedStart) {
            return std::unexpected{parsedStart.error()};
        }
        if (*parsedStart < 0) {
            return std::unexpected{"region start must be >= 0"};
        }
        beg = (*parsedStart > 0) ? (*parsedStart - 1) : 0;
    }

    // End: empty -> open; otherwise a 1-based inclusive end == 0-based exclusive end.
    std::int32_t end{OPEN_END};
    if (!std::empty(fields.EndText)) {
        const auto parsedEnd{ParseInteger<std::int32_t>(fields.EndText, "region end")};
        if (!parsedEnd) {
            return std::unexpected{parsedEnd.error()};
        }
        // beg is already 0-based, parsedEnd is the 1-based inclusive end, so the smallest
        // valid end equals beg + 1 (a single-base range). Reject end <= beg rather than
        // end < beg, which would let an end one base before the start become an empty range.
        if (*parsedEnd <= beg) {
            return std::unexpected{"region end must be >= start"};
        }
        end = *parsedEnd;
    }

    return ParsedRegion{.RefName = std::string{fields.RefName}, .Beg = beg, .End = end};
}

[[nodiscard]] inline std::expected<ParsedRegion, std::string> ParseRegion(
    std::string_view text, std::string_view expectedFormat)
{
    return SplitRegionFields(text, expectedFormat).and_then(BuildParsedRegion);
}

template <typename IntType>
IntType ParseIntegerOrThrow(std::string_view text, std::string_view name)
{
    const auto parsed{ParseInteger<IntType>(text, name)};
    if (!parsed) {
        throw std::runtime_error{parsed.error()};
    }
    return *parsed;
}

/// \brief Resolve a user-requested worker count to an actual thread count.
///
/// \param requested  User-requested count (< 0 means auto-detect)
/// \param explicitCap  Maximum when user provides explicit count
/// \param autoCap  Maximum when auto-detecting (default: 8)
inline std::size_t ResolveNumWorkers(std::int32_t requested, std::int32_t explicitCap,
                                     std::int32_t autoCap = 8)
{
    const std::uint32_t hardwareThreadCount{std::thread::hardware_concurrency()};
    const std::int32_t hardwareThreads{
        static_cast<std::int32_t>(std::ranges::max(hardwareThreadCount, 1U))};

    if (requested < 0) {
        return static_cast<std::size_t>(std::ranges::min(hardwareThreads, autoCap));
    }
    return static_cast<std::size_t>(std::ranges::min(requested, explicitCap));
}

}  // namespace Tools
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_PARSEUTILS_HPP
