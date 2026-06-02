#include <pbsamoa/core/Basemods.hpp>

#include <algorithm>
#include <array>
#include <charconv>
#include <iterator>
#include <ranges>
#include <string>
#include <string_view>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {
namespace {

std::string_view TakeDelimitedSegment(std::string_view& text, char delimiter)
{
    const std::size_t delimiterPos{text.find(delimiter)};
    if (delimiterPos == std::string_view::npos) {
        const std::string_view segment{text};
        text = {};
        return segment;
    }
    const std::string_view segment{text.substr(0, delimiterPos)};
    text = text.substr(delimiterPos + 1);
    return segment;
}

}  // namespace

std::vector<BasemodRecord> ParseBasemodString(std::string_view mm)
{
    std::vector<BasemodRecord> records;

    while (!std::empty(mm)) {
        const std::string_view segment{TakeDelimitedSegment(mm, ';')};

        if (std::size(segment) < 3) {
            // Minimum valid: "C+m" (3 chars)
            continue;
        }

        // Prefix: BASE OP MOD_CODE [.?]. MOD_CODE is either a numeric ChEBI id
        // (e.g. "C+76792") or a letter code (e.g. "C+m"); per htslib (sam_mods.c:298-312)
        // a digit run is a single ChEBI code rather than a base + skip count.
        std::size_t prefixLen{2};  // BASE + OP
        if ((prefixLen < std::size(segment)) && (segment[prefixLen] >= '0') &&
            (segment[prefixLen] <= '9')) {
            while ((prefixLen < std::size(segment)) && (segment[prefixLen] >= '0') &&
                   (segment[prefixLen] <= '9')) {
                ++prefixLen;
            }
        } else if (prefixLen < std::size(segment)) {
            ++prefixLen;  // single-letter mod code
        }
        if ((prefixLen < std::size(segment)) &&
            ((segment[prefixLen] == '?') || (segment[prefixLen] == '.'))) {
            ++prefixLen;
        }

        BasemodRecord rec;
        rec.Prefix = std::string{segment.substr(0, prefixLen)};

        std::size_t pos{prefixLen};
        while (pos < std::size(segment)) {
            if (segment[pos] == ',') {
                ++pos;
            }
            std::int32_t num{0};
            const char* first{std::data(segment) + pos};
            const char* last{std::data(segment) + std::size(segment)};
            const auto [ptr, ec]{std::from_chars(first, last, num)};
            if (ec == std::errc{}) {
                rec.Skips.push_back(num);
                pos += static_cast<std::size_t>(ptr - first);
            } else {
                break;
            }
        }

        records.push_back(std::move(rec));
    }

    return records;
}

std::string WriteBasemodString(std::span<const BasemodRecord> records)
{
    std::string out;
    std::array<char, 16> numBuf{};
    for (const BasemodRecord& rec : records) {
        out.append(rec.Prefix);
        for (const std::int32_t s : rec.Skips) {
            out += ',';
            const std::to_chars_result toCharsResult{
                std::to_chars(std::data(numBuf), std::data(numBuf) + std::size(numBuf), s)};
            out.append(std::data(numBuf), toCharsResult.ptr);
        }
        out += ';';
    }
    return out;
}

BasemodClipWindow ClipBasemodRecord(const BasemodRecord& record, std::string_view sequence,
                                    std::size_t clipOffset, std::size_t clipLength)
{
    BasemodClipWindow window;
    if (std::empty(record.Prefix)) {
        return window;
    }
    const char canonicalBase{record.Prefix[0]};

    // Count canonical bases before and within the clip window. The 'N' wildcard matches
    // every base (htslib freq[15] = l_qseq), not the literal character 'N'.
    const auto countBases{[canonicalBase](std::string_view s) -> std::int32_t {
        if (canonicalBase == 'N') {
            return static_cast<std::int32_t>(std::ssize(s));
        }
        return static_cast<std::int32_t>(std::ranges::count(s, canonicalBase));
    }};
    const std::int32_t basesBeforeClip{countBases(sequence.substr(0, clipOffset))};
    const std::int32_t basesInClip{countBases(sequence.substr(clipOffset, clipLength))};

    // prefixSum[i] = total canonical bases consumed up to and including site i
    // (skip value s consumes s+1 canonical bases to reach the next mod).
    std::vector<std::int32_t> prefixSum;
    prefixSum.reserve(std::size(record.Skips));
    std::int32_t pSum{0};
    for (const std::int32_t s : record.Skips) {
        pSum += (s + 1);
        prefixSum.push_back(pSum);
    }

    const auto startIt{std::ranges::lower_bound(prefixSum, basesBeforeClip + 1)};
    const auto endIt{std::ranges::upper_bound(prefixSum, basesBeforeClip + basesInClip)};
    const std::size_t startIdx{
        static_cast<std::size_t>(std::ranges::distance(std::begin(prefixSum), startIt))};
    const std::size_t endIdx{
        static_cast<std::size_t>(std::ranges::distance(std::begin(prefixSum), endIt))};

    window.FrontRemoved = startIdx;
    window.Retained = endIdx - startIdx;
    // Canonical bases consumed by the front-removed sites; the difference from
    // basesBeforeClip is what restore adds back to the first retained skip.
    const std::int32_t prefixAtLastFront{(startIdx > 0) ? prefixSum[startIdx - 1] : 0};
    window.PrefixLost = basesBeforeClip - prefixAtLastFront;
    if (window.Retained > 0) {
        window.RetainedSkips.assign(std::cbegin(record.Skips) + startIdx,
                                    std::cbegin(record.Skips) + endIdx);
        // The first retained skip counts canonical bases between the clip start
        // and the first retained mod site.
        window.RetainedSkips[0] = prefixSum[startIdx] - basesBeforeClip - 1;
    }
    return window;
}

}  // namespace Samoa
}  // namespace PacBio
