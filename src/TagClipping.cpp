#include <pbsamoa/core/TagClipping.hpp>

#include <pbsamoa/core/Tags.hpp>

#include "PulseBitset.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <functional>
#include <memory>
#include <numeric>
#include <ranges>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

namespace {

bool ClipSubstring(TagValue& value, std::size_t offset, std::size_t length)
{
    if (auto* str = std::get_if<std::string>(&value)) {
        *str = str->substr(offset, length);
        return true;
    }
    if (auto* arr = std::get_if<TagArray>(&value)) {
        const std::size_t elemSize{arr->ElementSize()};
        if (elemSize == 0) {
            return false;
        }
        const std::size_t byteOffset{offset * elemSize};
        const std::size_t byteLength{length * elemSize};
        auto& data = arr->MutableData();
        std::ranges::copy_n(std::data(data) + byteOffset, byteLength, std::data(data));
        arr->Resize(length);
        return true;
    }
    return false;
}

/// \brief A single modification type parsed from an MM string.
///
/// E.g., "C+m?,1,3,0" -> prefix = "C+m?", skips = {1, 3, 0}
struct BasemodRecord
{
    std::string prefix;  // e.g. "C+m", "C+m?", "A+a."
    std::vector<std::int32_t> skips;
};

/// \brief Parse an MM tag string into ordered modification records.
///
/// Format: "C+m,1,3,0;A+a,2;" -> [{prefix="C+m", skips={1,3,0}}, {prefix="A+a", skips={2}}]
std::vector<BasemodRecord> ParseBasemodString(std::string_view mm)
{
    std::vector<BasemodRecord> records;

    while (!std::empty(mm)) {
        // Find the semicolon that terminates this modification record
        const std::size_t semi{mm.find(';')};
        const std::string_view segment{(semi != std::string_view::npos) ? mm.substr(0, semi) : mm};

        if (std::size(segment) < 3) {
            // Minimum valid: "C+m" (3 chars)
            mm = (semi != std::string_view::npos) ? mm.substr(semi + 1) : std::string_view{};
            continue;
        }

        // Parse prefix: BASE OP MOD_CODE [.?]
        // First char is base, second is +/-, third+ is mod code(s)
        // Optional trailing '.' or '?' before the first comma
        std::size_t prefixLen{3};
        if ((prefixLen < std::size(segment)) &&
            ((segment[prefixLen] == '?') || (segment[prefixLen] == '.'))) {
            ++prefixLen;
        }

        BasemodRecord rec;
        rec.prefix = std::string{segment.substr(0, prefixLen)};

        // Parse comma-separated skip counts after the prefix
        std::size_t pos{prefixLen};
        while (pos < std::size(segment)) {
            if (segment[pos] == ',') {
                ++pos;
            }
            std::int32_t num{0};
            const auto* first{std::data(segment) + pos};
            const auto* last{std::data(segment) + std::size(segment)};
            const auto [ptr, ec]{std::from_chars(first, last, num)};
            if (ec == std::errc{}) {
                rec.skips.push_back(num);
                pos += static_cast<std::size_t>(ptr - first);
            } else {
                break;
            }
        }

        records.push_back(std::move(rec));
        mm = (semi != std::string_view::npos) ? mm.substr(semi + 1) : std::string_view{};
    }

    return records;
}

/// \brief For a single modification type, compute how many modification
/// sites fall before, within, and after a clip window.
///
/// Uses the same prefix-sum algorithm as pbbam's ClipBasemodsTag.
///
/// \param[in] canonicalBase  the base to count (e.g., 'C')
/// \param[in] skips          original skip counts
/// \param[in] sequence       full original sequence
/// \param[in] clipOffset     start of clip window
/// \param[in] clipLength     length of clip window
/// \param[out] frontRemoved  number of mod sites before clip window
/// \param[out] retained      number of mod sites within clip window
/// \param[out] newSkips      rewritten skip counts for retained sites
void ClipSingleModType(char canonicalBase, std::span<const std::int32_t> skips,
                       std::string_view sequence, std::size_t clipOffset, std::size_t clipLength,
                       std::size_t& frontRemoved, std::size_t& retained,
                       std::vector<std::int32_t>& newSkips)
{
    // Count canonical bases before and within the clip window
    const std::int32_t basesBeforeClip =
        std::ranges::count(sequence.substr(0, clipOffset), canonicalBase);
    const std::int32_t basesInClip =
        std::ranges::count(sequence.substr(clipOffset, clipLength), canonicalBase);

    // Build prefix sums: prefixSum[i] = total canonical bases seen up to
    // and including modification site i.
    // For skip value s, we pass (s+1) canonical bases to reach the next mod.
    std::vector<std::int32_t> prefixSum;
    prefixSum.reserve(std::size(skips));
    std::int32_t pSum{0};
    for (const std::int32_t s : skips) {
        pSum += (s + 1);
        prefixSum.push_back(pSum);
    }

    // Find the range of modification sites within the clip window.
    // A mod site at prefixSum[i] is "in" the clip if:
    //   basesBeforeClip < prefixSum[i] <= basesBeforeClip + basesInClip
    const auto startIt{std::ranges::lower_bound(prefixSum, basesBeforeClip + 1)};
    const auto endIt{std::ranges::upper_bound(prefixSum, basesBeforeClip + basesInClip)};

    const std::size_t startIdx = std::ranges::distance(prefixSum.begin(), startIt);
    const std::size_t endIdx = std::ranges::distance(prefixSum.begin(), endIt);

    frontRemoved = startIdx;
    retained = endIdx - startIdx;

    newSkips.clear();
    if (retained > 0) {
        // Copy the retained skip counts
        newSkips.assign(std::cbegin(skips) + startIdx, std::cbegin(skips) + endIdx);

        // Adjust the first retained skip count: it should reflect the number
        // of canonical bases between the clip window start and the first
        // retained modification site.
        newSkips[0] = prefixSum[startIdx] - basesBeforeClip - 1;
    }
}

/// \brief Reconstruct an MM string from clipped modification records.
std::string RewriteBasemodString(std::span<const BasemodRecord> records,
                                 std::span<const std::vector<std::int32_t>> clippedSkips)
{
    std::string out;
    std::array<char, 16> numBuf{};
    for (std::size_t i{0}; i < std::size(records); ++i) {
        out.append(records[i].prefix);
        for (const std::int32_t s : clippedSkips[i]) {
            out += ',';
            const auto [ptr, ec]{
                std::to_chars(std::data(numBuf), std::data(numBuf) + std::size(numBuf), s)};
            out.append(std::data(numBuf), ptr);
        }
        out += ';';
    }
    return out;
}

}  // namespace

// --- SubstringClipStrategy ---

bool SubstringClipStrategy::Clip(TagValue& value, std::size_t clipOffset, std::size_t clipLength,
                                 std::size_t /*seqLength*/, const ClipContext& /*ctx*/) const
{
    return ClipSubstring(value, clipOffset, clipLength);
}

// --- ReverseSubstringClipStrategy ---

bool ReverseSubstringClipStrategy::Clip(TagValue& value, std::size_t clipOffset,
                                        std::size_t clipLength, std::size_t seqLength,
                                        const ClipContext& /*ctx*/) const
{
    const std::size_t reverseOffset{seqLength - (clipOffset + clipLength)};
    return ClipSubstring(value, reverseOffset, clipLength);
}

// --- PulseClipStrategy ---

bool PulseClipStrategy::Clip(TagValue& value, std::size_t clipOffset, std::size_t clipLength,
                             std::size_t /*seqLength*/, const ClipContext& ctx) const
{
    if (clipLength == 0) {
        if (auto* str = std::get_if<std::string>(&value)) {
            str->clear();
            return true;
        }
        if (auto* arr = std::get_if<TagArray>(&value)) {
            arr->Resize(0);
            return true;
        }
        return false;
    }

    if (std::empty(ctx.pulseCalls)) {
        return false;
    }

    const PulseBitset bitset{ctx.pulseCalls};

    const std::size_t startPulse{bitset.FindNthBase(clipOffset)};
    const std::size_t endPulse{bitset.FindNthBase(clipOffset + clipLength - 1)};
    if ((startPulse == PulseBitset::NPOS) || (endPulse == PulseBitset::NPOS)) {
        return false;
    }

    const std::size_t pulseLength{endPulse - startPulse + 1};
    return ClipSubstring(value, startPulse, pulseLength);
}

// --- BasemodClipStrategy ---

bool BasemodClipStrategy::Clip(TagValue& value, std::size_t clipOffset, std::size_t clipLength,
                               std::size_t /*seqLength*/, const ClipContext& ctx) const
{
    // --- Handle MM tag (string) ---
    if (auto* str = std::get_if<std::string>(&value)) {
        if (std::empty(ctx.sequence)) {
            return false;
        }

        const std::vector<BasemodRecord> records{ParseBasemodString(*str)};
        std::vector<std::vector<std::int32_t>> allClippedSkips;
        allClippedSkips.reserve(std::size(records));

        for (const auto& rec : records) {
            const char canonicalBase{rec.prefix[0]};
            std::size_t frontRemoved{0};
            std::size_t retained{0};
            std::vector<std::int32_t> newSkips;

            ClipSingleModType(canonicalBase, rec.skips, ctx.sequence, clipOffset, clipLength,
                              frontRemoved, retained, newSkips);
            allClippedSkips.push_back(std::move(newSkips));
        }

        *str = RewriteBasemodString(records, allClippedSkips);
        return true;
    }

    // --- Handle ML tag (B:C uint8 array) ---
    if (auto* arr = std::get_if<TagArray>(&value)) {
        if (std::empty(ctx.sequence) || std::empty(ctx.basemodString)) {
            return false;
        }

        const std::size_t elemSize{arr->ElementSize()};
        if (elemSize == 0) {
            return false;
        }

        // Parse the ORIGINAL MM string to determine modification layout.
        // ML entries are stored in modification-type order: all QVs for the
        // first type, then all for the second, etc.
        const std::vector<BasemodRecord> records{ParseBasemodString(ctx.basemodString)};

        // Build the retained ML values per type, then concatenate.
        std::vector<std::byte> retainedBytes;
        std::size_t qvOffset{0};  // running offset into the original ML array

        for (const auto& rec : records) {
            const char canonicalBase{rec.prefix[0]};
            std::size_t frontRemoved{0};
            std::size_t retained{0};
            std::vector<std::int32_t> newSkips;

            ClipSingleModType(canonicalBase, rec.skips, ctx.sequence, clipOffset, clipLength,
                              frontRemoved, retained, newSkips);

            // Copy retained QVs for this modification type
            const std::size_t srcByteOffset{(qvOffset + frontRemoved) * elemSize};
            const std::size_t srcByteLength{retained * elemSize};
            const auto srcData{arr->Data()};
            retainedBytes.insert(std::end(retainedBytes), std::cbegin(srcData) + srcByteOffset,
                                 std::cbegin(srcData) + srcByteOffset + srcByteLength);

            // Advance past all QVs for this modification type
            qvOffset += std::size(rec.skips);
        }

        // Write retained bytes back
        const std::size_t totalRetained{std::size(retainedBytes) / elemSize};
        if (totalRetained > 0) {
            std::ranges::copy_n(std::data(retainedBytes), std::size(retainedBytes),
                                std::data(arr->MutableData()));
        }
        arr->Resize(static_cast<std::uint32_t>(totalRetained));
        return true;
    }

    return false;
}

// --- PileupClipStrategy ---

bool PileupClipStrategy::Clip(TagValue& value, std::size_t clipOffset, std::size_t clipLength,
                              std::size_t seqLength, const ClipContext& /*ctx*/) const
{
    auto* arr = std::get_if<TagArray>(&value);
    if (arr == nullptr) {
        return false;
    }

    const std::size_t elemSize{arr->ElementSize()};
    if (elemSize != 1) {
        return false;
    }

    const std::uint32_t count{arr->Count()};

    // Empty clip -> empty result
    if (clipLength == 0) {
        arr->Resize(0);
        return true;
    }

    // Must have an even number of elements (pairs of runLen, coverage)
    if ((count % 2) != 0) {
        return false;
    }

    const std::size_t numRuns{count / 2};
    const auto rawData{arr->Data()};

    // Extract run lengths into a vector for prefix-sum computation
    std::vector<std::size_t> lengths;
    lengths.reserve(numRuns);
    for (std::size_t i{0}; i < numRuns; ++i) {
        lengths.push_back(static_cast<std::uint8_t>(rawData[i * 2]));
    }

    // Validate invariants — corrupt tags get removed rather than causing OOB access
    if (std::ranges::fold_left(lengths, std::size_t{0}, std::plus{}) != seqLength) {
        return false;
    }
    if ((clipOffset > seqLength) || ((clipOffset + clipLength) > seqLength)) {
        return false;
    }

    // Build prefix sums over run lengths
    std::vector<std::size_t> prefixSum(numRuns);
    std::inclusive_scan(std::cbegin(lengths), std::cend(lengths), std::begin(prefixSum));

    // Build suffix sums (reverse inclusive scan)
    std::vector<std::size_t> suffixSum(numRuns);
    std::inclusive_scan(std::crbegin(lengths), std::crend(lengths), std::begin(suffixSum));

    const std::size_t clipEnd{clipOffset + clipLength};
    const std::size_t prefixSize{clipOffset};
    const std::size_t suffixSize{seqLength - clipEnd};

    // Find the first run that extends past the prefix (i.e., into the clip window)
    const auto prefixIt{std::ranges::upper_bound(prefixSum, prefixSize)};
    // Find the first suffix run that extends past the suffix
    const auto suffixIt{std::ranges::upper_bound(suffixSum, suffixSize)};

    // Convert to pair indices in the original array
    const std::size_t beginRun = std::ranges::distance(prefixSum.begin(), prefixIt);
    const std::size_t endRun{
        numRuns - static_cast<std::size_t>(std::ranges::distance(suffixSum.begin(), suffixIt))};

    // Compute how many bases from the first retained run are clipped off the front
    const std::size_t lostPrefixBases{prefixSize -
                                      ((prefixIt != prefixSum.begin()) ? *std::prev(prefixIt) : 0)};

    // Compute how many bases from the last retained run are clipped off the back
    const std::size_t lostSuffixBases{suffixSize -
                                      ((suffixIt != suffixSum.begin()) ? *std::prev(suffixIt) : 0)};

    // Build the new RLE pairs
    TagArray result{'C'};
    for (std::size_t i{beginRun}; i < endRun; ++i) {
        std::uint8_t runLen{static_cast<std::uint8_t>(rawData[i * 2])};
        const std::uint8_t coverage{static_cast<std::uint8_t>(rawData[i * 2 + 1])};

        // Adjust first run: trim prefix
        if ((i == beginRun) && (lostPrefixBases != 0)) {
            runLen -= static_cast<std::uint8_t>(lostPrefixBases);
        }
        // Adjust last run: trim suffix
        if ((i == (endRun - 1)) && (lostSuffixBases != 0)) {
            runLen -= static_cast<std::uint8_t>(lostSuffixBases);
        }

        // Skip zero-length runs (can happen if clip aligns exactly on run boundary
        // and both prefix and suffix eat the same single-pair run to zero)
        if (runLen == 0) {
            continue;
        }

        result.AppendUInt8(runLen);
        result.AppendUInt8(coverage);
    }

    *arr = std::move(result);
    return true;
}

// --- TagClipper ---

void TagClipper::Register(std::initializer_list<TagKey> tags,
                          std::shared_ptr<TagClipStrategy> strategy)
{
    for (const auto& key : tags) {
        registrations_.push_back(Registration{key, strategy});
    }
}

void TagClipper::ClipTags(TagMap& tags, std::size_t clipOffset, std::size_t clipLength,
                          std::size_t seqLength, std::string_view sequence) const
{
    // Snapshot the pulse-call string before the loop so that all pulse
    // tags see the unmodified value regardless of iteration order.
    std::string pulseCalls;
    const TagValue* pcVal{tags.Get(TagKey{'p', 'c'})};
    if (pcVal != nullptr) {
        const auto* pcStr = std::get_if<std::string>(pcVal);
        if (pcStr != nullptr) {
            pulseCalls = *pcStr;
        }
    }

    // Snapshot the MM tag string before the loop so that ML processing
    // can reconstruct the original modification layout.
    std::string basemodString;
    const TagValue* mmVal{tags.Get(TagKey{'M', 'M'})};
    if (mmVal != nullptr) {
        const auto* mmStr = std::get_if<std::string>(mmVal);
        if (mmStr != nullptr) {
            basemodString = *mmStr;
        }
    }

    const ClipContext ctx{&tags, std::move(pulseCalls), sequence, std::move(basemodString)};
    for (const auto& [key, strategy] : registrations_) {
        const TagValue* existing{tags.Get(key)};
        if (existing == nullptr) {
            continue;
        }
        TagValue val{*existing};
        if (strategy->Clip(val, clipOffset, clipLength, seqLength, ctx)) {
            tags.Set(key, std::move(val));
        } else {
            tags.Remove(key);
        }
    }
}

TagClipper TagClipper::PacBioDefault()
{
    TagClipper clipper;

    auto substringStrategy{std::make_shared<SubstringClipStrategy>()};
    clipper.Register(
        {
            TagKey{'d', 'q'},
            TagKey{'i', 'q'},
            TagKey{'m', 'q'},
            TagKey{'s', 'q'},
            TagKey{'d', 't'},
            TagKey{'s', 't'},
            TagKey{'i', 'p'},
            TagKey{'p', 'w'},
            TagKey{'f', 'i'},
            TagKey{'f', 'p'},
        },
        substringStrategy);

    auto reverseStrategy{std::make_shared<ReverseSubstringClipStrategy>()};
    clipper.Register(
        {
            TagKey{'r', 'i'},
            TagKey{'r', 'p'},
        },
        reverseStrategy);

    auto pulseStrategy{std::make_shared<PulseClipStrategy>()};
    clipper.Register(
        {
            TagKey{'p', 'c'},
            TagKey{'p', 't'},
            TagKey{'p', 'q'},
            TagKey{'p', 'v'},
            TagKey{'p', 'g'},
            TagKey{'p', 'a'},
            TagKey{'p', 'm'},
            TagKey{'p', 's'},
            TagKey{'p', 'i'},
            TagKey{'p', 'd'},
            TagKey{'p', 'x'},
            TagKey{'p', 'e'},
            TagKey{'s', 'f'},
        },
        pulseStrategy);

    auto basemodStrategy{std::make_shared<BasemodClipStrategy>()};
    clipper.Register(
        {
            TagKey{'M', 'M'},
            TagKey{'M', 'L'},
        },
        basemodStrategy);

    // Subread pileup: sm/sx are per-position arrays (same as substring)
    clipper.Register(
        {
            TagKey{'s', 'm'},
            TagKey{'s', 'x'},
        },
        substringStrategy);

    // Subread pileup: sa is RLE-encoded coverage
    auto pileupStrategy{std::make_shared<PileupClipStrategy>()};
    clipper.Register(
        {
            TagKey{'s', 'a'},
        },
        pileupStrategy);

    return clipper;
}

}  // namespace Samoa
}  // namespace PacBio
