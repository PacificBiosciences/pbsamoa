#include <pbsamoa/core/TagClipping.hpp>

#include <pbsamoa/core/Basemods.hpp>
#include <pbsamoa/core/Tags.hpp>

#include "PulseBitset.hpp"

#include <algorithm>
#include <array>
#include <iterator>
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

std::string SnapshotStringTag(const TagMap& tags, TagKey key)
{
    const TagValue* tagValue{tags.Get(key)};
    if (!tagValue) {
        return {};
    }

    const auto* tagString{std::get_if<std::string>(tagValue)};
    if (!tagString) {
        return {};
    }
    return *tagString;
}

template <typename Iterator>
std::size_t PrefixTotalBefore(Iterator begin, Iterator position)
{
    if (position == begin) {
        return 0;
    }
    return *std::prev(position);
}

bool ClearClipValue(TagValue& value)
{
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
        auto data = arr->MutableData();
        std::ranges::copy_n(std::data(data) + byteOffset, byteLength, std::data(data));
        arr->Resize(length);
        return true;
    }
    return false;
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
    // Guard against unsigned wrap when the clip window exceeds the sequence (mirrors the
    // bounds check in PileupClipStrategy); a corrupt window removes the tag.
    if ((clipOffset + clipLength) > seqLength) {
        return false;
    }
    const std::size_t reverseOffset{seqLength - (clipOffset + clipLength)};
    return ClipSubstring(value, reverseOffset, clipLength);
}

// --- PulseClipStrategy ---

bool PulseClipStrategy::Clip(TagValue& value, std::size_t clipOffset, std::size_t clipLength,
                             std::size_t /*seqLength*/, const ClipContext& ctx) const
{
    if (clipLength == 0) {
        return ClearClipValue(value);
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
        std::vector<BasemodRecord> clippedRecords;
        clippedRecords.reserve(std::size(records));

        for (const BasemodRecord& rec : records) {
            const BasemodClipWindow window{
                ClipBasemodRecord(rec, ctx.sequence, clipOffset, clipLength, ctx.isReverse)};
            clippedRecords.push_back(BasemodRecord{rec.Prefix, window.RetainedSkips});
        }

        *str = WriteBasemodString(clippedRecords);
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

        // ML must carry ModCodeCount(prefix) values per site across every record. A short
        // or mismatched ML would make the slice below run past the buffer (OOB read), so
        // validate the total up front and remove the tag rather than risk it.
        std::size_t expectedValues{0};
        for (const BasemodRecord& rec : records) {
            expectedValues += std::size(rec.Skips) * ModCodeCount(rec.Prefix);
        }
        if (arr->Count() < expectedValues) {
            return false;
        }

        // Build the retained ML values per type, then concatenate.
        std::vector<std::byte> retainedBytes;
        std::size_t qvOffset{0};  // running offset into the original ML array

        for (const BasemodRecord& rec : records) {
            const BasemodClipWindow window{
                ClipBasemodRecord(rec, ctx.sequence, clipOffset, clipLength, ctx.isReverse)};
            // ML stores ModCodeCount values per site (interleaved); each retained site
            // therefore carries `stride` contiguous ML values.
            const std::size_t stride{ModCodeCount(rec.Prefix)};
            // Copy retained QVs for this modification type
            const std::size_t srcByteOffset{(qvOffset + window.FrontRemoved * stride) * elemSize};
            const std::size_t srcByteLength{window.Retained * stride * elemSize};
            const auto srcData{arr->Data()};
            retainedBytes.insert(std::end(retainedBytes), std::cbegin(srcData) + srcByteOffset,
                                 std::cbegin(srcData) + srcByteOffset + srcByteLength);

            // Advance past all QVs for this modification type
            qvOffset += std::size(rec.Skips) * stride;
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
    if (!arr) {
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

    // Validate invariants — corrupt tags get removed rather than causing OOB
    // access
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

    // Find the first run that extends past the prefix (i.e., into the clip
    // window)
    const auto prefixIt{std::ranges::upper_bound(prefixSum, prefixSize)};
    // Find the first suffix run that extends past the suffix
    const auto suffixIt{std::ranges::upper_bound(suffixSum, suffixSize)};

    // Convert to pair indices in the original array
    const std::size_t beginRun{
        static_cast<std::size_t>(std::ranges::distance(prefixSum.begin(), prefixIt))};
    const std::size_t endRun{
        numRuns - static_cast<std::size_t>(std::ranges::distance(suffixSum.begin(), suffixIt))};

    // Compute how many bases from the first retained run are clipped off the
    // front
    const std::size_t lostPrefixBases{prefixSize - PrefixTotalBefore(prefixSum.begin(), prefixIt)};

    // Compute how many bases from the last retained run are clipped off the back
    const std::size_t lostSuffixBases{suffixSize - PrefixTotalBefore(suffixSum.begin(), suffixIt)};

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
                          std::size_t seqLength, std::string_view sequence, bool isReverse) const
{
    // Snapshot the pulse-call string before the loop so that all pulse
    // tags see the unmodified value regardless of iteration order.
    std::string pulseCalls{SnapshotStringTag(tags, TagKey{'p', 'c'})};

    // Snapshot the MM tag string before the loop so that ML processing
    // can reconstruct the original modification layout.
    std::string basemodString{SnapshotStringTag(tags, TagKey{'M', 'M'})};

    const ClipContext ctx{&tags, std::move(pulseCalls), sequence, std::move(basemodString),
                          isReverse};
    for (const auto& [key, strategy] : registrations_) {
        const TagValue* existing{tags.Get(key)};
        if (!existing) {
            continue;
        }
        TagValue val{*existing};
        if (strategy->Clip(val, clipOffset, clipLength, seqLength, ctx)) {
            tags.Set(key, std::move(val));
        } else {
            tags.Remove(key);
        }
    }

    // MN:i records the SEQ length MM/ML were produced against; htslib rejects records
    // whose MN != l_qseq (sam_mods.c). Keep it consistent with the clipped sequence.
    static constexpr TagKey MN_TAG{'M', 'N'};
    if (tags.Contains(MN_TAG)) {
        tags.Set(MN_TAG, static_cast<std::int64_t>(clipLength));
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
