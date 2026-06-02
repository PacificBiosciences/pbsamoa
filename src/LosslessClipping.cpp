#include <pbsamoa/core/LosslessClipping.hpp>

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/Basemods.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <pbcopper/json/JSON.h>

#include <algorithm>
#include <array>
#include <iterator>
#include <map>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include <cstdint>
#include <cstring>

namespace PacBio {
namespace Samoa {
namespace {

constexpr TagKey LS_TAG{'l', 's'};

struct KineticTag
{
    TagKey Key;
    const char* Field;  ///< `ls` JSON field name (== tag name)
    bool Reverse;       ///< stored mirrored relative to the sequence (rp/ri)
};

constexpr std::array<KineticTag, 6> KINETIC_TAGS{{
    {TagKey{'p', 'w'}, "pw", false},
    {TagKey{'i', 'p'}, "ip", false},
    {TagKey{'f', 'p'}, "fp", false},
    {TagKey{'f', 'i'}, "fi", false},
    {TagKey{'r', 'p'}, "rp", true},
    {TagKey{'r', 'i'}, "ri", true},
}};

// Demultiplexing tags removed when a record is restored to its pre-clip state.
constexpr std::array<TagKey, 8> DEMUX_TAGS{{
    TagKey{'b', 'c'},
    TagKey{'b', 'q'},
    TagKey{'b', 'x'},
    TagKey{'b', 'l'},
    TagKey{'b', 't'},
    TagKey{'q', 'l'},
    TagKey{'q', 't'},
    LS_TAG,
}};

const TagArray* GetArray(const TagMap& tags, TagKey key)
{
    const TagValue* value{tags.Get(key)};
    if ((value == nullptr) || !std::holds_alternative<TagArray>(*value)) {
        return nullptr;
    }
    return &std::get<TagArray>(*value);
}

// Decode a numeric tag array to int64 values, widening from its on-disk element
// type. memcpy is the well-defined way to read a typed value out of raw bytes.
std::vector<std::int64_t> DecodeIntArray(const TagArray& array)
{
    const std::span<const std::byte> data{array.Data()};
    const std::uint32_t count{array.Count()};
    std::vector<std::int64_t> out;
    out.reserve(count);
    const auto decode = [&](auto sample) {
        using Element = decltype(sample);
        for (std::uint32_t i = 0; i < count; ++i) {
            Element value{};
            std::memcpy(&value, std::data(data) + (i * sizeof(Element)), sizeof(Element));
            out.push_back(static_cast<std::int64_t>(value));
        }
    };
    switch (array.ElementType()) {
        case 'c':
            decode(std::int8_t{});
            break;
        case 'C':
            decode(std::uint8_t{});
            break;
        case 's':
            decode(std::int16_t{});
            break;
        case 'S':
            decode(std::uint16_t{});
            break;
        case 'i':
            decode(std::int32_t{});
            break;
        case 'I':
            decode(std::uint32_t{});
            break;
        default:
            break;
    }
    return out;
}

TagArray EncodeIntArray(char elementType, const std::vector<std::int64_t>& values)
{
    TagArray array{elementType};
    for (const std::int64_t value : values) {
        switch (elementType) {
            case 'c':
                array.AppendInt8(static_cast<std::int8_t>(value));
                break;
            case 'C':
                array.AppendUInt8(static_cast<std::uint8_t>(value));
                break;
            case 's':
                array.AppendInt16(static_cast<std::int16_t>(value));
                break;
            case 'S':
                array.AppendUInt16(static_cast<std::uint16_t>(value));
                break;
            case 'i':
                array.AppendInt32(static_cast<std::int32_t>(value));
                break;
            case 'I':
                array.AppendUInt32(static_cast<std::uint32_t>(value));
                break;
            default:
                break;
        }
    }
    return array;
}

std::vector<std::uint8_t> RawBytes(const TagArray& array)
{
    const std::span<const std::byte> data{array.Data()};
    std::vector<std::uint8_t> out;
    out.reserve(std::size(data));
    for (const std::byte b : data) {
        out.push_back(static_cast<std::uint8_t>(b));
    }
    return out;
}

std::string FastqFromQualities(std::span<const std::uint8_t> quals)
{
    std::string fastq;
    fastq.reserve(std::size(quals));
    for (const std::uint8_t qv : quals) {
        fastq.push_back(static_cast<char>(qv + 33));
    }
    return fastq;
}

std::vector<std::uint8_t> QualitiesFromFastq(std::string_view fastq)
{
    std::vector<std::uint8_t> quals;
    quals.reserve(std::size(fastq));
    for (const char c : fastq) {
        quals.push_back(static_cast<std::uint8_t>(c - 33));
    }
    return quals;
}

// Per-base tags that the tag clipper trims but lossless storage does not yet
// capture (raw-subread pulse + subread-pileup tags). Rejected up front so a clip
// can never silently drop data that undo would need.
constexpr std::array<TagKey, 16> UNSUPPORTED_CLIP_TAGS{{
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
    TagKey{'s', 'm'},
    TagKey{'s', 'x'},
    TagKey{'s', 'a'},
}};

// Capture the removed base-modification calls (MM/ML) into the `ls` lead/trail
// maps, plus the per-mod prefix-lost-base count (pMM) needed to recover the first
// retained skip on undo. Mirrors lima's format. MM/ML are per-mod-type maps keyed
// by the modification prefix (e.g. "C+m?"); ML quals follow MM-record order.
void CaptureBasemods(const TagMap& tags, std::string_view sequence, std::int32_t clipLeft,
                     std::int32_t clipRight, JSON::Json& lead, JSON::Json& trail, JSON::Json& pMM)
{
    const TagValue* mmValue{tags.Get(TagKey{'M', 'M'})};
    const std::string* mm{mmValue ? std::get_if<std::string>(mmValue) : nullptr};
    if (mm == nullptr) {
        return;
    }
    const std::vector<BasemodRecord> records{ParseBasemodString(*mm)};

    std::vector<std::int64_t> mlValues;
    if (const TagArray* mlArray{GetArray(tags, TagKey{'M', 'L'})}; mlArray != nullptr) {
        mlValues = DecodeIntArray(*mlArray);
    }
    std::size_t totalMlValues{0};
    for (const BasemodRecord& rec : records) {
        // ML carries ModCodeCount values per site (interleaved); multi-letter mods (e.g.
        // "C+mh") therefore consume more ML entries than skip counts.
        totalMlValues += std::size(rec.Skips) * ModCodeCount(rec.Prefix);
    }
    // Require ML to cover every MM site before slicing it; a short or malformed ML
    // is treated as absent rather than risking an out-of-range read.
    const bool hasMl{!std::empty(mlValues) && (std::size(mlValues) >= totalMlValues)};
    const std::size_t clipLength{static_cast<std::size_t>(clipRight - clipLeft)};

    JSON::Json leadMM = JSON::Json::object();
    JSON::Json leadML = JSON::Json::object();
    JSON::Json trailMM = JSON::Json::object();
    JSON::Json trailML = JSON::Json::object();

    std::size_t qvOffset{0};
    for (const BasemodRecord& rec : records) {
        const BasemodClipWindow window{
            ClipBasemodRecord(rec, sequence, static_cast<std::size_t>(clipLeft), clipLength)};
        const std::size_t total{std::size(rec.Skips)};
        const std::size_t frontN{window.FrontRemoved};
        const std::size_t trailStart{window.FrontRemoved + window.Retained};
        // ML stride: each site occupies this many (interleaved) ML values.
        const std::size_t stride{ModCodeCount(rec.Prefix)};

        if (frontN > 0) {
            leadMM[rec.Prefix] =
                std::vector<std::int32_t>(std::cbegin(rec.Skips), std::cbegin(rec.Skips) + frontN);
            if (hasMl) {
                leadML[rec.Prefix] =
                    std::vector<std::int64_t>(std::cbegin(mlValues) + qvOffset,
                                              std::cbegin(mlValues) + qvOffset + frontN * stride);
            }
        }
        if (trailStart < total) {
            trailMM[rec.Prefix] = std::vector<std::int32_t>(std::cbegin(rec.Skips) + trailStart,
                                                            std::cend(rec.Skips));
            if (hasMl) {
                trailML[rec.Prefix] = std::vector<std::int64_t>(
                    std::cbegin(mlValues) + qvOffset + trailStart * stride,
                    std::cbegin(mlValues) + qvOffset + total * stride);
            }
        }

        // prefix-lost-bases = canonical bases before the clip minus canonical bases
        // consumed by the front-removed mods. Added back to the first retained skip
        // on restore.
        const char canonical{std::empty(rec.Prefix) ? '\0' : rec.Prefix[0]};
        const std::int32_t basesBeforeClip{
            static_cast<std::int32_t>(std::ranges::count(sequence.substr(0, clipLeft), canonical))};
        std::int32_t prefixAtLastFront{0};
        for (std::size_t j = 0; j < frontN; ++j) {
            prefixAtLastFront += rec.Skips[j] + 1;
        }
        if (const std::int32_t lost{basesBeforeClip - prefixAtLastFront}; lost > 0) {
            pMM[rec.Prefix] = lost;
        }

        qvOffset += total * stride;
    }

    if (!leadMM.empty()) {
        lead["MM"] = std::move(leadMM);
    }
    if (!leadML.empty()) {
        lead["ML"] = std::move(leadML);
    }
    if (!trailMM.empty()) {
        trail["MM"] = std::move(trailMM);
    }
    if (!trailML.empty()) {
        trail["ML"] = std::move(trailML);
    }
}

// Apply one nesting level's MM/ML lead/trail (and pMM fixup) to the running
// per-mod-type skip/qual maps. Inverse of CaptureBasemods.
void RestoreBasemodLevel(const JSON::Json& node,
                         std::map<std::string, std::vector<std::int32_t>>& sepMM,
                         std::map<std::string, std::vector<std::int64_t>>& sepML)
{
    std::map<std::string, std::int32_t> prefixLost;
    if (node.contains("pMM")) {
        for (const auto& [mod, value] : node["pMM"].items()) {
            prefixLost[mod] = value.get<std::int32_t>();
        }
    }

    if (node.contains("lead")) {
        const JSON::Json& cut{node["lead"]};
        // Recover the first retained skip. If nothing was retained, the lead and
        // trail skips concatenate to the original list unchanged, so pMM is unused.
        for (const auto& [mod, lost] : prefixLost) {
            if (lost <= 0) {
                continue;
            }
            if (auto it{sepMM.find(mod)}; (it != std::end(sepMM)) && !std::empty(it->second)) {
                it->second[0] += lost;
            }
        }
        if (cut.contains("MM")) {
            for (const auto& [mod, value] : cut["MM"].items()) {
                const std::vector<std::int32_t> skips = value.get<std::vector<std::int32_t>>();
                std::vector<std::int32_t>& dst{sepMM[mod]};
                dst.insert(std::cbegin(dst), std::cbegin(skips), std::cend(skips));
            }
        }
        if (cut.contains("ML")) {
            for (const auto& [mod, value] : cut["ML"].items()) {
                const std::vector<std::int64_t> quals = value.get<std::vector<std::int64_t>>();
                std::vector<std::int64_t>& dst{sepML[mod]};
                dst.insert(std::cbegin(dst), std::cbegin(quals), std::cend(quals));
            }
        }
    }

    if (node.contains("trail")) {
        const JSON::Json& cut{node["trail"]};
        if (cut.contains("MM")) {
            for (const auto& [mod, value] : cut["MM"].items()) {
                const std::vector<std::int32_t> skips = value.get<std::vector<std::int32_t>>();
                std::vector<std::int32_t>& dst{sepMM[mod]};
                dst.insert(std::cend(dst), std::cbegin(skips), std::cend(skips));
            }
        }
        if (cut.contains("ML")) {
            for (const auto& [mod, value] : cut["ML"].items()) {
                const std::vector<std::int64_t> quals = value.get<std::vector<std::int64_t>>();
                std::vector<std::int64_t>& dst{sepML[mod]};
                dst.insert(std::cend(dst), std::cbegin(quals), std::cend(quals));
            }
        }
    }
}

}  // namespace

void ClipToQueryLossless(BamRecord& record, std::int32_t clipLeft, std::int32_t clipRight)
{
    const std::string sequence{record.Sequence()};
    const std::int32_t seqLen{static_cast<std::int32_t>(std::ssize(sequence))};
    if ((clipLeft < 0) || (clipRight > seqLen) || (clipRight <= clipLeft)) {
        return;
    }
    if ((clipLeft == 0) && (clipRight == seqLen)) {
        return;  // nothing trimmed
    }

    // Base modifications (MM/ML) are captured below. Other per-base tags that the
    // tag clipper trims (raw-subread pulse + subread-pileup) are not yet captured;
    // reject them rather than emit a restore blob that drops them on undo.
    for (const TagKey key : UNSUPPORTED_CLIP_TAGS) {
        if (record.Tags().Contains(key)) {
            throw std::runtime_error{
                "ClipToQueryLossless: lossless clipping of records with pulse or subread-pileup "
                "tags is not supported"};
        }
    }

    const std::string fastq{FastqFromQualities(record.Qualities())};
    const bool hasQual{!std::empty(fastq)};

    JSON::Json lead = JSON::Json::object();
    JSON::Json trail = JSON::Json::object();

    if (clipLeft > 0) {
        lead["sq"] = sequence.substr(0, clipLeft);
        if (hasQual) {
            lead["ql"] = fastq.substr(0, clipLeft);
        }
    }
    if (clipRight < seqLen) {
        trail["sq"] = sequence.substr(clipRight);
        if (hasQual) {
            trail["ql"] = fastq.substr(clipRight);
        }
    }

    const TagMap& tags{record.Tags()};
    for (const KineticTag& kt : KINETIC_TAGS) {
        const TagArray* array{GetArray(tags, kt.Key)};
        if (array == nullptr) {
            continue;
        }
        const std::vector<std::int64_t> values{DecodeIntArray(*array)};
        if (std::ssize(values) != seqLen) {
            continue;  // not a per-base track
        }
        if (clipLeft > 0) {
            // Forward tags: leading flank is the front; reverse tags are mirrored,
            // so their leading flank is the tail of the array.
            if (kt.Reverse) {
                lead[kt.Field] = std::vector<std::int64_t>(
                    std::cbegin(values) + (seqLen - clipLeft), std::cend(values));
            } else {
                lead[kt.Field] =
                    std::vector<std::int64_t>(std::cbegin(values), std::cbegin(values) + clipLeft);
            }
        }
        if (clipRight < seqLen) {
            if (kt.Reverse) {
                trail[kt.Field] = std::vector<std::int64_t>(
                    std::cbegin(values), std::cbegin(values) + (seqLen - clipRight));
            } else {
                trail[kt.Field] =
                    std::vector<std::int64_t>(std::cbegin(values) + clipRight, std::cend(values));
            }
        }
    }

    JSON::Json basemodPrefix = JSON::Json::object();
    CaptureBasemods(tags, sequence, clipLeft, clipRight, lead, trail, basemodPrefix);

    JSON::Json root = JSON::Json::object();
    root["lead"] = std::move(lead);
    root["trail"] = std::move(trail);
    if (!basemodPrefix.empty()) {
        root["pMM"] = std::move(basemodPrefix);
    }
    // Chain any prior lossless storage so repeated clips can all be undone.
    if (const TagArray* existing{GetArray(tags, LS_TAG)}; existing != nullptr) {
        root["nested"] = JSON::Json::from_msgpack(RawBytes(*existing));
    } else {
        root["nested"] = nullptr;
    }

    const std::vector<std::uint8_t> packed{JSON::Json::to_msgpack(root)};
    TagArray lsArray{'C'};
    for (const std::uint8_t byte : packed) {
        lsArray.AppendUInt8(byte);
    }
    record.MutableTags().Set(LS_TAG, lsArray);

    record.Clip(ClipType::CLIP_TO_QUERY, clipLeft, clipRight);
}

bool RestoreFromLossless(BamRecord& record)
{
    const TagArray* lsArray{GetArray(record.Tags(), LS_TAG)};
    if (lsArray == nullptr) {
        return false;
    }
    const std::vector<std::uint8_t> packed{RawBytes(*lsArray)};

    std::string sequence{record.Sequence()};
    std::string fastq{FastqFromQualities(record.Qualities())};
    const bool hasQual{!std::empty(fastq)};

    struct Track
    {
        char ElementType;
        bool Reverse;
        std::vector<std::int64_t> Values;
    };

    std::map<std::string, Track> tracks;  // keyed by `ls` field name
    for (const KineticTag& kt : KINETIC_TAGS) {
        const TagArray* array{GetArray(record.Tags(), kt.Key)};
        if (array == nullptr) {
            continue;
        }
        tracks.emplace(kt.Field, Track{array->ElementType(), kt.Reverse, DecodeIntArray(*array)});
    }

    // Base modifications: split the clipped record's MM/ML per modification type so
    // each nesting level can prepend/append its removed calls.
    std::map<std::string, std::vector<std::int32_t>> sepMM;
    std::map<std::string, std::vector<std::int64_t>> sepML;
    // Detect basemods from the clipped record's MM tag. The tag clipper keeps an
    // emptied MM as "prefix;" rather than dropping it, so this stays true whenever
    // the original had MM and the `ls` lead/trail calls can be reattached.
    bool hasBasemods{false};
    char mlElementType{'C'};
    if (const TagValue* mmValue{record.Tags().Get(TagKey{'M', 'M'})}; mmValue != nullptr) {
        if (const std::string* mm{std::get_if<std::string>(mmValue)}; mm != nullptr) {
            hasBasemods = true;
            const std::vector<BasemodRecord> records{ParseBasemodString(*mm)};
            std::vector<std::int64_t> ml;
            if (const TagArray* mlArray{GetArray(record.Tags(), TagKey{'M', 'L'})};
                mlArray != nullptr) {
                ml = DecodeIntArray(*mlArray);
                mlElementType = mlArray->ElementType();
            }
            std::size_t qvOffset{0};
            for (const BasemodRecord& rec : records) {
                sepMM[rec.Prefix] = rec.Skips;
                const std::size_t n{std::size(rec.Skips)};
                if (qvOffset + n <= std::size(ml)) {
                    sepML[rec.Prefix] = std::vector<std::int64_t>(std::cbegin(ml) + qvOffset,
                                                                  std::cbegin(ml) + qvOffset + n);
                }
                qvOffset += n;
            }
        }
    }

    const auto applyFlank = [&](const JSON::Json& cut, bool isLead) {
        if (cut.contains("sq")) {
            const std::string piece = cut["sq"];
            sequence = isLead ? (piece + sequence) : (sequence + piece);
        }
        if (hasQual && cut.contains("ql")) {
            const std::string piece = cut["ql"];
            fastq = isLead ? (piece + fastq) : (fastq + piece);
        }
        for (auto& [field, track] : tracks) {
            if (!cut.contains(field)) {
                continue;
            }
            const std::vector<std::int64_t> piece = cut[field];
            // Forward tags prepend the leading flank; reverse tags are mirrored.
            const bool prepend{isLead != track.Reverse};
            if (prepend) {
                track.Values.insert(std::cbegin(track.Values), std::cbegin(piece),
                                    std::cend(piece));
            } else {
                track.Values.insert(std::cend(track.Values), std::cbegin(piece), std::cend(piece));
            }
        }
    };

    // NOTE: copy-init, not brace-init: `Json node{parsedObject}` would invoke the
    // initializer_list constructor and wrap the object in a one-element array.
    JSON::Json node = JSON::Json::from_msgpack(packed);
    while (!node.is_null()) {
        if (node.contains("lead")) {
            applyFlank(node["lead"], true);
        }
        if (node.contains("trail")) {
            applyFlank(node["trail"], false);
        }
        if (hasBasemods) {
            RestoreBasemodLevel(node, sepMM, sepML);
        }
        node = node.contains("nested") ? node["nested"] : JSON::Json{};
    }

    record.Sequence(sequence);
    if (hasQual) {
        record.Qualities(QualitiesFromFastq(fastq));
    }
    for (const KineticTag& kt : KINETIC_TAGS) {
        if (const auto it{tracks.find(kt.Field)}; it != std::cend(tracks)) {
            record.MutableTags().Set(kt.Key,
                                     EncodeIntArray(it->second.ElementType, it->second.Values));
        }
    }
    if (hasBasemods) {
        // Rebuild MM/ML in sorted modification-type order (ML follows MM order).
        std::vector<BasemodRecord> outRecords;
        std::vector<std::int64_t> outMl;
        for (const auto& [mod, skips] : sepMM) {
            outRecords.push_back(BasemodRecord{mod, skips});
            if (const auto it{sepML.find(mod)}; it != std::cend(sepML)) {
                outMl.insert(std::cend(outMl), std::cbegin(it->second), std::cend(it->second));
            }
        }
        record.MutableTags().Set(TagKey{'M', 'M'}, WriteBasemodString(outRecords));
        if (!std::empty(outMl)) {
            record.MutableTags().Set(TagKey{'M', 'L'}, EncodeIntArray(mlElementType, outMl));
        }
    }
    for (const TagKey key : DEMUX_TAGS) {
        record.MutableTags().Remove(key);
    }
    return true;
}

}  // namespace Samoa
}  // namespace PacBio
