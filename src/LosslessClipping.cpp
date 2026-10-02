#include <pbsamoa/core/LosslessClipping.hpp>

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/Basemods.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <pbcopper/data/LocalContextFlags.h>
#include <pbcopper/json/JSON.h>

#include <algorithm>
#include <array>
#include <iterator>
#include <map>
#include <optional>
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
constexpr TagKey BC_TAG{'b', 'c'};
constexpr TagKey BQ_TAG{'b', 'q'};

struct KineticTag
{
    TagKey Key;
    const char* Field;  ///< `ls` JSON field name (== tag name)
    bool Reverse;       ///< stored mirrored relative to the sequence (rp/ri)
};

// Per-base tracks: kinetics plus the subread-pileup match/mismatch counts (sm/sx).
constexpr std::array<KineticTag, 8> KINETIC_TAGS{{
    {TagKey{'p', 'w'}, "pw", false},
    {TagKey{'i', 'p'}, "ip", false},
    {TagKey{'f', 'p'}, "fp", false},
    {TagKey{'f', 'i'}, "fi", false},
    {TagKey{'r', 'p'}, "rp", true},
    {TagKey{'r', 'i'}, "ri", true},
    {TagKey{'s', 'm'}, "sm", false},
    {TagKey{'s', 'x'}, "sx", false},
}};

// Subread-pileup coverage: run-length encoded [runLen, cov, ...] pairs summing to seqLen.
constexpr TagKey SA_TAG{'s', 'a'};
constexpr TagKey CX_TAG{'c', 'x'};
constexpr TagKey QS_TAG{'q', 's'};
constexpr TagKey QE_TAG{'q', 'e'};
constexpr std::int64_t CX_BARCODE_BITS{Data::LocalContextFlags::BARCODE_BEFORE |
                                       Data::LocalContextFlags::BARCODE_AFTER};

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
    std::vector<std::uint8_t> out(std::size(data));
    std::memcpy(std::data(out), std::data(data), std::size(data));
    return out;
}

std::string FastqFromQualities(std::span<const std::uint8_t> quals)
{
    std::string fastq(std::size(quals), '\0');
    std::ranges::transform(quals, std::begin(fastq),
                           [](std::uint8_t qv) { return static_cast<char>(qv + 33); });
    return fastq;
}

std::vector<std::uint8_t> QualitiesFromFastq(std::string_view fastq)
{
    std::vector<std::uint8_t> quals(std::size(fastq));
    std::ranges::transform(fastq, std::begin(quals),
                           [](char c) { return static_cast<std::uint8_t>(c - 33); });
    return quals;
}

// Per-base tags that the PacBio default clipper trims but lossless storage does not
// capture: legacy per-base QV/tag strings (dq/iq/mq/sq/dt/st), raw-subread pulse, and
// raw-subread pulse frames. Only the per-base tracks, sa, and MM/ML are captured;
// everything else is rejected up front so a clip can never silently drop data that
// undo would need.
constexpr std::array<TagKey, 19> UNSUPPORTED_CLIP_TAGS{{
    TagKey{'d', 'q'}, TagKey{'i', 'q'}, TagKey{'m', 'q'}, TagKey{'s', 'q'}, TagKey{'d', 't'},
    TagKey{'s', 't'}, TagKey{'p', 'c'}, TagKey{'p', 't'}, TagKey{'p', 'q'}, TagKey{'p', 'v'},
    TagKey{'p', 'g'}, TagKey{'p', 'a'}, TagKey{'p', 'm'}, TagKey{'p', 's'}, TagKey{'p', 'i'},
    TagKey{'p', 'd'}, TagKey{'p', 'x'}, TagKey{'p', 'e'}, TagKey{'s', 'f'},
}};

// Split the run-length encoded `sa` coverage at the clip window, lima-style: runs
// fully inside a flank move to lead/trail; a run straddling a boundary stays in the
// retained array with its length reduced, and the removed base count is recorded as
// psa/ssa so undo can grow it back.
struct PileupSplit
{
    std::vector<std::int64_t> Lead;
    std::vector<std::int64_t> Retained;
    std::vector<std::int64_t> Trail;
    std::int64_t LostPrefix{0};
    std::int64_t LostSuffix{0};
};

PileupSplit SplitPileupCoverage(const std::vector<std::int64_t>& sa, std::int64_t seqLen,
                                std::int64_t clipLeft, std::int64_t clipRight)
{
    PileupSplit split;
    std::int64_t runStart{0};
    for (std::size_t i = 0; i + 1 < std::size(sa); i += 2) {
        const std::int64_t runLen{sa[i]};
        const std::int64_t runEnd{runStart + runLen};
        std::vector<std::int64_t>* dst{&split.Retained};
        if (runEnd <= clipLeft) {
            dst = &split.Lead;
        } else if (runStart >= clipRight) {
            dst = &split.Trail;
        }
        dst->push_back(runLen);
        dst->push_back(sa[i + 1]);
        runStart = runEnd;
    }
    if (runStart != seqLen) {
        throw std::runtime_error{"ClipToQueryLossless: sa run lengths do not sum to SEQ length"};
    }
    std::int64_t leadBases{0};
    for (std::size_t i = 0; i < std::size(split.Lead); i += 2) {
        leadBases += split.Lead[i];
    }
    std::int64_t trailBases{0};
    for (std::size_t i = 0; i < std::size(split.Trail); i += 2) {
        trailBases += split.Trail[i];
    }
    split.LostPrefix = clipLeft - leadBases;
    split.LostSuffix = (seqLen - clipRight) - trailBases;
    split.Retained[0] -= split.LostPrefix;
    split.Retained[std::size(split.Retained) - 2] -= split.LostSuffix;
    return split;
}

// Capture the removed base-modification calls (MM/ML) into the `ls` lead/trail
// maps, plus the per-mod prefix-lost-base count (pMM) needed to recover the first
// retained skip on undo. Mirrors lima's format. MM/ML are per-mod-type maps keyed
// by the modification prefix (e.g. "C+m?"); ML quals follow MM-record order.
// Lossless clipping targets CCS/demultiplexed reads, which are unmapped (origQueryStart
// is 0 and the strand flag is unset), so clipLeft is the SEQ-space offset and basemods are
// walked in stored orientation. Reverse-strand handling lives in the lossy clip path
// (ClipBasemodRecord's reverse flag), where the aligned SEQ is reverse-complemented.
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
    // If ML is present it must cover every MM site. A short/mismatched ML cannot be
    // sliced losslessly (undo would lose probabilities), so reject the clip rather than
    // silently dropping data — consistent with the UNSUPPORTED_CLIP_TAGS contract.
    if (!std::empty(mlValues) && (std::size(mlValues) < totalMlValues)) {
        throw std::runtime_error{
            "ClipToQueryLossless: ML array shorter than the MM modification sites it must cover"};
    }
    const bool hasMl{!std::empty(mlValues)};
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

        // Like lima, every modification type is written into both flanks, as an empty list when
        // the flank holds none of its sites. lima-undo depends on it: it rewrites MM with the
        // types sorted by name and reorders ML to match only when a flank carries an ML map.
        // Without the maps, a record whose MM is not already sorted comes back from lima-undo
        // with its ML probabilities assigned to the wrong modification types.
        leadMM[rec.Prefix] =
            std::vector<std::int32_t>(std::cbegin(rec.Skips), std::cbegin(rec.Skips) + frontN);
        trailMM[rec.Prefix] =
            std::vector<std::int32_t>(std::cbegin(rec.Skips) + trailStart, std::cend(rec.Skips));
        if (hasMl) {
            leadML[rec.Prefix] =
                std::vector<std::int64_t>(std::cbegin(mlValues) + qvOffset,
                                          std::cbegin(mlValues) + qvOffset + frontN * stride);
            trailML[rec.Prefix] =
                std::vector<std::int64_t>(std::cbegin(mlValues) + qvOffset + trailStart * stride,
                                          std::cbegin(mlValues) + qvOffset + total * stride);
        }

        // prefix-lost-bases = canonical bases before the clip not consumed by the
        // front-removed mods, added back to the first retained skip on restore.
        // ClipBasemodRecord already computes this with the correct wildcard-aware count
        // (canonical 'N' matches every base, not the literal character 'N').
        // With no retained site there is no skip to add it to: the lead and trail lists
        // concatenate to the original. Write 0 then, as lima-undo takes a positive count for
        // a type with nothing retained as a cue to rewrite the trail skips.
        pMM[rec.Prefix] = (window.Retained > 0) ? window.PrefixLost : 0;

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

    // The default tag clipper only handles a uint8 `sa`; clip it here and reinstall the
    // retained runs after the query clip.
    std::optional<TagArray> retainedSa;
    std::int64_t lostPrefixSa{0};
    std::int64_t lostSuffixSa{0};
    if (const TagArray* saArray{GetArray(tags, SA_TAG)}; saArray != nullptr) {
        const PileupSplit split{
            SplitPileupCoverage(DecodeIntArray(*saArray), seqLen, clipLeft, clipRight)};
        if (clipLeft > 0) {
            lead["sa"] = split.Lead;
        }
        if (clipRight < seqLen) {
            trail["sa"] = split.Trail;
        }
        retainedSa = EncodeIntArray(saArray->ElementType(), split.Retained);
        lostPrefixSa = split.LostPrefix;
        lostSuffixSa = split.LostSuffix;
    }

    JSON::Json root = JSON::Json::object();
    root["lead"] = std::move(lead);
    root["trail"] = std::move(trail);
    if (!basemodPrefix.empty()) {
        root["pMM"] = std::move(basemodPrefix);
    }
    // lima writes the pileup fixups for every record, 0 when it has no `sa`.
    root["psa"] = lostPrefixSa;
    root["ssa"] = lostSuffixSa;
    root["csa"] = -1;  // lima: coverage of an emptied window; impossible here
    // lima also notes the barcode call that led to the clip.
    if (const TagArray* call{GetArray(tags, BC_TAG)}; (call != nullptr) && (call->Count() == 2)) {
        const std::vector<std::int64_t> pair{DecodeIntArray(*call)};
        root["bc"] = std::to_string(pair.front()) + "--" + std::to_string(pair.back());
    }
    if (const TagValue* quality{tags.Get(BQ_TAG)}; quality != nullptr) {
        if (const std::int64_t* value{std::get_if<std::int64_t>(quality)}; value != nullptr) {
            root["bq"] = *value;
        }
    }
    // Chain any prior lossless storage so repeated clips can all be undone.
    if (const TagArray* existing{GetArray(tags, LS_TAG)}; existing != nullptr) {
        root["nested"] = JSON::Json::from_msgpack(RawBytes(*existing));
    } else {
        root["nested"] = nullptr;
    }

    const std::vector<std::uint8_t> packed{JSON::Json::to_msgpack(root)};
    TagArray lsArray{'C'};
    lsArray.Resize(static_cast<std::uint32_t>(std::size(packed)));
    std::memcpy(std::data(lsArray.MutableData()), std::data(packed), std::size(packed));
    record.MutableTags().Set(LS_TAG, lsArray);

    record.MutableTags().Remove(SA_TAG);
    record.Clip(ClipType::CLIP_TO_QUERY, clipLeft, clipRight);
    if (retainedSa) {
        record.MutableTags().Set(SA_TAG, std::move(*retainedSa));
    }
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
                // ML carries ModCodeCount values per site (interleaved); multi-code mods
                // (e.g. "C+mh") therefore occupy stride entries per skip. Seed by the same
                // stride capture used so retained multi-code sites keep every probability.
                const std::size_t n{std::size(rec.Skips) * ModCodeCount(rec.Prefix)};
                if (qvOffset + n <= std::size(ml)) {
                    sepML[rec.Prefix] = std::vector<std::int64_t>(std::cbegin(ml) + qvOffset,
                                                                  std::cbegin(ml) + qvOffset + n);
                }
                qvOffset += n;
            }
        }
    }

    std::optional<Track> saTrack;
    if (const TagArray* saArray{GetArray(record.Tags(), SA_TAG)}; saArray != nullptr) {
        saTrack = Track{saArray->ElementType(), false, DecodeIntArray(*saArray)};
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
        if (saTrack && cut.contains("sa")) {
            const std::vector<std::int64_t> piece = cut["sa"];
            std::vector<std::int64_t>& sa{saTrack->Values};
            sa.insert(isLead ? std::cbegin(sa) : std::cend(sa), std::cbegin(piece),
                      std::cend(piece));
        }
    };

    // Grow the boundary runs back before re-attaching the flank runs (lima-undo order).
    const auto restorePileupRuns = [&](const JSON::Json& node) {
        if (!saTrack) {
            return;
        }
        std::vector<std::int64_t>& sa{saTrack->Values};
        if (node.contains("csa") && (node["csa"].get<std::int64_t>() != -1) && std::empty(sa)) {
            sa = {0, node["csa"].get<std::int64_t>()};
        }
        if (node.contains("psa") && !std::empty(sa)) {
            sa[0] += node["psa"].get<std::int64_t>();
        }
        if (node.contains("ssa") && !std::empty(sa)) {
            sa[std::size(sa) - 2] += node["ssa"].get<std::int64_t>();
        }
    };

    // NOTE: copy-init, not brace-init: `Json node{parsedObject}` would invoke the
    // initializer_list constructor and wrap the object in a one-element array.
    JSON::Json node = JSON::Json::from_msgpack(packed);
    while (!node.is_null()) {
        restorePileupRuns(node);
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
        // MN tracks SEQ length for MM/ML; the clipped record carries the clipped length,
        // so restore it to the rebuilt full-length sequence.
        if (record.Tags().Contains(TagKey{'M', 'N'})) {
            record.MutableTags().Set(TagKey{'M', 'N'}, std::ssize(sequence));
        }
    }
    if (saTrack) {
        record.MutableTags().Set(SA_TAG, EncodeIntArray(saTrack->ElementType, saTrack->Values));
    }
    for (const TagKey key : DEMUX_TAGS) {
        record.MutableTags().Remove(key);
    }

    // Demultiplexing side effects beyond the clip itself, mirroring lima-undo.
    TagMap& tags{record.MutableTags()};
    // CCS names end in "/ccs[/fwd|/rev]"; segment reads append "/qs_qe" and keep qs/qe.
    std::string_view name{record.Name()};
    if (name.ends_with("/fwd") || name.ends_with("/rev")) {
        name.remove_suffix(4);
    }
    const bool isCcs{name.ends_with("/ccs")};
    if (isCcs) {
        tags.Remove(QS_TAG);
        tags.Remove(QE_TAG);
    }
    if (const TagValue* cx{tags.Get(CX_TAG)}; cx != nullptr) {
        if (const std::int64_t* value{std::get_if<std::int64_t>(cx)}; value != nullptr) {
            const std::int64_t cleared{*value & ~CX_BARCODE_BITS};
            if ((cleared == 0) && isCcs) {
                tags.Remove(CX_TAG);
            } else {
                tags.Set(CX_TAG, cleared);
            }
        }
    }
    if (const TagValue* rg{tags.Get(RG_TAG)}; rg != nullptr) {
        if (const std::string* id{std::get_if<std::string>(rg)}; id != nullptr) {
            tags.Set(RG_TAG, std::string{ReadGroupBaseId(*id)});
        }
    }
    return true;
}

}  // namespace Samoa
}  // namespace PacBio
