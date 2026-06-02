#include <pbsamoa/core/BamRecord.hpp>

#include "BinaryUtils.hpp"

#include <pbsamoa/core/CigarClipping.hpp>
#include <pbsamoa/core/Sequence.hpp>
#include <pbsamoa/core/TagClipping.hpp>

#include <algorithm>
#include <bit>
#include <format>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <utility>
#include <variant>
#include <vector>

#include <cstring>

namespace PacBio {
namespace Samoa {

namespace {

constexpr TagKey QS_TAG{'q', 's'};
constexpr TagKey QE_TAG{'q', 'e'};
constexpr std::uint16_t BAM_UNMAPPED_FLAG{0x4};
constexpr std::uint16_t BAM_REVERSE_STRAND_FLAG{0x10};
constexpr std::uint16_t BAM_NON_PRIMARY_FLAGS{0x900};
constexpr std::uint16_t BAM_UNMAPPED_BIN{4680};

void WriteU16LEAt(std::byte* dst, std::uint16_t v)
{
    dst[0] = static_cast<std::byte>(v & 0xFFU);
    dst[1] = static_cast<std::byte>((v >> 8U) & 0xFFU);
}

void WriteU32LEAt(std::byte* dst, std::uint32_t v)
{
    dst[0] = static_cast<std::byte>(v & 0xFFU);
    dst[1] = static_cast<std::byte>((v >> 8U) & 0xFFU);
    dst[2] = static_cast<std::byte>((v >> 16U) & 0xFFU);
    dst[3] = static_cast<std::byte>((v >> 24U) & 0xFFU);
}

void WriteI32LEAt(std::byte* dst, std::int32_t v)
{
    WriteU32LEAt(dst, std::bit_cast<std::uint32_t>(v));
}

std::int32_t IntTagOr(const TagMap& tags, TagKey key, std::int32_t fallback)
{
    if (const auto* tagValue{tags.Get(key)}) {
        if (const auto* integral{std::get_if<std::int64_t>(tagValue)}) {
            return static_cast<std::int32_t>(*integral);
        }
        if (const auto* character{std::get_if<char>(tagValue)}) {
            return static_cast<std::int32_t>(*character);
        }
        throw std::runtime_error{
            std::format("Tag '{}{}' has unexpected type", key.First(), key.Second())};
    }
    return fallback;
}

std::uint16_t ComputeBamBin(std::int32_t pos, CigarView cigar)
{
    if (pos < 0) {
        return BAM_UNMAPPED_BIN;
    }
    if (std::empty(cigar)) {
        return Reg2Bin(pos, pos + 1);
    }

    const std::int32_t referenceEnd{static_cast<std::int32_t>(pos + ReferenceLength(cigar))};
    return Reg2Bin(pos, NonEmptyAlignmentEnd(pos, referenceEnd));
}

}  // namespace

// --- accessors ---

std::string_view BamRecord::Name() const { return name_; }

std::uint16_t BamRecord::Flag() const { return flag_; }

std::int32_t BamRecord::RefId() const { return refId_; }

std::int32_t BamRecord::Pos() const { return pos_; }

std::uint8_t BamRecord::MapQ() const { return mapQ_; }

CigarView BamRecord::Cigar() const { return cigar_; }

std::int32_t BamRecord::NextRefId() const { return nextRefId_; }

std::int32_t BamRecord::NextPos() const { return nextPos_; }

std::int32_t BamRecord::Tlen() const { return tlen_; }

std::string_view BamRecord::Sequence() const { return sequence_; }

std::span<const std::uint8_t> BamRecord::Qualities() const { return qualities_; }

const TagMap& BamRecord::Tags() const { return tags_; }

// --- mutators ---

BamRecord& BamRecord::Name(std::string name)
{
    // BAM spec: l_read_name is uint8_t, includes NUL terminator => max 254 chars
    if (std::size(name) > 254) {
        throw std::invalid_argument{
            "BamRecord::Name: read name exceeds BAM limit of 254 characters"};
    }
    name_ = std::move(name);
    return *this;
}

BamRecord& BamRecord::Flag(std::uint16_t flag)
{
    flag_ = flag;
    return *this;
}

BamRecord& BamRecord::RefId(std::int32_t refId)
{
    refId_ = refId;
    return *this;
}

BamRecord& BamRecord::Pos(std::int32_t pos)
{
    pos_ = pos;
    return *this;
}

BamRecord& BamRecord::MapQ(std::uint8_t mapq)
{
    mapQ_ = mapq;
    return *this;
}

BamRecord& BamRecord::Cigar(std::vector<CigarOp> cigar)
{
    cigar_ = std::move(cigar);
    return *this;
}

BamRecord& BamRecord::NextRefId(std::int32_t nextRefId)
{
    nextRefId_ = nextRefId;
    return *this;
}

BamRecord& BamRecord::NextPos(std::int32_t nextPos)
{
    nextPos_ = nextPos;
    return *this;
}

BamRecord& BamRecord::Tlen(std::int32_t tlen)
{
    tlen_ = tlen;
    return *this;
}

BamRecord& BamRecord::Sequence(std::string seq)
{
    sequence_ = std::move(seq);
    return *this;
}

BamRecord& BamRecord::Qualities(std::vector<std::uint8_t> qual)
{
    qualities_ = std::move(qual);
    return *this;
}

BamRecord& BamRecord::Tags(TagMap tags)
{
    tags_ = std::move(tags);
    return *this;
}

// --- derived fields ---

bool BamRecord::IsMapped() const { return (flag_ & BAM_UNMAPPED_FLAG) == 0; }

bool BamRecord::IsReverseStrand() const { return (flag_ & BAM_REVERSE_STRAND_FLAG) != 0; }

bool BamRecord::IsPrimary() const { return (flag_ & BAM_NON_PRIMARY_FLAGS) == 0; }

std::int32_t BamRecord::ReferenceEnd() const { return pos_ + ReferenceLength(cigar_); }

// --- serialization ---

std::vector<std::byte> BamRecord::SerializeToBam() const
{
    // Compute variable-length field sizes
    if (std::size(name_) > 254) {
        throw std::invalid_argument{
            "BamRecord::SerializeToBam: read name exceeds "
            "BAM limit of 254 characters"};
    }
    const std::uint8_t nameLen{static_cast<std::uint8_t>(std::size(name_) + 1)};  // +NUL

    // A CIGAR with >65535 ops cannot fit n_cigar_op; BAM stores a 2-op 'kSmN' placeholder
    // and the real ops in a CG:B,I tag (SAMv1 §4.2.2), matching htslib.
    const std::size_t realOpCount{std::size(cigar_)};
    const bool longCigar{realOpCount > std::numeric_limits<std::uint16_t>::max()};
    const std::int64_t cigRefLen{ReferenceLength(cigar_)};
    if (longCigar && (cigRefLen >= (std::int64_t{1} << 28))) {
        // The placeholder's N op length is the reference span, which must fit the 28-bit
        // CIGAR field; htslib falls back to SAM/CRAM when it does not.
        throw std::invalid_argument{
            "BamRecord::SerializeToBam: CIGAR reference span too large for the CG-tag "
            "long-CIGAR encoding; write SAM or CRAM instead"};
    }
    const std::uint16_t nCigarOp{longCigar ? std::uint16_t{2}
                                           : static_cast<std::uint16_t>(realOpCount)};
    const std::uint32_t seqLen{static_cast<std::uint32_t>(std::size(sequence_))};
    if (!std::empty(qualities_) && (std::size(qualities_) != seqLen)) {
        throw std::invalid_argument{
            std::format("BamRecord::SerializeToBam: quality length ({}) != "
                        "sequence length ({})",
                        std::size(qualities_), seqLen)};
    }
    const std::uint32_t packedSeqLen{(seqLen + 1) / 2};

    // Compute bin
    const std::uint16_t bin{ComputeBamBin(pos_, cigar_)};

    // Compute tag size without allocating
    const std::size_t tagSize{SerializedBamSize(tags_)};

    // CG spillover: 'C','G','B','I' (4) + uint32 count (4) + 4 bytes per real op.
    const std::size_t cgTagBytes{longCigar ? (8U + 4U * realOpCount) : 0U};

    // Total size: 32 (fixed) + nameLen + 4*nCigarOp + packedSeqLen + seqLen + tagSize + CG tag
    const std::size_t totalSize{32U + nameLen + 4U * nCigarOp + packedSeqLen + seqLen + tagSize +
                                cgTagBytes};

    std::vector<std::byte> result(totalSize);
    std::byte* p{std::data(result)};

    // Fixed fields (32 bytes)
    WriteI32LEAt(p + 0, refId_);
    WriteI32LEAt(p + 4, pos_);
    p[8] = static_cast<std::byte>(nameLen);
    p[9] = static_cast<std::byte>(mapQ_);
    WriteU16LEAt(p + 10, bin);
    WriteU16LEAt(p + 12, nCigarOp);
    WriteU16LEAt(p + 14, flag_);
    WriteU32LEAt(p + 16, seqLen);
    WriteI32LEAt(p + 20, nextRefId_);
    WriteI32LEAt(p + 24, nextPos_);
    WriteI32LEAt(p + 28, tlen_);

    std::size_t offset{32};

    // Read name (NUL-terminated)
    std::ranges::copy_n(reinterpret_cast<const std::byte*>(std::data(name_)), std::size(name_),
                        p + offset);
    p[offset + std::size(name_)] = std::byte{0};
    offset += nameLen;

    // CIGAR (array of uint32). A long CIGAR writes the 2-op 'kSmN' placeholder here; the
    // real ops follow in a CG tag appended after the regular tags.
    if (longCigar) {
        WriteU32LEAt(p + offset, (seqLen << 4U) | std::to_underlying(CigarOpType::S));
        offset += 4;
        WriteU32LEAt(p + offset, (static_cast<std::uint32_t>(cigRefLen) << 4U) |
                                     std::to_underlying(CigarOpType::N));
        offset += 4;
    } else {
        for (const CigarOp& op : cigar_) {
            WriteU32LEAt(p + offset, op.RawValue());
            offset += 4;
        }
    }

    // Packed sequence (write directly into result buffer)
    if (seqLen > 0) {
        PackSequenceInto(sequence_, p + offset);
    }
    offset += packedSeqLen;

    // Quality scores
    if (!std::empty(qualities_)) {
        std::ranges::copy_n(reinterpret_cast<const std::byte*>(std::data(qualities_)), seqLen,
                            p + offset);
    } else if (seqLen > 0) {
        // Missing quality: fill with 0xFF
        std::ranges::fill_n(p + offset, seqLen, std::byte{0xFF});
    }
    offset += seqLen;

    // Tags (write directly into result buffer)
    if (tagSize > 0) {
        AppendTagsToBam(tags_, p + offset);
    }
    offset += tagSize;

    // CG tag carrying the real long CIGAR as a B:I array of packed ops.
    if (longCigar) {
        p[offset + 0] = std::byte{'C'};
        p[offset + 1] = std::byte{'G'};
        p[offset + 2] = std::byte{'B'};
        p[offset + 3] = std::byte{'I'};
        WriteU32LEAt(p + offset + 4, static_cast<std::uint32_t>(realOpCount));
        offset += 8;
        for (const CigarOp& op : cigar_) {
            WriteU32LEAt(p + offset, op.RawValue());
            offset += 4;
        }
    }

    return result;
}

// --- clipping ---

BamRecord& BamRecord::Clip(ClipType type, std::int32_t start, std::int32_t end,
                           bool exciseFlankingInserts)
{
    static const TagClipper DEFAULT_CLIPPER{TagClipper::PacBioDefault()};
    return Clip(type, start, end, DEFAULT_CLIPPER, exciseFlankingInserts);
}

BamRecord& BamRecord::Clip(ClipType type, std::int32_t start, std::int32_t end,
                           const TagClipper& clipper, bool exciseFlankingInserts)
{
    // Unmapped + reference clip => no-op
    if (type == ClipType::CLIP_TO_REFERENCE && !IsMapped()) {
        return *this;
    }

    const bool isReverse{IsReverseStrand()};
    const std::size_t origSeqLen{std::size(sequence_)};

    // Get query start/end from tags (default to 0/seqLen for CCS/transcript)
    const std::int32_t origQStart{IntTagOr(tags_, QS_TAG, 0)};
    const std::int32_t origQEnd{IntTagOr(tags_, QE_TAG, static_cast<std::int32_t>(origSeqLen))};

    // Compute ClipResult
    ClipResult clipResult;
    if (type == ClipType::CLIP_TO_QUERY) {
        clipResult = ClipCigarToQuery(cigar_, start, end, origQStart, origQEnd, pos_, isReverse);
    } else {
        clipResult =
            ClipCigarToReference(cigar_, start, end, pos_, isReverse, exciseFlankingInserts);
    }

    // Clip tags BEFORE modifying sequence (basemods reads old sequence)
    clipper.ClipTags(tags_, clipResult.clipOffset, clipResult.clipLength, origSeqLen, sequence_);

    // Clip sequence and qualities (in-place to avoid extra allocation)
    sequence_.erase(0, clipResult.clipOffset);
    sequence_.resize(clipResult.clipLength);
    if (!std::empty(qualities_)) {
        qualities_.erase(
            std::begin(qualities_),
            std::begin(qualities_) + static_cast<std::ptrdiff_t>(clipResult.clipOffset));
        qualities_.resize(clipResult.clipLength);
    }

    // Update CIGAR and position
    cigar_ = std::move(clipResult.cigar);
    pos_ = clipResult.newPos;

    // Update qs/qe tags
    if (!tags_.Contains(QS_TAG)) {
        return *this;
    }

    if (type == ClipType::CLIP_TO_QUERY) {
        tags_.Set(QS_TAG, std::int64_t{start});
        tags_.Set(QE_TAG, std::int64_t{end});
    } else {
        // For reference clipping, qs/qe shift by clipOffset
        const std::int32_t newQs{origQStart + static_cast<std::int32_t>(clipResult.clipOffset)};
        tags_.Set(QS_TAG, std::int64_t{newQs});
        tags_.Set(QE_TAG, std::int64_t{newQs + static_cast<std::int32_t>(clipResult.clipLength)});
    }

    return *this;
}

BamRecord BamRecord::Clipped(ClipType type, std::int32_t start, std::int32_t end,
                             bool exciseFlankingInserts) const
{
    BamRecord copy{*this};
    copy.Clip(type, start, end, exciseFlankingInserts);
    return copy;
}

BamRecord BamRecord::Clipped(ClipType type, std::int32_t start, std::int32_t end,
                             const TagClipper& clipper, bool exciseFlankingInserts) const
{
    BamRecord copy{*this};
    copy.Clip(type, start, end, clipper, exciseFlankingInserts);
    return copy;
}

TagMap& BamRecord::MutableTags() { return tags_; }

// --- PacBio BAM accessors ---

namespace {

constexpr TagKey CX_TAG{'c', 'x'};
constexpr TagKey SN_TAG{'s', 'n'};
constexpr TagKey RQ_TAG{'r', 'q'};
constexpr TagKey PW_TAG{'p', 'w'};
constexpr TagKey IP_TAG{'i', 'p'};
constexpr TagKey WS_TAG{'w', 's'};
constexpr TagKey WE_TAG{'w', 'e'};

std::int32_t TagToInt32(const TagValue& value)
{
    if (const auto* integral{std::get_if<std::int64_t>(&value)}; integral) {
        return static_cast<std::int32_t>(*integral);
    }
    if (const auto* character{std::get_if<char>(&value)}; character) {
        return static_cast<std::int32_t>(*character);
    }
    throw std::runtime_error{"Expected integral PacBio BAM tag"};
}

float TagToFloat(const TagValue& value)
{
    if (const auto* floating{std::get_if<float>(&value)}; floating) {
        return *floating;
    }
    throw std::runtime_error{"Expected floating-point PacBio BAM tag"};
}

const TagValue& RequiredTag(const TagMap& tags, TagKey key, std::string_view missingMessage)
{
    const auto* tag{tags.Get(key)};
    if (!tag) {
        throw std::runtime_error{std::string{missingMessage}};
    }
    return *tag;
}

const TagArray& RequiredTagArray(const TagValue& value, std::string_view malformedMessage)
{
    const auto* array{std::get_if<TagArray>(&value)};
    if (!array) {
        throw std::runtime_error{std::string{malformedMessage}};
    }
    return *array;
}

std::optional<std::int32_t> OptionalIntTag(const TagMap& tags, TagKey key)
{
    if (const auto* tag{tags.Get(key)}; tag) {
        return TagToInt32(*tag);
    }
    return std::nullopt;
}

std::string_view QueryIntervalText(std::string_view fullName)
{
    const std::size_t lastSlash{fullName.rfind('/')};
    if (lastSlash == std::string_view::npos) {
        throw std::runtime_error{"Malformed PacBio BAM read name: " + std::string{fullName}};
    }
    return fullName.substr(lastSlash + 1);
}

std::pair<std::int32_t, std::int32_t> ParseQueryInterval(std::string_view fullName)
{
    const std::string_view interval{QueryIntervalText(fullName)};
    const std::size_t underscore{interval.find('_')};
    if (underscore == std::string_view::npos) {
        throw std::runtime_error{"Malformed PacBio BAM query interval: " + std::string{interval}};
    }
    const std::string_view startText{interval.substr(0, underscore)};
    const std::string_view endText{interval.substr(underscore + 1)};
    return {
        std::stoi(std::string{startText}),
        std::stoi(std::string{endText}),
    };
}

std::vector<std::uint8_t> ToUInt8Vector(const TagArray& array)
{
    std::vector<std::uint8_t> result;
    result.reserve(array.Count());
    const auto data{array.Data()};
    switch (array.ElementType()) {
        case 'C':
            for (const std::byte value : data) {
                result.push_back(static_cast<std::uint8_t>(value));
            }
            return result;
        case 'c':
            for (const std::byte value : data) {
                result.push_back(static_cast<std::uint8_t>(static_cast<std::int8_t>(value)));
            }
            return result;
        default:
            throw std::runtime_error{"Expected uint8/int8 BAM tag array"};
    }
}

template <typename T>
std::vector<T> DecodePodArray(const TagArray& array)
{
    std::vector<T> result(array.Count());
    std::memcpy(result.data(), array.Data().data(), std::size(result) * sizeof(T));
    return result;
}

Data::Frames FramesFromTagArray(const TagArray& array)
{
    switch (array.ElementType()) {
        case 'C':
        case 'c':
            return Data::Frames::Decode(ToUInt8Vector(array));
        case 'S':
            return Data::Frames{DecodePodArray<std::uint16_t>(array)};
        case 's': {
            const auto values{DecodePodArray<std::int16_t>(array)};
            return Data::Frames{std::vector<std::uint16_t>(values.begin(), values.end())};
        }
        case 'I': {
            const auto values{DecodePodArray<std::uint32_t>(array)};
            return Data::Frames{std::vector<std::uint16_t>(values.begin(), values.end())};
        }
        case 'i': {
            const auto values{DecodePodArray<std::int32_t>(array)};
            return Data::Frames{std::vector<std::uint16_t>(values.begin(), values.end())};
        }
        default:
            throw std::runtime_error{"Unsupported PacBio BAM frame array type"};
    }
}

std::optional<Data::Frames> OptionalFramesTag(const TagMap& tags, TagKey key,
                                              std::string_view malformedMessage)
{
    const auto* tag{tags.Get(key)};
    if (!tag) {
        return std::nullopt;
    }
    return FramesFromTagArray(RequiredTagArray(*tag, malformedMessage));
}

}  // namespace

std::string BamRecord::FullName() const { return std::string{name_}; }

std::string BamRecord::MovieName() const
{
    const std::size_t firstSlash{name_.find('/')};
    return std::string{name_.substr(0, firstSlash)};
}

std::int32_t BamRecord::HoleNumber() const
{
    const std::string_view name{name_};
    const std::size_t firstSlash{name.find('/')};
    if (firstSlash == std::string::npos) {
        throw std::runtime_error{"Malformed PacBio BAM read name: " + name_};
    }
    const std::size_t secondSlash{name.find('/', firstSlash + 1)};
    const std::string_view holeField{name.substr(firstSlash + 1, secondSlash - (firstSlash + 1))};
    return std::stoi(std::string{holeField});
}

std::int32_t BamRecord::QueryStart() const
{
    if (const auto* tag{tags_.Get(QS_TAG)}; tag) {
        return TagToInt32(*tag);
    }
    return ParseQueryInterval(name_).first;
}

std::int32_t BamRecord::QueryEnd() const
{
    if (const auto* tag{tags_.Get(QE_TAG)}; tag) {
        return TagToInt32(*tag);
    }
    return ParseQueryInterval(name_).second;
}

std::string BamRecord::ReadGroupId() const
{
    const TagValue& tag{RequiredTag(tags_, RG_TAG, "PacBio BAM record is missing RG tag")};
    if (const auto* readGroupId{std::get_if<std::string>(&tag)}; readGroupId) {
        return *readGroupId;
    }
    throw std::runtime_error{"PacBio BAM RG tag is not a string"};
}

std::optional<Data::LocalContextFlags> BamRecord::LocalContextFlags() const
{
    if (const auto value{OptionalIntTag(tags_, CX_TAG)}; value) {
        return static_cast<Data::LocalContextFlags>(*value);
    }
    return std::nullopt;
}

void BamRecord::LocalContextFlags(Data::LocalContextFlags flags)
{
    tags_.Set(CX_TAG, static_cast<std::int64_t>(flags));
}

Data::SNR BamRecord::SignalToNoise() const
{
    const TagValue& tag{RequiredTag(tags_, SN_TAG, "PacBio BAM record is missing sn tag")};
    const TagArray& array{RequiredTagArray(tag, "PacBio BAM sn tag is malformed")};
    if (array.ElementType() != 'f') {
        throw std::runtime_error{"PacBio BAM sn tag is malformed"};
    }
    return Data::SNR{DecodePodArray<float>(array)};
}

Data::Accuracy BamRecord::ReadAccuracy() const
{
    return TagToFloat(RequiredTag(tags_, RQ_TAG, "PacBio BAM record is missing rq tag"));
}

std::optional<Data::Frames> BamRecord::PulseWidth() const
{
    return OptionalFramesTag(tags_, PW_TAG, "PacBio BAM pw tag is malformed");
}

std::optional<Data::Frames> BamRecord::IPD() const
{
    return OptionalFramesTag(tags_, IP_TAG, "PacBio BAM ip tag is malformed");
}

std::optional<std::int32_t> BamRecord::WallStart() const { return OptionalIntTag(tags_, WS_TAG); }

std::optional<std::int32_t> BamRecord::WallEnd() const { return OptionalIntTag(tags_, WE_TAG); }

}  // namespace Samoa
}  // namespace PacBio
