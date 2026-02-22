#include <pbsamoa/core/BamRecord.hpp>

#include <pbsamoa/core/CigarClipping.hpp>
#include <pbsamoa/core/Sequence.hpp>
#include <pbsamoa/core/TagClipping.hpp>

#include <algorithm>
#include <format>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <utility>

namespace PacBio {
namespace Samoa {

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

bool BamRecord::IsMapped() const { return (flag_ & 0x4) == 0; }

bool BamRecord::IsReverseStrand() const { return (flag_ & 0x10) != 0; }

bool BamRecord::IsPrimary() const { return (flag_ & 0x900) == 0; }

std::int32_t BamRecord::ReferenceEnd() const { return pos_ + ReferenceLength(cigar_); }

// --- serialization ---

std::vector<std::byte> BamRecord::SerializeToBam() const
{
    // Compute variable-length field sizes
    if (std::size(name_) > 254) {
        throw std::invalid_argument{
            "BamRecord::SerializeToBam: read name exceeds BAM limit of 254 characters"};
    }
    const std::uint8_t nameLen{static_cast<std::uint8_t>(std::size(name_) + 1)};  // +NUL
    if (std::size(cigar_) > std::numeric_limits<std::uint16_t>::max()) {
        throw std::invalid_argument{
            "BamRecord::SerializeToBam: CIGAR exceeds BAM limit of 65535 operations"};
    }
    const std::uint16_t nCigarOp{static_cast<std::uint16_t>(std::size(cigar_))};
    const std::uint32_t seqLen{static_cast<std::uint32_t>(std::size(sequence_))};
    if (!std::empty(qualities_) && (std::size(qualities_) != seqLen)) {
        throw std::invalid_argument{
            std::format("BamRecord::SerializeToBam: quality length ({}) != sequence length ({})",
                        std::size(qualities_), seqLen)};
    }
    const std::uint32_t packedSeqLen{(seqLen + 1) / 2};

    // Compute bin
    std::uint16_t bin{0};
    if (pos_ >= 0 && !std::empty(cigar_)) {
        const std::int32_t refEnd{ReferenceEnd()};
        bin = Reg2Bin(pos_, (refEnd > pos_) ? refEnd : (pos_ + 1));
    } else if (pos_ >= 0) {
        bin = Reg2Bin(pos_, pos_ + 1);
    } else {
        bin = 4680;  // unmapped
    }

    // Compute tag size without allocating
    const std::size_t tagSize{SerializedBamSize(tags_)};

    // Total size: 32 (fixed) + nameLen + 4*nCigarOp + packedSeqLen + seqLen + tagSize
    const std::size_t totalSize{32U + nameLen + 4U * nCigarOp + packedSeqLen + seqLen + tagSize};

    std::vector<std::byte> result(totalSize);
    std::byte* p{std::data(result)};

    // Fixed fields (32 bytes)
    std::ranges::copy_n(reinterpret_cast<const std::byte*>(&refId_), sizeof(refId_), p + 0);
    std::ranges::copy_n(reinterpret_cast<const std::byte*>(&pos_), sizeof(pos_), p + 4);
    p[8] = static_cast<std::byte>(nameLen);
    p[9] = static_cast<std::byte>(mapQ_);
    std::ranges::copy_n(reinterpret_cast<const std::byte*>(&bin), sizeof(bin), p + 10);
    std::ranges::copy_n(reinterpret_cast<const std::byte*>(&nCigarOp), sizeof(nCigarOp), p + 12);
    std::ranges::copy_n(reinterpret_cast<const std::byte*>(&flag_), sizeof(flag_), p + 14);
    std::ranges::copy_n(reinterpret_cast<const std::byte*>(&seqLen), sizeof(seqLen), p + 16);
    std::ranges::copy_n(reinterpret_cast<const std::byte*>(&nextRefId_), sizeof(nextRefId_),
                        p + 20);
    std::ranges::copy_n(reinterpret_cast<const std::byte*>(&nextPos_), sizeof(nextPos_), p + 24);
    std::ranges::copy_n(reinterpret_cast<const std::byte*>(&tlen_), sizeof(tlen_), p + 28);

    std::size_t offset{32};

    // Read name (NUL-terminated)
    std::ranges::copy_n(reinterpret_cast<const std::byte*>(std::data(name_)), std::size(name_),
                        p + offset);
    p[offset + std::size(name_)] = std::byte{0};
    offset += nameLen;

    // CIGAR (array of uint32)
    for (const CigarOp& op : cigar_) {
        const std::uint32_t raw{op.RawValue()};
        std::ranges::copy_n(reinterpret_cast<const std::byte*>(&raw), sizeof(raw), p + offset);
        offset += 4;
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

    return result;
}

// --- clipping ---

BamRecord& BamRecord::Clip(ClipType type, std::int32_t start, std::int32_t end,
                           bool exciseFlankingInserts)
{
    static const auto DEFAULT_CLIPPER{TagClipper::PacBioDefault()};
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
    std::int32_t origQStart{0};
    std::int32_t origQEnd{static_cast<std::int32_t>(origSeqLen)};
    if (const auto* qs{tags_.Get(TagKey{'q', 's'})}) {
        origQStart = static_cast<std::int32_t>(std::get<std::int64_t>(*qs));
    }
    if (const auto* qe{tags_.Get(TagKey{'q', 'e'})}) {
        origQEnd = static_cast<std::int32_t>(std::get<std::int64_t>(*qe));
    }

    // Compute ClipResult
    ClipResult cr;
    if (type == ClipType::CLIP_TO_QUERY) {
        cr = ClipCigarToQuery(cigar_, start, end, origQStart, origQEnd, pos_, isReverse);
    } else {
        cr = ClipCigarToReference(cigar_, start, end, pos_, isReverse, exciseFlankingInserts);
    }

    // Clip tags BEFORE modifying sequence (basemods reads old sequence)
    clipper.ClipTags(tags_, cr.clipOffset, cr.clipLength, origSeqLen, sequence_);

    // Clip sequence and qualities (in-place to avoid extra allocation)
    sequence_.erase(0, cr.clipOffset);
    sequence_.resize(cr.clipLength);
    if (!std::empty(qualities_)) {
        qualities_.erase(std::begin(qualities_),
                         std::begin(qualities_) + static_cast<std::ptrdiff_t>(cr.clipOffset));
        qualities_.resize(cr.clipLength);
    }

    // Update CIGAR and position
    cigar_ = std::move(cr.cigar);
    pos_ = cr.newPos;

    // Update qs/qe tags
    if (type == ClipType::CLIP_TO_QUERY) {
        if (tags_.Contains(TagKey{'q', 's'})) {
            tags_.Set(TagKey{'q', 's'}, std::int64_t{start});
            tags_.Set(TagKey{'q', 'e'}, std::int64_t{end});
        }
    } else {
        // For reference clipping, qs/qe shift by clipOffset
        if (tags_.Contains(TagKey{'q', 's'})) {
            const std::int32_t newQs{origQStart + static_cast<std::int32_t>(cr.clipOffset)};
            tags_.Set(TagKey{'q', 's'}, std::int64_t{newQs});
            tags_.Set(TagKey{'q', 'e'},
                      std::int64_t{newQs + static_cast<std::int32_t>(cr.clipLength)});
        }
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

}  // namespace Samoa
}  // namespace PacBio
