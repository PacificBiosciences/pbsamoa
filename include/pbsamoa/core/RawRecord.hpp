#ifndef PBSAMOA_CORE_RAWRECORD_HPP
#define PBSAMOA_CORE_RAWRECORD_HPP

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/Endian.hpp>
#include <pbsamoa/core/Sequence.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <algorithm>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <cstddef>
#include <cstdint>
#include <cstring>

namespace PacBio {
namespace Samoa {

/// \brief Owning byte view into raw BAM record bytes with decode-on-demand.
///
/// Copies raw record bytes on construction. All accessors decode on-the-fly
/// from the owned buffer. The view owns its data and has independent lifetime
/// — safe to pass across thread boundaries.
///
/// The span covers the record bytes after block_size (starts at refID).
class RawRecord
{
public:
    explicit RawRecord(std::span<const std::byte> data);

    // --- fixed fields ---
    std::int32_t RefId() const;
    std::int32_t Pos() const;
    std::uint8_t NameLength() const;
    std::uint8_t MapQ() const;
    std::uint16_t Bin() const;
    std::uint16_t CigarOpCount() const;
    std::uint16_t Flag() const;
    std::uint32_t SeqLength() const;
    std::int32_t NextRefId() const;
    std::int32_t NextPos() const;
    std::int32_t Tlen() const;

    // --- variable-length fields ---
    std::string_view Name() const;
    CigarView CigarOps() const;
    SequenceView Seq() const;
    std::span<const std::uint8_t> Qual() const;
    std::span<const std::byte> AuxData() const;

    // --- derived ---
    bool IsMapped() const;
    bool IsReverseStrand() const;
    bool IsPrimary() const;
    std::int64_t ReferenceLength() const;
    std::int64_t QueryLength() const;

    TagMap ParseTags() const;
    BamRecord ToOwned() const;
    BamRecord ToOwned(const DropTags& filter) const;
    BamRecord ToOwned(const KeepTags& filter) const;
    std::span<const std::byte> RawData() const;

private:
    std::vector<std::byte> data_;
    std::vector<CigarOp> cigar_;
    bool cgExpanded_{false};  // true if a CG-tag long CIGAR was expanded into cigar_

    std::size_t CigarOffset() const;
    std::size_t SeqOffset() const;
    std::size_t QualOffset() const;
    std::size_t AuxOffset() const;

    /// \brief Expand a CG:B,I long-CIGAR placeholder into cigar_ (SAMv1 §4.2.2).
    /// Out-of-line: only invoked when the cheap placeholder guard in the ctor matches.
    void ExpandLongCigarFromCgTag();
};

// --- inline implementations ---

inline RawRecord::RawRecord(std::span<const std::byte> data) : data_{data.begin(), data.end()}
{
    // Validate minimum size for fixed fields (32 bytes)
    if (std::size(data_) < 32U) {
        throw std::invalid_argument{
            "RawRecord: data too small for BAM fixed fields (need >= 32 bytes)"};
    }

    // Validate name length (l_read_name includes NUL, must be >= 1)
    const std::uint8_t nameLen{NameLength()};
    if (nameLen == 0) {
        throw std::invalid_argument{"RawRecord: l_read_name is 0 (must be >= 1)"};
    }

    // Validate that all variable-length fields fit within the buffer
    const std::size_t auxOffset{AuxOffset()};
    if (auxOffset > std::size(data_)) {
        throw std::invalid_argument{
            "RawRecord: truncated record — variable-length fields exceed buffer size"};
    }
    if (data_[32 + nameLen - 1] != std::byte{0}) {
        throw std::invalid_argument{"RawRecord: read name is not NUL-terminated"};
    }

    // Copy CIGAR directly — CigarOp stores BAM-native uint32 layout and
    // Endian.hpp guarantees a little-endian platform.
    static_assert(sizeof(CigarOp) == sizeof(std::uint32_t));
    if (const std::uint16_t count{CigarOpCount()}; count > 0) {
        cigar_.resize(count);
        std::memcpy(cigar_.data(), std::data(data_) + CigarOffset(),
                    static_cast<std::size_t>(count) * sizeof(CigarOp));
    }

    // CG-tag long-CIGAR placeholder check (SAMv1 §4.2.2; htslib bam_tag2cigar). A record
    // with >65535 ops stores 'kSmN' (n_cigar_op=2, cigar[0] = soft-clip of the whole read)
    // and the real CIGAR in a CG:B,I tag. The guard is cheap and almost always false;
    // the actual expansion is out-of-line.
    if (cigar_.size() == 2U && RefId() >= 0 && Pos() >= 0 &&
        cigar_[0].RawValue() == ((SeqLength() << 4U) | std::to_underlying(CigarOpType::S))) {
        ExpandLongCigarFromCgTag();
    }
}

// --- fixed fields ---

inline std::int32_t RawRecord::RefId() const { return ReadI32LE(std::data(data_)); }

inline std::int32_t RawRecord::Pos() const { return ReadI32LE(std::data(data_) + 4); }

inline std::uint8_t RawRecord::NameLength() const { return static_cast<std::uint8_t>(data_[8]); }

inline std::uint8_t RawRecord::MapQ() const { return static_cast<std::uint8_t>(data_[9]); }

inline std::uint16_t RawRecord::Bin() const { return ReadU16LE(std::data(data_) + 10); }

inline std::uint16_t RawRecord::CigarOpCount() const { return ReadU16LE(std::data(data_) + 12); }

inline std::uint16_t RawRecord::Flag() const { return ReadU16LE(std::data(data_) + 14); }

inline std::uint32_t RawRecord::SeqLength() const { return ReadU32LE(std::data(data_) + 16); }

inline std::int32_t RawRecord::NextRefId() const { return ReadI32LE(std::data(data_) + 20); }

inline std::int32_t RawRecord::NextPos() const { return ReadI32LE(std::data(data_) + 24); }

inline std::int32_t RawRecord::Tlen() const { return ReadI32LE(std::data(data_) + 28); }

// --- variable-length field offsets ---

inline std::size_t RawRecord::CigarOffset() const { return 32U + NameLength(); }

inline std::size_t RawRecord::SeqOffset() const { return CigarOffset() + 4U * CigarOpCount(); }

inline std::size_t RawRecord::QualOffset() const { return SeqOffset() + (SeqLength() + 1) / 2; }

inline std::size_t RawRecord::AuxOffset() const { return QualOffset() + SeqLength(); }

// --- variable-length fields ---

inline std::string_view RawRecord::Name() const
{
    return {
        reinterpret_cast<const char*>(std::data(data_) + 32),
        static_cast<std::size_t>(NameLength() - 1),
    };
}

inline CigarView RawRecord::CigarOps() const { return cigar_; }

inline SequenceView RawRecord::Seq() const
{
    const std::size_t offset{SeqOffset()};
    const std::uint32_t seqLength{SeqLength()};
    const std::size_t packedLength{(seqLength + 1) / 2};
    return {std::span<const std::byte>{data_}.subspan(offset, packedLength), seqLength};
}

inline std::span<const std::uint8_t> RawRecord::Qual() const
{
    const std::size_t offset{QualOffset()};
    return {reinterpret_cast<const std::uint8_t*>(std::data(data_) + offset), SeqLength()};
}

inline std::span<const std::byte> RawRecord::AuxData() const
{
    const std::size_t offset{AuxOffset()};
    return std::span<const std::byte>{data_}.subspan(offset);
}

// --- derived ---

inline bool RawRecord::IsMapped() const { return (Flag() & 0x4) == 0; }

inline bool RawRecord::IsReverseStrand() const { return (Flag() & 0x10) != 0; }

inline bool RawRecord::IsPrimary() const { return (Flag() & 0x900) == 0; }

inline std::int64_t RawRecord::ReferenceLength() const
{
    return ::PacBio::Samoa::ReferenceLength(CigarOps());
}

inline std::int64_t RawRecord::QueryLength() const
{
    return ::PacBio::Samoa::QueryLength(CigarOps());
}

inline TagMap RawRecord::ParseTags() const
{
    TagMap tags{ParseTagsFromBam(AuxData())};
    if (cgExpanded_) {
        // The CG tag is now redundant: its contents live in the expanded cigar_.
        tags.Remove(TagKey{'C', 'G'});
    }
    return tags;
}

inline std::span<const std::byte> RawRecord::RawData() const { return data_; }

namespace detail {

constexpr bool HasStoredQuality(std::uint8_t quality);

template <typename Filter>
BamRecord ToOwnedFiltered(const RawRecord& raw, const Filter& filter,
                          bool (*keep)(const Filter&, TagKey));

bool KeepDroppedTag(const DropTags& filter, TagKey key);
bool KeepKeptTag(const KeepTags& filter, TagKey key);

}  // namespace detail

// --- ToOwned ---

inline BamRecord RawRecord::ToOwned() const
{
    BamRecord record;
    const std::span<const std::uint8_t> qualities{Qual()};
    const bool hasStoredQualities{std::ranges::any_of(qualities, detail::HasStoredQuality)};

    record.Name(std::string{Name()})
        .Flag(Flag())
        .RefId(RefId())
        .Pos(Pos())
        .MapQ(MapQ())
        .Cigar(std::vector<CigarOp>{CigarOps().begin(), CigarOps().end()})
        .NextRefId(NextRefId())
        .NextPos(NextPos())
        .Tlen(Tlen())
        .Sequence(Seq().ToString());

    if (hasStoredQualities) {
        record.Qualities(std::vector<std::uint8_t>{qualities.begin(), qualities.end()});
    }

    record.Tags(ParseTags());
    return record;
}

namespace detail {

constexpr bool HasStoredQuality(std::uint8_t quality) { return quality != 0xFF; }

template <typename Filter>
BamRecord ToOwnedFiltered(const RawRecord& raw, const Filter& filter,
                          bool (*keep)(const Filter&, TagKey))
{
    BamRecord record{raw.ToOwned()};
    TagMap filtered;
    const TagMap& tags{record.Tags()};
    for (const auto& [key, value] : tags.Entries()) {
        if (keep(filter, key)) {
            filtered.Append(key, value);
        }
    }
    record.Tags(std::move(filtered));
    return record;
}

inline bool KeepDroppedTag(const DropTags& filter, TagKey key) { return !filter.ShouldDrop(key); }

inline bool KeepKeptTag(const KeepTags& filter, TagKey key) { return filter.ShouldKeep(key); }

}  // namespace detail

inline BamRecord RawRecord::ToOwned(const DropTags& filter) const
{
    return detail::ToOwnedFiltered(*this, filter, detail::KeepDroppedTag);
}

inline BamRecord RawRecord::ToOwned(const KeepTags& filter) const
{
    return detail::ToOwnedFiltered(*this, filter, detail::KeepKeptTag);
}

/// \brief Owns decompressed buffer(s) and provides access to raw record data.
///
/// The batch is the unit of parallel processing. Users access raw record
/// bytes via RecordData(), constructing RawRecord on demand for
/// decode access. The batch is invalidated when destroyed.
class RawRecordBatch
{
public:
    /// \brief Describes where a record lives within the buffer.
    struct RecordExtent
    {
        std::uint32_t offset;  // byte offset of record data within buffer
        std::uint32_t size;    // byte size of record data (block_size value)
    };

    /// \brief Construct from a buffer and record extents.
    RawRecordBatch(std::vector<std::byte> buffer, std::vector<RecordExtent> extents);

    /// \brief Raw bytes for record at index \p i.
    std::span<const std::byte> RecordData(std::size_t i) const;

    std::size_t RecordCount() const;
    std::size_t BufferSize() const;

private:
    std::vector<std::byte> buffer_;
    std::vector<RecordExtent> extents_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_CORE_RAWRECORD_HPP
