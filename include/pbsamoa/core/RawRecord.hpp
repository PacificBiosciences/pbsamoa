#ifndef PBSAMOA_CORE_RAWRECORD_HPP
#define PBSAMOA_CORE_RAWRECORD_HPP

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/Endian.hpp>
#include <pbsamoa/core/Sequence.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <algorithm>
#include <memory>
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

/// \brief Non-owning view into raw BAM record bytes with decode-on-demand.
///
/// The source bytes must outlive the view and every span that the view returns. Each
/// accessor decodes its field at the time of the call, not in advance. The view copies
/// CIGAR words into aligned storage because BAM does not guarantee that their source
/// address is 4-byte aligned.
///
/// The span covers the record bytes after block_size (starts at refID).
class RawRecordView
{
public:
    explicit RawRecordView(std::span<const std::byte> data);

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

    /// \brief The record's CIGAR ops.
    ///
    /// The returned span points into storage owned by this view. If the caller moves or
    /// destroys the view, the span becomes invalid. Call CigarOps() at the point of use,
    /// and do not save the result across a move.
    CigarView CigarOps() const;
    SequenceView Seq() const;
    std::span<const std::uint8_t> Qual() const;
    std::span<const std::byte> AuxData() const;

    // --- derived ---
    bool IsMapped() const;
    bool IsReverseStrand() const;
    bool IsPrimary() const;
    bool IsSecondary() const;
    bool IsSupplementary() const;
    std::int64_t ReferenceLength() const;
    std::int64_t QueryLength() const;

    TagMap ParseTags() const;

    /// \brief Parse aux tags, keeping only keys accepted by \p keep.
    ///
    /// Walks the raw BAM aux bytes once, skipping tags that fail the predicate.
    /// Replicates ParseTags()'s CG-tag suppression: when cgExpanded_ is true the
    /// CG tag is never materialised (its data already lives in cigar_).
    template <typename Filter>
    TagMap ParseTagsFiltered(const Filter& filter, bool (*keep)(const Filter&, TagKey)) const;

    BamRecord ToOwned() const;
    BamRecord ToOwned(const DropTags& filter) const;
    BamRecord ToOwned(const KeepTags& filter) const;
    std::span<const std::byte> RawData() const;

protected:
    void Rebind(std::span<const std::byte> data) noexcept;

private:
    std::span<const std::byte> data_;
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

/// \brief Owning raw BAM record with decode-on-demand accessors.
///
/// Copies raw record bytes on construction and is safe to pass across thread
/// boundaries. Views returned by Name()/Seq()/Qual()/AuxData() remain valid after the
/// caller moves the record. CigarOps() follows RawRecordView's move-invalidation rule.
class RawRecord : private RawRecordView
{
public:
    explicit RawRecord(std::span<const std::byte> data);

    RawRecord(const RawRecord& other);
    RawRecord& operator=(const RawRecord& other);
    RawRecord(RawRecord&& other) noexcept;
    RawRecord& operator=(RawRecord&& other) noexcept;
    ~RawRecord() = default;

    using RawRecordView::AuxData;
    using RawRecordView::Bin;
    using RawRecordView::CigarOpCount;
    using RawRecordView::CigarOps;
    using RawRecordView::Flag;
    using RawRecordView::IsMapped;
    using RawRecordView::IsPrimary;
    using RawRecordView::IsReverseStrand;
    using RawRecordView::IsSecondary;
    using RawRecordView::IsSupplementary;
    using RawRecordView::MapQ;
    using RawRecordView::Name;
    using RawRecordView::NameLength;
    using RawRecordView::NextPos;
    using RawRecordView::NextRefId;
    using RawRecordView::ParseTags;
    using RawRecordView::ParseTagsFiltered;
    using RawRecordView::Pos;
    using RawRecordView::Qual;
    using RawRecordView::QueryLength;
    using RawRecordView::RawData;
    using RawRecordView::ReferenceLength;
    using RawRecordView::RefId;
    using RawRecordView::Seq;
    using RawRecordView::SeqLength;
    using RawRecordView::Tlen;
    using RawRecordView::ToOwned;

    /// \brief Borrow this record through the non-owning decoder interface.
    ///
    /// The returned view stays valid only while this record exists and the caller does
    /// not move it.
    const RawRecordView& View() const& noexcept;
    const RawRecordView& View() const&& = delete;

private:
    // A BAM record's byte count is known only at run time, so std::array cannot represent
    // it. ByteArray owns the bytes instead.
    using ByteArray = std::unique_ptr<
        std::byte[]>;  // NOLINT(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)

    static ByteArray Allocate(std::size_t size);

    ByteArray ownedData_;
};

// --- inline implementations ---

inline RawRecordView::RawRecordView(std::span<const std::byte> data) : data_{data}
{
    // Validate minimum size for fixed fields (32 bytes)
    if (std::size(data_) < 32U) {
        throw std::invalid_argument{
            "RawRecordView: data too small for BAM fixed fields (need >= 32 bytes)"};
    }

    // Validate name length (l_read_name includes NUL, must be >= 1)
    const std::uint8_t nameLen{NameLength()};
    if (nameLen == 0) {
        throw std::invalid_argument{"RawRecordView: l_read_name is 0 (must be >= 1)"};
    }

    // Validate that all variable-length fields fit within the buffer
    const std::size_t auxOffset{AuxOffset()};
    if (auxOffset > std::size(data_)) {
        throw std::invalid_argument{
            "RawRecordView: truncated record — variable-length fields exceed buffer size"};
    }
    if (data_[32 + nameLen - 1] != std::byte{0}) {
        throw std::invalid_argument{"RawRecordView: read name is not NUL-terminated"};
    }

    // Copy CIGAR directly — CigarOp stores BAM-native uint32 layout and
    // Endian.hpp guarantees a little-endian platform.
    static_assert(sizeof(CigarOp) == sizeof(std::uint32_t));
    if (const std::uint16_t count{CigarOpCount()}; count > 0) {
        const std::size_t byteCount{static_cast<std::size_t>(count) * sizeof(CigarOp)};
        cigar_.resize(count);
        std::memcpy(cigar_.data(), std::data(data_) + CigarOffset(), byteCount);
    }

    // CG-tag long-CIGAR placeholder check (SAMv1 §4.2.2; htslib bam_tag2cigar). A record
    // with >65535 ops stores 'kSmN' (n_cigar_op=2, cigar[0] = soft-clip of the whole read)
    // and the real CIGAR in a CG:B,I tag. The guard is cheap and almost always false;
    // the actual expansion is out-of-line.
    if (const CigarView placeholder{CigarOps()};
        (std::size(placeholder) == 2U) && (RefId() >= 0) && (Pos() >= 0) &&
        (placeholder[0].RawValue() == ((SeqLength() << 4U) | std::to_underlying(CigarOpType::S)))) {
        ExpandLongCigarFromCgTag();
    }

    // A mapped record's CIGAR must consume exactly l_seq query bases; htslib rejects a
    // mismatch (sam.c bam_read1). The check reads the CIGAR again through CigarOps(). The
    // check then sees the CG-tag expansion above, not the two-op placeholder that
    // ExpandLongCigarFromCgTag() replaced.
    const CigarView cigar{CigarOps()};
    if (!std::empty(cigar) && IsMapped() && (SeqLength() > 0) &&
        (::PacBio::Samoa::QueryLength(cigar) != static_cast<std::int64_t>(SeqLength()))) {
        throw std::invalid_argument{"RawRecordView: CIGAR query length does not match SEQ length"};
    }
}

inline void RawRecordView::Rebind(std::span<const std::byte> data) noexcept { data_ = data; }

inline RawRecord::ByteArray RawRecord::Allocate(std::size_t size)
{
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
    return std::make_unique_for_overwrite<std::byte[]>(size);
}

inline RawRecord::RawRecord(std::span<const std::byte> data)
    : RawRecordView{data}, ownedData_{Allocate(std::size(data))}
{
    std::copy_n(std::data(data), std::size(data), ownedData_.get());
    Rebind({ownedData_.get(), std::size(data)});
}

inline RawRecord::RawRecord(const RawRecord& other)
    : RawRecordView{other}, ownedData_{Allocate(std::size(other.RawData()))}
{
    const std::span<const std::byte> source{other.RawData()};
    std::copy_n(std::data(source), std::size(source), ownedData_.get());
    Rebind({ownedData_.get(), std::size(source)});
}

inline RawRecord& RawRecord::operator=(const RawRecord& other)
{
    if (this != &other) {
        RawRecord copy{other};
        *this = std::move(copy);
    }
    return *this;
}

inline RawRecord::RawRecord(RawRecord&& other) noexcept
    : RawRecordView{std::move(other)}, ownedData_{std::move(other.ownedData_)}
{
    // This move constructor does not call Rebind for `this`. Each RawRecord's data_
    // always points into its own ownedData_. The moves above keep the invariant for
    // `this`: the base class move copies the data_ span, and the unique_ptr move
    // transfers the same buffer. So this object's span and its buffer still match. The
    // moved-from `other` is different: its data_ span still points at the transferred
    // buffer, but its ownedData_ is now empty. The constructor clears that span with
    // other.Rebind({}). Otherwise, other.RawData() would return dangling bytes after this
    // object is destroyed.
    other.Rebind({});
}

inline RawRecord& RawRecord::operator=(RawRecord&& other) noexcept
{
    if (this != &other) {
        RawRecordView::operator=(std::move(other));
        ownedData_ = std::move(other.ownedData_);
        other.Rebind({});
    }
    return *this;
}

inline const RawRecordView& RawRecord::View() const& noexcept { return *this; }

// --- fixed fields ---

inline std::int32_t RawRecordView::RefId() const { return ReadI32LE(std::data(data_)); }

inline std::int32_t RawRecordView::Pos() const { return ReadI32LE(std::data(data_) + 4); }

inline std::uint8_t RawRecordView::NameLength() const
{
    return std::to_integer<std::uint8_t>(data_[8]);
}

inline std::uint8_t RawRecordView::MapQ() const { return std::to_integer<std::uint8_t>(data_[9]); }

inline std::uint16_t RawRecordView::Bin() const { return ReadU16LE(std::data(data_) + 10); }

inline std::uint16_t RawRecordView::CigarOpCount() const
{
    return ReadU16LE(std::data(data_) + 12);
}

inline std::uint16_t RawRecordView::Flag() const { return ReadU16LE(std::data(data_) + 14); }

inline std::uint32_t RawRecordView::SeqLength() const { return ReadU32LE(std::data(data_) + 16); }

inline std::int32_t RawRecordView::NextRefId() const { return ReadI32LE(std::data(data_) + 20); }

inline std::int32_t RawRecordView::NextPos() const { return ReadI32LE(std::data(data_) + 24); }

inline std::int32_t RawRecordView::Tlen() const { return ReadI32LE(std::data(data_) + 28); }

// --- variable-length field offsets ---

inline std::size_t RawRecordView::CigarOffset() const { return 32U + NameLength(); }

inline std::size_t RawRecordView::SeqOffset() const { return CigarOffset() + 4U * CigarOpCount(); }

inline std::size_t RawRecordView::QualOffset() const
{
    // This function widens the value before it adds 1. l_seq comes from untrusted input,
    // and at UINT32_MAX a uint32 computation of l_seq + 1 wraps to 0. That wrap would
    // produce an auxOffset small enough to pass the constructor's bounds check.
    // Endian.hpp's static assertion requires a 64-bit std::size_t, and this function's
    // correctness depends on that requirement.
    const std::size_t seqLength{SeqLength()};
    return SeqOffset() + ((seqLength + 1U) / 2U);
}

inline std::size_t RawRecordView::AuxOffset() const { return QualOffset() + SeqLength(); }

// --- variable-length fields ---

inline std::string_view RawRecordView::Name() const
{
    return {
        reinterpret_cast<const char*>(std::data(data_) + 32),
        static_cast<std::size_t>(NameLength() - 1),
    };
}

inline CigarView RawRecordView::CigarOps() const { return cigar_; }

inline SequenceView RawRecordView::Seq() const
{
    const std::size_t offset{SeqOffset()};
    const std::size_t seqLength{SeqLength()};
    const std::size_t packedLength{(seqLength + 1U) / 2U};
    return {std::span<const std::byte>{data_}.subspan(offset, packedLength), SeqLength()};
}

inline std::span<const std::uint8_t> RawRecordView::Qual() const
{
    const std::size_t offset{QualOffset()};
    return {reinterpret_cast<const std::uint8_t*>(std::data(data_) + offset), SeqLength()};
}

inline std::span<const std::byte> RawRecordView::AuxData() const
{
    const std::size_t offset{AuxOffset()};
    return std::span<const std::byte>{data_}.subspan(offset);
}

// --- derived ---

inline bool RawRecordView::IsMapped() const { return (Flag() & 0x4) == 0; }

inline bool RawRecordView::IsReverseStrand() const { return (Flag() & 0x10) != 0; }

inline bool RawRecordView::IsPrimary() const { return (Flag() & 0x900) == 0; }

inline bool RawRecordView::IsSecondary() const { return (Flag() & 0x100) != 0; }

inline bool RawRecordView::IsSupplementary() const { return (Flag() & 0x800) != 0; }

inline std::int64_t RawRecordView::ReferenceLength() const
{
    return ::PacBio::Samoa::ReferenceLength(CigarOps());
}

inline std::int64_t RawRecordView::QueryLength() const
{
    return ::PacBio::Samoa::QueryLength(CigarOps());
}

inline TagMap RawRecordView::ParseTags() const
{
    TagMap tags{ParseTagsFromBam(AuxData())};
    if (cgExpanded_) {
        // The CG tag is now redundant: its contents live in the expanded cigar_.
        tags.Remove(TagKey{'C', 'G'});
    }
    return tags;
}

inline std::span<const std::byte> RawRecordView::RawData() const { return data_; }

template <typename Filter>
TagMap RawRecordView::ParseTagsFiltered(const Filter& filter,
                                        bool (*keep)(const Filter&, TagKey)) const
{
    const TagKey cgKey{'C', 'G'};
    const TagMap parsed{ParseTagsFromBam(AuxData())};
    TagMap tags;
    for (const auto& [key, value] : parsed.Entries()) {
        // Replicate ParseTags() CG-tag suppression: when the CG tag's CIGAR was
        // expanded into cigar_, the tag is semantically stale. Never include it.
        if (cgExpanded_ && (key == cgKey)) {
            continue;
        }
        if (keep(filter, key)) {
            tags.Append(key, value);
        }
    }
    return tags;
}

namespace detail {

constexpr bool HasStoredQuality(std::uint8_t quality);

template <typename Filter>
BamRecord ToOwnedFiltered(const RawRecordView& raw, const Filter& filter,
                          bool (*keep)(const Filter&, TagKey));

bool KeepDroppedTag(const DropTags& filter, TagKey key);
bool KeepKeptTag(const KeepTags& filter, TagKey key);

}  // namespace detail

// --- ToOwned ---

inline BamRecord RawRecordView::ToOwned() const
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
BamRecord ToOwnedFiltered(const RawRecordView& raw, const Filter& filter,
                          bool (*keep)(const Filter&, TagKey))
{
    BamRecord record;
    const std::span<const std::uint8_t> qualities{raw.Qual()};
    const bool hasStoredQualities{std::ranges::any_of(qualities, HasStoredQuality)};

    record.Name(std::string{raw.Name()})
        .Flag(raw.Flag())
        .RefId(raw.RefId())
        .Pos(raw.Pos())
        .MapQ(raw.MapQ())
        .Cigar(std::vector<CigarOp>{raw.CigarOps().begin(), raw.CigarOps().end()})
        .NextRefId(raw.NextRefId())
        .NextPos(raw.NextPos())
        .Tlen(raw.Tlen())
        .Sequence(raw.Seq().ToString());

    if (hasStoredQualities) {
        record.Qualities(std::vector<std::uint8_t>{qualities.begin(), qualities.end()});
    }

    record.Tags(raw.ParseTagsFiltered(filter, keep));
    return record;
}

inline bool KeepDroppedTag(const DropTags& filter, TagKey key) { return !filter.ShouldDrop(key); }

inline bool KeepKeptTag(const KeepTags& filter, TagKey key) { return filter.ShouldKeep(key); }

}  // namespace detail

inline BamRecord RawRecordView::ToOwned(const DropTags& filter) const
{
    return detail::ToOwnedFiltered(*this, filter, detail::KeepDroppedTag);
}

inline BamRecord RawRecordView::ToOwned(const KeepTags& filter) const
{
    return detail::ToOwnedFiltered(*this, filter, detail::KeepKeptTag);
}

/// \brief Owns decompressed buffer(s) and provides access to raw record data.
///
/// The batch is the unit of parallel processing. Users access raw record
/// bytes through View() and do not copy them. When the caller destroys the batch, the
/// view and every span obtained from it become invalid.
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

    /// \brief Non-owning decode view for record at index \p i.
    ///
    /// View() is inline so the compiler can place its field validation and CIGAR copy
    /// directly in the caller's per-record loop. The rest of the decode path follows the
    /// same pattern.
    RawRecordView View(std::size_t i) const;

    std::size_t RecordCount() const;
    std::size_t BufferSize() const;

private:
    std::vector<std::byte> buffer_;
    std::vector<RecordExtent> extents_;
};

inline RawRecordView RawRecordBatch::View(std::size_t i) const
{
    return RawRecordView{RecordData(i)};
}

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_CORE_RAWRECORD_HPP
