#include <pbsamoa/core/RawRecord.hpp>

#include <span>
#include <utility>
#include <variant>
#include <vector>

#include <cstring>

namespace PacBio {
namespace Samoa {

// --- RawRecordView ---

void RawRecordView::ExpandLongCigarFromCgTag()
{
    // The ctor has already verified the kSmN placeholder shape. Expand only when a
    // well-formed CG:B,I tag is present, mirroring htslib bam_tag2cigar's guards so that
    // a record carrying an unrelated CG tag is left untouched.
    const TagMap tags{ParseTagsFromBam(AuxData())};
    const TagValue* value{tags.Get(TagKey{'C', 'G'})};
    if (value == nullptr) {
        return;
    }
    const TagArray* array{std::get_if<TagArray>(value)};
    if ((array == nullptr) || ((array->ElementType() != 'I') && (array->ElementType() != 'i'))) {
        return;
    }

    // Reject counts below the placeholder size or beyond the 28-bit CIGAR field (htslib
    // uses 1<<29 as the guard since the count itself is a uint32).
    const std::uint32_t count{array->Count()};
    if ((count < CigarOpCount()) || (count >= (1U << 29U))) {
        return;
    }

    const std::span<const std::byte> bytes{array->Data()};
    const std::size_t byteCount{static_cast<std::size_t>(count) * sizeof(CigarOp)};
    if (std::size(bytes) < byteCount) {
        return;
    }

    // CG payload is BAM-native little-endian uint32 ops, identical to CigarOp's layout.
    static_assert(sizeof(CigarOp) == sizeof(std::uint32_t));
    std::vector<CigarOp> expanded(count);
    std::memcpy(expanded.data(), std::data(bytes), byteCount);

    cigar_ = std::move(expanded);
    cgExpanded_ = true;
}

// --- RawRecordBatch ---

RawRecordBatch::RawRecordBatch(std::vector<std::byte> buffer, std::vector<RecordExtent> extents)
    : buffer_{std::move(buffer)}, extents_{std::move(extents)}
{
}

std::span<const std::byte> RawRecordBatch::RecordData(std::size_t i) const
{
    const RecordExtent& extent{extents_[i]};
    return std::span<const std::byte>{buffer_}.subspan(extent.offset, extent.size);
}

std::size_t RawRecordBatch::RecordCount() const { return std::size(extents_); }

std::size_t RawRecordBatch::BufferSize() const { return std::size(buffer_); }

}  // namespace Samoa
}  // namespace PacBio
