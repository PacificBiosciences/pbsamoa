#include <pbsamoa/core/RawRecord.hpp>

#include <utility>

namespace PacBio {
namespace Samoa {

// RawRecord: all implementations are inline in the header.

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
