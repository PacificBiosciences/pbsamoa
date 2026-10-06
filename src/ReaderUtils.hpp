#ifndef PBSAMOA_SRC_READERUTILS_HPP
#define PBSAMOA_SRC_READERUTILS_HPP

#include <pbsamoa/core/RawRecord.hpp>

#include <algorithm>
#include <functional>
#include <optional>
#include <utility>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {
namespace detail {

template <typename ReaderT, typename RecordT, typename ReadNextFn>
void AdvanceReaderIterator(ReaderT*& reader, std::optional<RecordT>& current, ReadNextFn readNext)
{
    current = std::invoke(readNext, reader);
    if (!current) {
        reader = nullptr;
    }
}

template <typename ReaderT, typename RecordT>
void AdvanceReaderIterator(ReaderT*& reader, std::optional<RecordT>& current)
{
    AdvanceReaderIterator(reader, current, &ReaderT::ReadRecord);
}

// The ByteLimit default. A batch with a larger limit grows as records come in.
inline constexpr std::size_t MAX_BATCH_RESERVE{std::size_t{256} * 1024 * 1024};

// Copies records from reader.ReadRecord() into one batch. The batch ends after the
// record that reaches maxBytes, or at maxRecords.
template <typename ReaderT>
std::optional<RawRecordBatch> ReadBatchFromReader(ReaderT& reader, const std::size_t maxBytes,
                                                  const std::size_t maxRecords)
{
    std::vector<std::byte> buffer;
    buffer.reserve(std::min(maxBytes, MAX_BATCH_RESERVE));
    std::vector<RawRecordBatch::RecordExtent> extents;

    while (std::size(extents) < maxRecords) {
        if (!std::empty(extents) && (std::size(buffer) >= maxBytes)) {
            break;
        }

        const auto record{reader.ReadRecord()};
        if (!record) {
            break;
        }
        const auto raw{record->RawData()};
        extents.emplace_back(static_cast<std::uint32_t>(std::size(buffer)),
                             static_cast<std::uint32_t>(std::size(raw)));
        buffer.insert(std::end(buffer), std::begin(raw), std::end(raw));
    }

    if (std::empty(extents)) {
        return std::nullopt;
    }
    return RawRecordBatch{std::move(buffer), std::move(extents)};
}

}  // namespace detail
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_SRC_READERUTILS_HPP
