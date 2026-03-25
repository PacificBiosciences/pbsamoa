#ifndef PBSAMOA_SRC_READERUTILS_HPP
#define PBSAMOA_SRC_READERUTILS_HPP

#include <functional>
#include <optional>

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

}  // namespace detail
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_SRC_READERUTILS_HPP
