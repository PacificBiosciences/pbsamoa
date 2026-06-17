#ifndef PBSAMOA_SRC_CRAM_STRUCTS_HPP
#define PBSAMOA_SRC_CRAM_STRUCTS_HPP

// Internal Cram helpers shared between CramReader.cpp and CramWriter.cpp.

#include <atomic>

#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief Build a work-stealing loop callable for use with Parallel::Dispatch.
///
/// The returned callable runs:
///   while (true) {
///     IndexType idx = next->fetch_add(1, relaxed);
///     if (idx >= total) break;
///     fn(idx);
///   }
///
/// \param next   Pointer to the shared atomic counter (stack-allocated at call site).
/// \param total  Total number of items to process.
/// \param fn     Callable invoked with each item index.
template <typename IndexType, typename Fn>
auto MakeWorkStealingTask(std::atomic<IndexType>* next, IndexType total, Fn fn)
{
    return [next, total, fn](std::int32_t) {
        while (true) {
            const IndexType idx = next->fetch_add(1, std::memory_order_relaxed);
            if (idx >= total) {
                break;
            }
            fn(idx);
        }
    };
}

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_SRC_CRAM_STRUCTS_HPP
