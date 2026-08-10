#ifndef PBSAMOA_SRC_PARALLELUTILS_HPP
#define PBSAMOA_SRC_PARALLELUTILS_HPP

#include <algorithm>

#include <cstddef>

namespace PacBio {
namespace Samoa {
namespace detail {

// Dispatch by chunk so workers reuse each chunk's output buffer.
inline constexpr std::size_t PARALLEL_CHUNKS_PER_WORKER{4};
inline constexpr std::size_t MAX_PARALLEL_CHUNK_SIZE{256};

// The hw.cachelinesize value is 128 on Apple Silicon and 64 on x86_64. This
// constant uses the larger size so that counter groups stay on separate
// lines on both platforms.
// (std::hardware_destructive_interference_size is unavailable on libc++.)
inline constexpr std::size_t PIPELINE_CACHE_LINE_SIZE{128};

[[nodiscard]] constexpr std::size_t ParallelChunkSize(std::size_t recordCount,
                                                      std::size_t workerCount) noexcept
{
    const std::size_t targetChunks{std::max(workerCount, std::size_t{1}) *
                                   PARALLEL_CHUNKS_PER_WORKER};
    const std::size_t chunkSize{(recordCount / targetChunks) + ((recordCount % targetChunks) != 0)};
    return std::clamp(chunkSize, std::size_t{1}, MAX_PARALLEL_CHUNK_SIZE);
}

}  // namespace detail
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_SRC_PARALLELUTILS_HPP
