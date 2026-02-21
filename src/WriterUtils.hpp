#ifndef PBSAMOA_WRITERUTILS_HPP
#define PBSAMOA_WRITERUTILS_HPP

#include <atomic>
#include <chrono>
#include <filesystem>
#include <format>
#include <random>
#include <stdexcept>
#include <string_view>
#include <system_error>

#include <cstdint>

namespace PacBio {
namespace Samoa {
namespace detail {

/// \brief Generate a unique temporary path adjacent to \p finalPath.
///
/// The name is constructed from a nanosecond timestamp, a 64-bit RNG nonce,
/// and an atomic counter to guarantee uniqueness across threads and processes.
inline std::filesystem::path TemporaryWritePath(const std::filesystem::path& finalPath)
{
    static std::atomic<std::uint64_t> counter{0};
    thread_local std::mt19937_64 rng{std::random_device{}()};

    const std::uint64_t nowNs{
        static_cast<std::uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                                       std::chrono::system_clock::now().time_since_epoch())
                                       .count())};
    const std::uint64_t nonce{rng() ^ counter.fetch_add(1, std::memory_order_relaxed)};

    std::filesystem::path candidate{finalPath};
    candidate += std::format(".tmp.{}.{}", nowNs, nonce);
    return candidate;
}

/// \brief Resolve actual output path for a writer.
///
/// Returns \p finalPath directly unless temp-file writes are enabled, in which
/// case this returns a unique adjacent temporary path.
inline std::filesystem::path ResolveWritePath(const std::filesystem::path& finalPath,
                                              bool useTempFile)
{
    if (!useTempFile) {
        return finalPath;
    }
    return TemporaryWritePath(finalPath);
}

/// \brief Atomically rename \p src to \p dst, with fallback remove+retry.
///
/// On failure the temp file at \p src is removed (best-effort) and an
/// exception is thrown that includes \p writerName for diagnostics.
inline void AtomicRename(const std::filesystem::path& src, const std::filesystem::path& dst,
                         std::string_view writerName)
{
    std::error_code ec;
    std::filesystem::rename(src, dst, ec);
    if (ec) {
        std::error_code removeEc;
        std::filesystem::remove(dst, removeEc);
        ec.clear();
        std::filesystem::rename(src, dst, ec);
    }
    if (ec) {
        // Best-effort cleanup of the temp file.
        std::error_code cleanupEc;
        std::filesystem::remove(src, cleanupEc);

        throw std::runtime_error{std::format("{}: failed to rename temporary output {} to {} ({})",
                                             writerName, src.string(), dst.string(), ec.message())};
    }
}

}  // namespace detail
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_WRITERUTILS_HPP
