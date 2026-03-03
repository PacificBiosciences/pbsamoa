#ifndef PBSAMOA_CORE_ENDIAN_HPP
#define PBSAMOA_CORE_ENDIAN_HPP

#include <bit>
#include <concepts>

#include <cstddef>
#include <cstdint>
#include <cstring>

namespace PacBio {
namespace Samoa {

static_assert(std::endian::native == std::endian::little, "Little-endian platform required");

/// \brief Read a little-endian integral value from raw bytes (unaligned-safe).
template <std::integral T>
inline T ReadLE(const std::byte* p)
{
    T v;
    std::memcpy(&v, p, sizeof(v));
    return v;
}

/// \brief Read a little-endian uint16 from raw bytes (unaligned-safe).
inline std::uint16_t ReadU16LE(const std::byte* p) { return ReadLE<std::uint16_t>(p); }

/// \brief Read a little-endian uint32 from raw bytes (unaligned-safe).
inline std::uint32_t ReadU32LE(const std::byte* p) { return ReadLE<std::uint32_t>(p); }

/// \brief Read a little-endian int32 from raw bytes (unaligned-safe).
inline std::int32_t ReadI32LE(const std::byte* p) { return ReadLE<std::int32_t>(p); }

/// \brief Read a little-endian uint64 from raw bytes (unaligned-safe).
inline std::uint64_t ReadU64LE(const std::byte* p) { return ReadLE<std::uint64_t>(p); }

/// \brief Read a little-endian int64 from raw bytes (unaligned-safe).
inline std::int64_t ReadI64LE(const std::byte* p) { return ReadLE<std::int64_t>(p); }

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_CORE_ENDIAN_HPP
