#ifndef PBSAMOA_CORE_ENDIAN_HPP
#define PBSAMOA_CORE_ENDIAN_HPP

#include <bit>

#include <cstddef>
#include <cstdint>
#include <cstring>

namespace PacBio {
namespace Samoa {

static_assert(std::endian::native == std::endian::little, "Little-endian platform required");

/// \brief Read a little-endian uint16 from raw bytes (unaligned-safe).
inline std::uint16_t ReadU16LE(const std::byte* p)
{
    std::uint16_t v;
    std::memcpy(&v, p, sizeof(v));
    return v;
}

/// \brief Read a little-endian uint32 from raw bytes (unaligned-safe).
inline std::uint32_t ReadU32LE(const std::byte* p)
{
    std::uint32_t v;
    std::memcpy(&v, p, sizeof(v));
    return v;
}

/// \brief Read a little-endian int32 from raw bytes (unaligned-safe).
inline std::int32_t ReadI32LE(const std::byte* p)
{
    return std::bit_cast<std::int32_t>(ReadU32LE(p));
}

/// \brief Read a little-endian uint64 from raw bytes (unaligned-safe).
inline std::uint64_t ReadU64LE(const std::byte* p)
{
    std::uint64_t v;
    std::memcpy(&v, p, sizeof(v));
    return v;
}

/// \brief Read a little-endian int64 from raw bytes (unaligned-safe).
inline std::int64_t ReadI64LE(const std::byte* p)
{
    return std::bit_cast<std::int64_t>(ReadU64LE(p));
}

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_CORE_ENDIAN_HPP
