#ifndef PBSAMOA_ZMIINTERNAL_HPP
#define PBSAMOA_ZMIINTERNAL_HPP

#include <array>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {
namespace detail {

inline constexpr std::size_t ZMI_HEADER_SIZE{64};
inline constexpr std::size_t ZMI_ENTRY_SIZE{16};  // rgId(4) + zmw(4) + virtualOffset(8)
inline constexpr std::uint16_t ZMI_ENTRY_SIZE_FIELD{16};

// Magic bytes: "ZMI\1"
inline constexpr std::array<std::byte, 4> ZMI_MAGIC{std::byte{'Z'}, std::byte{'M'}, std::byte{'I'},
                                                    std::byte{'\1'}};

// Version 1.0.0 encoded as 0x010000 (major=1, minor=0, patch=0)
inline constexpr std::uint32_t ZMI_VERSION{0x010000};

}  // namespace detail
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_ZMIINTERNAL_HPP
