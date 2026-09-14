#ifndef PBSAMOA_ZMIINTERNAL_HPP
#define PBSAMOA_ZMIINTERNAL_HPP

#include <array>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {
namespace detail {

inline constexpr std::size_t ZMI_HEADER_SIZE{64};
inline constexpr std::size_t ZMI_BASE_ENTRY_SIZE{16};
inline constexpr std::size_t ZMI_ENTRY_SIZE{24};  // base(16) + forward(4) + reverse(4)
inline constexpr std::uint16_t ZMI_ENTRY_SIZE_FIELD{24};

// Magic bytes: "ZMI\1"
inline constexpr std::array<std::byte, 4> ZMI_MAGIC{std::byte{'Z'}, std::byte{'M'}, std::byte{'I'},
                                                    std::byte{'\1'}};

// Version 1.1.0 adds barcode IDs after the backward-compatible 16-byte prefix.
inline constexpr std::uint32_t ZMI_VERSION{0x010100};

}  // namespace detail
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_ZMIINTERNAL_HPP
