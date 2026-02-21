#ifndef PBSAMOA_BINARYUTILS_HPP
#define PBSAMOA_BINARYUTILS_HPP

#include <algorithm>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief Read a little-endian uint16 from raw bytes.
inline std::uint16_t ReadU16LE(const std::byte* p)
{
    const std::uint16_t low{std::to_integer<std::uint16_t>(p[0])};
    const std::uint16_t high{std::to_integer<std::uint16_t>(p[1])};
    return low | (high << 8U);
}

/// \brief Read a little-endian uint32 from raw bytes.
inline std::uint32_t ReadU32LE(const std::byte* p)
{
    std::uint32_t v{};
    std::ranges::copy_n(p, sizeof(v), reinterpret_cast<std::byte*>(&v));
    return v;
}

/// \brief Read a little-endian int32 from raw bytes.
inline std::int32_t ReadI32LE(const std::byte* p)
{
    std::int32_t v{};
    std::ranges::copy_n(p, sizeof(v), reinterpret_cast<std::byte*>(&v));
    return v;
}

/// \brief Compute total BAM header size from binary data.
///
/// Returns the total number of bytes from BAM magic through end of reference
/// dictionary, or 0 if insufficient data.
inline std::size_t ComputeHeaderSize(const std::byte* data, std::size_t available)
{
    // Need at least: magic(4) + l_text(4) = 8
    if (available < 8) {
        return 0;
    }

    const std::uint32_t lText{ReadU32LE(data + 4)};
    // Need: magic(4) + l_text(4) + text(lText) + n_ref(4)
    const std::size_t minSize{8 + lText + 4};
    if (available < minSize) {
        return 0;
    }

    const std::uint32_t nRef{ReadU32LE(data + 8 + lText)};
    std::size_t offset{8 + lText + 4};

    for (std::uint32_t i{0}; i < nRef; ++i) {
        if (offset + 4 > available) {
            return 0;
        }
        const std::uint32_t lName{ReadU32LE(data + offset)};
        offset += 4;
        // name(lName) + l_ref(4)
        if (offset + lName + 4 > available) {
            return 0;
        }
        offset += lName + 4;
    }

    return offset;
}

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_BINARYUTILS_HPP
