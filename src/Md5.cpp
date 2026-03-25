// RFC 1321 MD5 implementation for PacBio read group ID generation.

#include "Md5.hpp"

#include <array>
#include <span>
#include <string>
#include <string_view>

#include <cstddef>
#include <cstdint>
#include <cstring>

namespace PacBio {
namespace Samoa {
namespace detail {
namespace {

// clang-format off
constexpr std::array<std::uint32_t, 64> K{
    0xd76aa478, 0xe8c7b756, 0x242070db, 0xc1bdceee,
    0xf57c0faf, 0x4787c62a, 0xa8304613, 0xfd469501,
    0x698098d8, 0x8b44f7af, 0xffff5bb1, 0x895cd7be,
    0x6b901122, 0xfd987193, 0xa679438e, 0x49b40821,
    0xf61e2562, 0xc040b340, 0x265e5a51, 0xe9b6c7aa,
    0xd62f105d, 0x02441453, 0xd8a1e681, 0xe7d3fbc8,
    0x21e1cde6, 0xc33707d6, 0xf4d50d87, 0x455a14ed,
    0xa9e3e905, 0xfcefa3f8, 0x676f02d9, 0x8d2a4c8a,
    0xfffa3942, 0x8771f681, 0x6d9d6122, 0xfde5380c,
    0xa4beea44, 0x4bdecfa9, 0xf6bb4b60, 0xbebfbc70,
    0x289b7ec6, 0xeaa127fa, 0xd4ef3085, 0x04881d05,
    0xd9d4d039, 0xe6db99e5, 0x1fa27cf8, 0xc4ac5665,
    0xf4292244, 0x432aff97, 0xab9423a7, 0xfc93a039,
    0x655b59c3, 0x8f0ccc92, 0xffeff47d, 0x85845dd1,
    0x6fa87e4f, 0xfe2ce6e0, 0xa3014314, 0x4e0811a1,
    0xf7537e82, 0xbd3af235, 0x2ad7d2bb, 0xeb86d391,
};

constexpr std::array<std::uint32_t, 64> S{
    7, 12, 17, 22,  7, 12, 17, 22,  7, 12, 17, 22,  7, 12, 17, 22,
    5,  9, 14, 20,  5,  9, 14, 20,  5,  9, 14, 20,  5,  9, 14, 20,
    4, 11, 16, 23,  4, 11, 16, 23,  4, 11, 16, 23,  4, 11, 16, 23,
    6, 10, 15, 21,  6, 10, 15, 21,  6, 10, 15, 21,  6, 10, 15, 21,
};
// clang-format on

constexpr std::uint32_t LeftRotate(std::uint32_t x, std::uint32_t c)
{
    return (x << c) | (x >> (32 - c));
}

void Md5Transform(std::array<std::uint32_t, 4>& state, std::span<const std::uint8_t, 64> block)
{
    std::array<std::uint32_t, 16> M{};
    for (std::size_t i{0}; i < 16; ++i) {
        std::memcpy(&M[i], std::data(block) + i * 4, 4);
    }

    std::uint32_t a{state[0]};
    std::uint32_t b{state[1]};
    std::uint32_t c{state[2]};
    std::uint32_t d{state[3]};

    for (std::uint32_t i{0}; i < 64; ++i) {
        std::uint32_t f{0};
        std::uint32_t g{0};
        if (i < 16) {
            f = (b & c) | (~b & d);
            g = i;
        } else if (i < 32) {
            f = (d & b) | (~d & c);
            g = (5 * i + 1) % 16;
        } else if (i < 48) {
            f = b ^ c ^ d;
            g = (3 * i + 5) % 16;
        } else {
            f = c ^ (b | ~d);
            g = (7 * i) % 16;
        }
        f = f + a + K[i] + M[g];
        a = d;
        d = c;
        c = b;
        b = b + LeftRotate(f, S[i]);
    }

    state[0] += a;
    state[1] += b;
    state[2] += c;
    state[3] += d;
}

}  // namespace

std::string Md5DigestToHex(std::span<const std::byte, 16> digest)
{
    static constexpr std::string_view HEX{"0123456789abcdef"};

    std::string hex;
    hex.resize(32);
    for (std::size_t i{0}; i < 16; ++i) {
        const std::uint8_t value{std::to_integer<std::uint8_t>(digest[i])};
        hex[2 * i] = HEX[value >> 4];
        hex[2 * i + 1] = HEX[value & 0x0F];
    }
    return hex;
}

std::array<std::byte, 16> ComputeMd5(std::span<const std::byte> data)
{
    std::array<std::uint32_t, 4> state{0x67452301, 0xefcdab89, 0x98badcfe, 0x10325476};

    const std::uint64_t bitLen{static_cast<std::uint64_t>(std::size(data)) * 8};

    // Process complete 64-byte blocks.
    const std::span<const std::uint8_t> bytes{reinterpret_cast<const std::uint8_t*>(data.data()),
                                              std::size(data)};
    std::size_t offset{0};
    while ((offset + 64) <= std::size(bytes)) {
        Md5Transform(state, bytes.subspan(offset).first<64>());
        offset += 64;
    }

    // Pad: remaining bytes + 0x80 + zeros + 8-byte length.
    std::array<std::uint8_t, 128> buffer{};
    const std::size_t remaining{std::size(bytes) - offset};
    std::memcpy(std::data(buffer), std::data(bytes) + offset, remaining);
    buffer[remaining] = 0x80;

    std::size_t padded{64U};
    if (remaining >= 56) {
        padded = 128U;
    }
    std::memcpy(std::data(buffer) + padded - 8, &bitLen, 8);

    const std::span<const std::uint8_t, 128> padSpan{buffer};
    Md5Transform(state, padSpan.first<64>());
    if (padded == 128U) {
        Md5Transform(state, padSpan.last<64>());
    }

    std::array<std::byte, 16> md5{};
    for (std::size_t i{0}; i < std::size(state); ++i) {
        for (std::size_t j{0}; j < 4; ++j) {
            md5[i * 4 + j] = static_cast<std::byte>((state[i] >> (8 * j)) & 0xFFu);
        }
    }
    return md5;
}

std::string Md5Hex(std::span<const std::byte> data) { return Md5DigestToHex(ComputeMd5(data)); }

std::string Md5Hex(std::string_view text) { return Md5Hex(std::as_bytes(std::span{text})); }

}  // namespace detail
}  // namespace Samoa
}  // namespace PacBio
