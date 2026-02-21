#ifndef PBSAMOA_CRAM_MD5_HPP
#define PBSAMOA_CRAM_MD5_HPP

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace PacBio {
namespace Samoa {

inline constexpr std::array<std::uint32_t, 64> MD5_SHIFT_AMOUNTS{
    7,  12, 17, 22, 7,  12, 17, 22, 7,  12, 17, 22, 7,  12, 17, 22, 5,  9,  14, 20, 5,  9,
    14, 20, 5,  9,  14, 20, 5,  9,  14, 20, 4,  11, 16, 23, 4,  11, 16, 23, 4,  11, 16, 23,
    4,  11, 16, 23, 6,  10, 15, 21, 6,  10, 15, 21, 6,  10, 15, 21, 6,  10, 15, 21};

inline constexpr std::array<std::uint32_t, 64> MD5_TABLE{
    0xd76aa478u, 0xe8c7b756u, 0x242070dbu, 0xc1bdceeeu, 0xf57c0fafu, 0x4787c62au, 0xa8304613u,
    0xfd469501u, 0x698098d8u, 0x8b44f7afu, 0xffff5bb1u, 0x895cd7beu, 0x6b901122u, 0xfd987193u,
    0xa679438eu, 0x49b40821u, 0xf61e2562u, 0xc040b340u, 0x265e5a51u, 0xe9b6c7aau, 0xd62f105du,
    0x02441453u, 0xd8a1e681u, 0xe7d3fbc8u, 0x21e1cde6u, 0xc33707d6u, 0xf4d50d87u, 0x455a14edu,
    0xa9e3e905u, 0xfcefa3f8u, 0x676f02d9u, 0x8d2a4c8au, 0xfffa3942u, 0x8771f681u, 0x6d9d6122u,
    0xfde5380cu, 0xa4beea44u, 0x4bdecfa9u, 0xf6bb4b60u, 0xbebfbc70u, 0x289b7ec6u, 0xeaa127fau,
    0xd4ef3085u, 0x04881d05u, 0xd9d4d039u, 0xe6db99e5u, 0x1fa27cf8u, 0xc4ac5665u, 0xf4292244u,
    0x432aff97u, 0xab9423a7u, 0xfc93a039u, 0x655b59c3u, 0x8f0ccc92u, 0xffeff47du, 0x85845dd1u,
    0x6fa87e4fu, 0xfe2ce6e0u, 0xa3014314u, 0x4e0811a1u, 0xf7537e82u, 0xbd3af235u, 0x2ad7d2bbu,
    0xeb86d391u};

inline constexpr std::uint32_t LeftRotate(std::uint32_t value, std::uint32_t amount)
{
    return std::rotl(value, static_cast<int>(amount));
}

inline constexpr std::uint32_t Md5ByteValue(const std::byte value)
{
    return static_cast<std::uint32_t>(std::to_integer<std::uint8_t>(value));
}

inline std::uint32_t LoadMd5Word(std::span<const std::byte> bytes, std::size_t offset)
{
    return Md5ByteValue(bytes[offset]) | (Md5ByteValue(bytes[offset + 1]) << 8U) |
           (Md5ByteValue(bytes[offset + 2]) << 16U) | (Md5ByteValue(bytes[offset + 3]) << 24U);
}

inline std::array<std::byte, 16> ComputeMd5(std::span<const std::byte> data)
{
    std::vector<std::byte> padded;
    padded.reserve(std::size(data) + 72);
    padded.insert(std::end(padded), std::begin(data), std::end(data));
    padded.push_back(std::byte{0x80});

    while ((std::size(padded) % 64) != 56) {
        padded.push_back(std::byte{0});
    }

    const std::uint64_t bitLength{std::uint64_t{std::size(data)} * 8U};
    for (std::size_t i{0}; i < 8; ++i) {
        const std::uint64_t shifted{(bitLength >> (8 * i)) & 0xFFU};
        padded.push_back(static_cast<std::byte>(shifted));
    }

    std::uint32_t a0{0x67452301U};
    std::uint32_t b0{0xefcdab89U};
    std::uint32_t c0{0x98badcfeU};
    std::uint32_t d0{0x10325476U};

    for (std::size_t chunkOffset{0}; chunkOffset < std::size(padded); chunkOffset += 64) {
        std::array<std::uint32_t, 16> words{};
        for (std::size_t i{0}; i < 16; ++i) {
            const std::size_t base{chunkOffset + (i * 4)};
            words[i] = LoadMd5Word(padded, base);
        }

        std::uint32_t a = a0;
        std::uint32_t b = b0;
        std::uint32_t c = c0;
        std::uint32_t d = d0;

        for (std::size_t i = 0; i < 64; ++i) {
            std::uint32_t f = 0;
            std::size_t g = 0;

            if (i < 16) {
                f = (b & c) | ((~b) & d);
                g = i;
            } else if (i < 32) {
                f = (d & b) | ((~d) & c);
                g = (5 * i + 1) % 16;
            } else if (i < 48) {
                f = b ^ c ^ d;
                g = (3 * i + 5) % 16;
            } else {
                f = c ^ (b | (~d));
                g = (7 * i) % 16;
            }

            const std::uint32_t temp = d;
            d = c;
            c = b;
            const std::uint32_t sum = a + f + MD5_TABLE[i] + words[g];
            b = b + LeftRotate(sum, MD5_SHIFT_AMOUNTS[i]);
            a = temp;
        }

        a0 += a;
        b0 += b;
        c0 += c;
        d0 += d;
    }

    std::array<std::byte, 16> digest{};
    const std::array<std::uint32_t, 4> state{a0, b0, c0, d0};
    for (std::size_t i{0}; i < std::size(state); ++i) {
        for (std::size_t j{0}; j < 4; ++j) {
            digest[i * 4 + j] = static_cast<std::byte>((state[i] >> (8 * j)) & 0xFFu);
        }
    }
    return digest;
}

inline std::string Md5ToHex(std::span<const std::byte, 16> digest)
{
    static constexpr std::array<char, 16> HEX{'0', '1', '2', '3', '4', '5', '6', '7',
                                              '8', '9', 'a', 'b', 'c', 'd', 'e', 'f'};

    std::string hex;
    hex.resize(32);
    for (std::size_t i{0}; i < 16; ++i) {
        const std::uint8_t value{std::to_integer<std::uint8_t>(digest[i])};
        hex[2 * i] = HEX[value >> 4];
        hex[2 * i + 1] = HEX[value & 0x0F];
    }
    return hex;
}

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_CRAM_MD5_HPP
