#ifndef PBSAMOA_CRAM_MD5_HPP
#define PBSAMOA_CRAM_MD5_HPP

#include "Md5.hpp"

#include <array>
#include <cstddef>
#include <span>
#include <string>

namespace PacBio {
namespace Samoa {

inline std::array<std::byte, 16> ComputeMd5(std::span<const std::byte> data)
{
    return detail::ComputeMd5(data);
}

inline std::string Md5ToHex(std::span<const std::byte, 16> digest)
{
    return detail::Md5DigestToHex(digest);
}

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_CRAM_MD5_HPP
