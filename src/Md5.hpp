#ifndef PBSAMOA_MD5_HPP
#define PBSAMOA_MD5_HPP

#include <array>
#include <functional>
#include <span>
#include <string>
#include <string_view>

#include <cstddef>

namespace PacBio {
namespace Samoa {
namespace detail {

/// \brief Compute the raw 16-byte MD5 digest of input bytes.
std::array<std::byte, 16> ComputeMd5(std::span<const std::byte> data);

/// \brief Format a raw 16-byte MD5 digest as lowercase hex.
std::string Md5DigestToHex(std::span<const std::byte, 16> digest);

/// \brief Compute the MD5 hex digest of input bytes.
///
/// Returns a 32-character lowercase hex string (RFC 1321).
std::string Md5Hex(std::span<const std::byte> data);

/// \brief Convenience overload for string data.
std::string Md5Hex(std::string_view text);

/// \brief Compute the raw 16-byte MD5 digest of the bases in upper case, without a copy of
///        them (the M5 of a CRAM reference and the reference MD5 of a CRAM slice).
std::array<std::byte, 16> ComputeUppercaseMd5(std::string_view bases);

/// \brief Compute the raw 16-byte MD5 digest of bases in upper case that come in parts:
///        nextPart gives the next part, an empty view after the last.
std::array<std::byte, 16> ComputeUppercaseMd5(const std::function<std::string_view()>& nextPart);

}  // namespace detail
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_MD5_HPP
