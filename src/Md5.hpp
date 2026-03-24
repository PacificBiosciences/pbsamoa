#ifndef PBSAMOA_MD5_HPP
#define PBSAMOA_MD5_HPP

#include <span>
#include <string>
#include <string_view>

#include <cstddef>

namespace PacBio {
namespace Samoa {
namespace detail {

/// \brief Compute the MD5 hex digest of input bytes.
///
/// Returns a 32-character lowercase hex string (RFC 1321).
std::string Md5Hex(std::span<const std::byte> data);

/// \brief Convenience overload for string data.
std::string Md5Hex(std::string_view text);

}  // namespace detail
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_MD5_HPP
