#ifndef PBSAMOA_ZMWUTILS_HPP
#define PBSAMOA_ZMWUTILS_HPP

#include <pbsamoa/index/ZmwIndex.hpp>

#include <string_view>

#include <cstdint>

namespace PacBio {
namespace Samoa {

class TagMap;

/// \brief Parse ZMW hole number from PacBio read name format "movie/<zmw>/...".
/// \returns hole number, or 0 on parse failure.
std::int32_t ParseZmwFromName(std::string_view name);

/// \brief Parse PacBio read-group ID string to its numeric int32 form.
/// \returns numeric ID, or 0 on parse failure.
std::int32_t ParseReadGroupId(std::string_view readGroupId);

/// \brief Parse the ZMW identity carried by a record name and tags.
/// \returns `(rgId, zmw)` with either field set to 0 when parsing fails.
ZmwIdentity ParseZmwIdentity(std::string_view name, const TagMap& tags);

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_ZMWUTILS_HPP
