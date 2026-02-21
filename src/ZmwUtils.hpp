#ifndef PBSAMOA_ZMWUTILS_HPP
#define PBSAMOA_ZMWUTILS_HPP

#include <string_view>

#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief Parse ZMW hole number from PacBio read name format "movie/<zmw>/...".
/// \returns hole number, or 0 on parse failure.
std::int32_t ParseZmwFromName(std::string_view name);

/// \brief Extract the "movie/zmw" prefix from a PacBio read name.
/// \returns "movie/zmw" substring, or the full name if no slashes found.
std::string_view MovieZmwPrefix(std::string_view name);

/// \brief Extract the base read-group ID before any barcode suffix.
std::string_view ReadGroupBaseId(std::string_view readGroupId);

/// \brief Parse PacBio read-group ID string to its numeric int32 form.
/// \returns numeric ID, or 0 on parse failure.
std::int32_t ParseReadGroupId(std::string_view readGroupId);

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_ZMWUTILS_HPP
