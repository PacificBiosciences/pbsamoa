#ifndef PBSAMOA_PBSAMOALIBRARYINFO_HPP
#define PBSAMOA_PBSAMOALIBRARYINFO_HPP

#include <pbcopper/LibraryInfo.h>

#include <string>

namespace PacBio {
namespace Samoa {

Library::Info LibraryInfo();
std::string LibraryFormattedVersion();

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_PBSAMOALIBRARYINFO_HPP
