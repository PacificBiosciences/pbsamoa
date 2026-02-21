#ifndef PBSAMOA_PBSAMOALIBRARYINFO_HPP
#define PBSAMOA_PBSAMOALIBRARYINFO_HPP

#include <string>
#include <string_view>

namespace PacBio {
namespace Samoa {

struct LibraryInfo
{
    std::string_view Name;
    std::string_view Release;
    std::string_view GitSha1;
};

LibraryInfo GetLibraryInfo();
std::string LibraryFormattedVersion();

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_PBSAMOALIBRARYINFO_HPP
