#ifndef PBSAMOA_PATHUTILS_HPP
#define PBSAMOA_PATHUTILS_HPP

#include <filesystem>
#include <string>
#include <string_view>

namespace PacBio {
namespace Samoa {

inline std::filesystem::path SidecarPath(const std::filesystem::path& path, std::string_view suffix)
{
    std::filesystem::path sidecarPath{path};
    sidecarPath += std::string{suffix};
    return sidecarPath;
}

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_PATHUTILS_HPP
