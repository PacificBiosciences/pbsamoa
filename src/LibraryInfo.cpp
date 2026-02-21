#include <pbsamoa/PbSamoaLibraryInfo.hpp>

#include "LibraryGitHash.hpp"
#include "LibraryVersion.hpp"

#include <pbcopper/LibraryInfo.h>

namespace PacBio {
namespace Samoa {

Library::Info LibraryInfo()
{
    return {
        .Name = "pbsamoa",
        .Release = std::string{RELEASE_VERSION},
        .GitSha1 = std::string{LIBRARY_GIT_SHA1},
    };
}

std::string LibraryFormattedVersion()
{
    const Library::Info info{LibraryInfo()};
    return info.Release + " (commit " + info.GitSha1 + ')';
}

}  // namespace Samoa
}  // namespace PacBio
