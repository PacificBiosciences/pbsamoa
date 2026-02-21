#include <pbsamoa/PbSamoaLibraryInfo.hpp>

#include "LibraryGitHash.hpp"
#include "LibraryVersion.hpp"

#include <string>

namespace PacBio {
namespace Samoa {

LibraryInfo GetLibraryInfo()
{
    return {
        .Name = "pbsamoa",
        .Release = RELEASE_VERSION,
        .GitSha1 = LIBRARY_GIT_SHA1,
    };
}

std::string LibraryFormattedVersion()
{
    const LibraryInfo info{GetLibraryInfo()};
    return std::string{info.Release} + " (commit " + std::string{info.GitSha1} + ')';
}

}  // namespace Samoa
}  // namespace PacBio
