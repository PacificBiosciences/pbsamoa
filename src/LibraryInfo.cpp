#include <pbsamoa/PbSamoaLibraryInfo.hpp>

#include "LibraryGitHash.hpp"
#include "LibraryVersion.hpp"

#include <format>

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
    const LibraryInfo libraryInfo{GetLibraryInfo()};
    return std::format("{} (commit {})", libraryInfo.Release, libraryInfo.GitSha1);
}

}  // namespace Samoa
}  // namespace PacBio
