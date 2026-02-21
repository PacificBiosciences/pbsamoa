#include <pbsamoa/PbSamoaLibraryInfo.hpp>

#include <gtest/gtest.h>

namespace PacBio {
namespace Samoa {

TEST(Smoke, LibraryInfoDefined)
{
    const auto info = LibraryInfo();
    EXPECT_EQ(info.Name, "pbsamoa");
    EXPECT_FALSE(info.Release.empty());
    EXPECT_FALSE(info.GitSha1.empty());
}

}  // namespace Samoa
}  // namespace PacBio
