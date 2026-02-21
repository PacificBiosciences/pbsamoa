#include "BaiBuild.hpp"

#include <pbsamoa/index/BaiIndex.hpp>

#include <filesystem>
#include <string_view>

#include <cstdio>
#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace BaiBuild {

int Runner(int argc, char* argv[])
{
    if (argc < 1) {
        std::fprintf(stderr, "Usage: pbsamoa bai-build IN.bam\n");
        return EXIT_FAILURE;
    }

    const std::filesystem::path bamPath{argv[0]};
    const auto index{BaiIndex::Build(bamPath)};
    const std::filesystem::path outPath{bamPath.string() + ".bai"};
    index.ToFile(outPath);
    std::fprintf(stderr, "Index written to %s\n", outPath.c_str());
    return EXIT_SUCCESS;
}

}  // namespace BaiBuild
}  // namespace Samoa
}  // namespace PacBio
