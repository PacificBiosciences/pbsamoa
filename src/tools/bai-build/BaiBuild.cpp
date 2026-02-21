#include "BaiBuild.hpp"

#include <pbsamoa/index/BaiIndex.hpp>

#include <filesystem>
#include <print>

#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace BaiBuild {

int Runner(int argc, char** argv)
{
    if (argc < 1) {
        std::println(stderr, "Usage: pbsamoa bai-build IN.bam");
        return EXIT_FAILURE;
    }

    const std::filesystem::path bamPath{argv[0]};
    const BaiIndex index{BaiIndex::Build(bamPath)};
    const std::filesystem::path outPath{bamPath.string() + ".bai"};
    index.ToFile(outPath);
    std::println(stderr, "Index written to {}", outPath.string());
    return EXIT_SUCCESS;
}

}  // namespace BaiBuild
}  // namespace Samoa
}  // namespace PacBio
