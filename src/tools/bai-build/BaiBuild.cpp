#include "BaiBuild.hpp"

#include <pbsamoa/index/BaiIndex.hpp>

#include <pbcopper/cli2/CLI.h>
#include <pbcopper/logging/Logging.h>

#include <filesystem>
#include <string>

#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace BaiBuild {
namespace {

// clang-format off
const PacBio::CLI_v2::PositionalArgument BaiBuildInput{
R"({
    "name" : "IN.bam",
    "description" : "Input BAM file to index.",
    "type" : "file"
})"};
// clang-format on

}  // namespace

PacBio::CLI_v2::Interface CreateInterface()
{
    PacBio::CLI_v2::Interface iface{"bai-build", "Build BAI index for a BAM file", "0.1.0"};
    iface.DisableLogFileOption();
    iface.DisableNumThreadsOption();
    iface.AddPositionalArgument(BaiBuildInput);
    return iface;
}

int Runner(const PacBio::CLI_v2::Results& results)
{
    const std::filesystem::path bamPath{results[BaiBuildInput]};
    const auto index{BaiIndex::Build(bamPath)};
    const std::filesystem::path outPath{bamPath.string() + ".bai"};
    index.ToFile(outPath);
    PBLOG_INFO << "Index written to " << outPath.string();
    return EXIT_SUCCESS;
}

}  // namespace BaiBuild
}  // namespace Samoa
}  // namespace PacBio
