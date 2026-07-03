#include "BaiBuild.hpp"

#include "../../PathUtils.hpp"

#include <pbsamoa/PbSamoaLibraryInfo.hpp>
#include <pbsamoa/index/BaiIndex.hpp>

#include <pbcopper/cli2/Interface.h>
#include <pbcopper/cli2/Option.h>
#include <pbcopper/cli2/PositionalArgument.h>
#include <pbcopper/cli2/Results.h>

#include <algorithm>
#include <filesystem>
#include <print>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <cstdint>
#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace BaiBuild {
namespace {

const CLI_v2::Option Threads{
    R"({
    "names" : ["threads"],
    "description" : "BGZF inflate workers; 0 = auto (min of hardware concurrency and 8).",
    "type" : "unsigned integer",
    "default" : 0
})"};

const CLI_v2::PositionalArgument Input{
    R"({
    "name" : "input",
    "description" : "Input BAM file.",
    "type" : "file"
})"};

/// Resolve the worker count: 0 means auto = min(hardware_concurrency, 8).
std::size_t ResolveThreads(std::size_t requested)
{
    if (requested != 0) {
        return requested;
    }
    const std::uint32_t hardware{std::thread::hardware_concurrency()};
    const std::size_t available{(hardware == 0) ? std::size_t{1} : hardware};
    return std::min<std::size_t>(available, 8);
}

}  // namespace

CLI_v2::Interface CreateInterface()
{
    CLI_v2::Interface interface{"pbsamoa bai-build", "Build BAI index for a BAM file.",
                                LibraryFormattedVersion()};
    // `bai-build` keeps its own --threads (0 = auto, capped at min(hw, 8)); the built-in
    // --num-threads resolves 0 to the raw hardware count, which would drop that cap.
    interface.DisableNumThreadsOption();
    interface.AddOptions({Threads});
    interface.AddPositionalArguments({Input});
    return interface;
}

int Runner(const CLI_v2::Results& results)
{
    // CLIv2 does not enforce the required positional-argument count, so guard it
    // here before indexing (matches the tool's prior "exactly one operand" rule).
    const std::vector<std::string>& positional{results.PositionalArguments()};
    if (positional.size() != 1) {
        throw std::runtime_error{"bai-build requires exactly one argument: <input>"};
    }

    // Read as fixed-width type; std::size_t has no exact conversion operator.
    const std::uint32_t requestedThreads{results[Threads]};
    const std::size_t numWorkers{ResolveThreads(requestedThreads)};

    const std::filesystem::path bamPath{positional[0]};
    const BaiIndex index{BaiIndex::Build(bamPath, numWorkers)};
    const std::filesystem::path outPath{SidecarPath(bamPath, ".bai")};
    index.ToFile(outPath);
    std::println(stderr, "Index written to {}", outPath.string());
    return EXIT_SUCCESS;
}

}  // namespace BaiBuild
}  // namespace Samoa
}  // namespace PacBio
