#include "ZmiBuild.hpp"
#include "../../PathUtils.hpp"

#include <pbsamoa/PbSamoaLibraryInfo.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/ZmiBamWriter.hpp>

#include <pbcopper/cli2/Interface.h>
#include <pbcopper/cli2/PositionalArgument.h>
#include <pbcopper/cli2/Results.h>

#include <filesystem>
#include <print>
#include <stdexcept>
#include <string>
#include <vector>

#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace ZmiBuild {
namespace {

const CLI_v2::PositionalArgument Input{
    R"({
    "name" : "input",
    "description" : "Input BAM file.",
    "type" : "file"
})"};

const CLI_v2::PositionalArgument Output{
    R"({
    "name" : "output",
    "description" : "Output BAM file.",
    "type" : "file"
})"};

}  // namespace

CLI_v2::Interface CreateInterface()
{
    CLI_v2::Interface interface{"pbsamoa zmi-build", "Copy BAM and build ZMI index alongside.",
                                LibraryFormattedVersion()};
    // No custom thread flag; disable the built-in --num-threads to keep the
    // interface clean (this tool performs no threaded work).
    interface.DisableNumThreadsOption();
    interface.AddPositionalArguments({Input, Output});
    return interface;
}

int Runner(const CLI_v2::Results& results)
{
    // CLIv2 does not enforce the required positional-argument count, so guard it
    // here before indexing (matches the tool's prior "exactly two operands" rule).
    const std::vector<std::string>& positional{results.PositionalArguments()};
    if (positional.size() != 2) {
        throw std::runtime_error{"zmi-build requires exactly two arguments: <input> <output>"};
    }

    const std::filesystem::path inputPath{positional[0]};
    const std::filesystem::path outputPath{positional[1]};

    BamRawReader reader{inputPath};
    ZmiBamWriter writer{outputPath, reader.Header()};
    for (const auto& view : reader.Records()) {
        writer.Write(view);
    }

    const std::filesystem::path zmiPath{SidecarPath(outputPath, ".zmi")};
    std::println(stderr, "ZMI index written to {}", zmiPath.string());
    return EXIT_SUCCESS;
}

}  // namespace ZmiBuild
}  // namespace Samoa
}  // namespace PacBio
