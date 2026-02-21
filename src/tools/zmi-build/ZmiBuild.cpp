#include "ZmiBuild.hpp"

#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/ZmiBamWriter.hpp>

#include <pbcopper/cli2/CLI.h>
#include <pbcopper/logging/Logging.h>

#include <filesystem>
#include <string>

#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace ZmiBuild {
namespace {

// clang-format off
const PacBio::CLI_v2::PositionalArgument ZmiBuildInput{
R"({
    "name" : "IN.bam",
    "description" : "Input BAM file.",
    "type" : "file"
})"};

const PacBio::CLI_v2::PositionalArgument ZmiBuildOutput{
R"({
    "name" : "OUT.bam",
    "description" : "Output BAM file (ZMI created alongside as OUT.bam.zmi).",
    "type" : "file"
})"};
// clang-format on

}  // namespace

PacBio::CLI_v2::Interface CreateInterface()
{
    PacBio::CLI_v2::Interface iface{"zmi-build", "Copy BAM and build ZMI index alongside", "0.1.0"};
    iface.DisableLogFileOption();
    iface.DisableNumThreadsOption();
    iface.AddPositionalArguments({ZmiBuildInput, ZmiBuildOutput});
    return iface;
}

int Runner(const PacBio::CLI_v2::Results& results)
{
    const std::filesystem::path inputPath{results[ZmiBuildInput]};
    const std::filesystem::path outputPath{results[ZmiBuildOutput]};

    BamRawReader reader{inputPath};
    ZmiBamWriter writer{outputPath, reader.Header()};
    for (const auto& view : reader.Records()) {
        writer.Write(view);
    }

    PBLOG_INFO << "ZMI index written to " << outputPath.string() + ".zmi";
    return EXIT_SUCCESS;
}

}  // namespace ZmiBuild
}  // namespace Samoa
}  // namespace PacBio
