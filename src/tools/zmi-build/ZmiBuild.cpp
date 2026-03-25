#include "ZmiBuild.hpp"
#include "../../PathUtils.hpp"

#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/ZmiBamWriter.hpp>

#include <filesystem>
#include <print>

#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace ZmiBuild {

int Runner(int argc, char** argv)
{
    if (argc < 2) {
        std::println(stderr, "Usage: pbsamoa zmi-build IN.bam OUT.bam");
        return EXIT_FAILURE;
    }

    const std::filesystem::path inputPath{argv[0]};
    const std::filesystem::path outputPath{argv[1]};

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
