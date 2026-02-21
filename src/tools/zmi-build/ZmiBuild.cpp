#include "ZmiBuild.hpp"

#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/ZmiBamWriter.hpp>

#include <filesystem>

#include <cstdio>
#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace ZmiBuild {

int Runner(int argc, char* argv[])
{
    if (argc < 2) {
        std::fprintf(stderr, "Usage: pbsamoa zmi-build IN.bam OUT.bam\n");
        return EXIT_FAILURE;
    }

    const std::filesystem::path inputPath{argv[0]};
    const std::filesystem::path outputPath{argv[1]};

    BamRawReader reader{inputPath};
    ZmiBamWriter writer{outputPath, reader.Header()};
    for (const auto& view : reader.Records()) {
        writer.Write(view);
    }

    std::fprintf(stderr, "ZMI index written to %s.zmi\n", outputPath.c_str());
    return EXIT_SUCCESS;
}

}  // namespace ZmiBuild
}  // namespace Samoa
}  // namespace PacBio
