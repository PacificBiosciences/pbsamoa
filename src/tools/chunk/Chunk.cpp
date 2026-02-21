#include "Chunk.hpp"
#include "../ParseUtils.hpp"

#include "../SamOutput.hpp"

#include <pbsamoa/io/BamRawReader.hpp>

#include <filesystem>
#include <print>

#include <cstdint>
#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace ChunkTool {

int Runner(int argc, char** argv)
{
    if (argc < 3) {
        std::println(stderr, "Usage: pbsamoa chunk IN.bam CHUNK TOTAL");
        return EXIT_FAILURE;
    }

    const std::filesystem::path bamPath{argv[0]};
    const BamRawReaderConfig readerConfig{
        .ChunkNum = Tools::ParseIntegerOrThrow<std::int32_t>(argv[1], "CHUNK"),
        .TotalChunks = Tools::ParseIntegerOrThrow<std::int32_t>(argv[2], "TOTAL"),
    };
    BamRawReader reader{bamPath, readerConfig};
    const auto& header{reader.Header()};
    std::print("{}", header.ToText());

    while (const auto view = reader.ReadRecord()) {
        WriteViewAsSam(header, *view);
    }

    return EXIT_SUCCESS;
}

}  // namespace ChunkTool
}  // namespace Samoa
}  // namespace PacBio
