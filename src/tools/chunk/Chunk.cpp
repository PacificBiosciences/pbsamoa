#include "Chunk.hpp"

#include "../SamOutput.hpp"

#include <pbsamoa/io/BamRawReader.hpp>

#include <charconv>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>

#include <cstdint>
#include <cstdio>
#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace ChunkTool {

int Runner(int argc, char* argv[])
{
    if (argc < 3) {
        std::fprintf(stderr, "Usage: pbsamoa chunk IN.bam CHUNK TOTAL\n");
        return EXIT_FAILURE;
    }

    const std::filesystem::path bamPath{argv[0]};

    std::int32_t chunkNum{0};
    std::int32_t totalChunks{0};

    {
        const std::string_view chunkStr{argv[1]};
        auto r1{std::from_chars(std::data(chunkStr), std::data(chunkStr) + std::size(chunkStr),
                                chunkNum)};
        if (r1.ec != std::errc{}) {
            throw std::runtime_error{"invalid CHUNK: " + std::string{chunkStr}};
        }

        const std::string_view totalStr{argv[2]};
        auto r2{std::from_chars(std::data(totalStr), std::data(totalStr) + std::size(totalStr),
                                totalChunks)};
        if (r2.ec != std::errc{}) {
            throw std::runtime_error{"invalid TOTAL: " + std::string{totalStr}};
        }
    }

    BamRawReader reader{bamPath,
                        BamRawReaderConfig{.ChunkNum = chunkNum, .TotalChunks = totalChunks}};
    const auto& header{reader.Header()};
    std::fputs(header.ToText().c_str(), stdout);

    while (const auto view = reader.ReadRecord()) {
        WriteViewAsSam(header, *view);
    }

    return EXIT_SUCCESS;
}

}  // namespace ChunkTool
}  // namespace Samoa
}  // namespace PacBio
