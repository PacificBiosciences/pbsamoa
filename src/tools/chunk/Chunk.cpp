#include "Chunk.hpp"

#include "../SamOutput.hpp"

#include <pbsamoa/io/BamRawReader.hpp>

#include <pbcopper/cli2/CLI.h>

#include <filesystem>
#include <string>

#include <cstdint>
#include <cstdio>
#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace ChunkTool {
namespace {

// clang-format off
const PacBio::CLI_v2::PositionalArgument ChunkInput{
R"({
    "name" : "IN.bam",
    "description" : "Input BAM file (requires .zmi index alongside).",
    "type" : "file"
})"};

const PacBio::CLI_v2::PositionalArgument ChunkNum{
R"({
    "name" : "CHUNK",
    "description" : "Chunk number (1-based).",
    "type" : "integer"
})"};

const PacBio::CLI_v2::PositionalArgument ChunkTotal{
R"({
    "name" : "TOTAL",
    "description" : "Total number of chunks.",
    "type" : "integer"
})"};
// clang-format on

}  // namespace

PacBio::CLI_v2::Interface CreateInterface()
{
    PacBio::CLI_v2::Interface iface{"chunk", "Dump a chunk of BAM records as SAM text", "0.1.0"};
    iface.DisableLogFileOption();
    iface.DisableNumThreadsOption();
    iface.AddPositionalArguments({ChunkInput, ChunkNum, ChunkTotal});
    return iface;
}

int Runner(const PacBio::CLI_v2::Results& results)
{
    const std::filesystem::path bamPath{results[ChunkInput]};
    const std::int32_t chunkNum{std::stoi(std::string{results[ChunkNum]})};
    const std::int32_t totalChunks{std::stoi(std::string{results[ChunkTotal]})};

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
