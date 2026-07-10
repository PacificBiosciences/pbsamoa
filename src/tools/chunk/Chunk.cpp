#include "Chunk.hpp"

#include "../ParseUtils.hpp"
#include "../SamOutput.hpp"

#include <pbsamoa/PbSamoaLibraryInfo.hpp>
#include <pbsamoa/io/BamRawReader.hpp>

#include <pbcopper/cli2/Interface.h>
#include <pbcopper/cli2/Option.h>
#include <pbcopper/cli2/PositionalArgument.h>
#include <pbcopper/cli2/Results.h>

#include <filesystem>
#include <print>
#include <stdexcept>
#include <string>
#include <vector>

#include <cstdint>
#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace ChunkTool {
namespace {

const CLI_v2::Option Mode{
    R"({
    "names" : ["mode"],
    "description" : "Chunking mode: contiguous (default) selects one contiguous ZMW range per chunk; scatter uses a seeded balanced shuffle of tile-ZMW groups sampled across the whole file (chaotic but deterministic).",
    "type" : "string",
    "choices" : ["contiguous", "scatter"],
    "default" : "contiguous"
})"};

const CLI_v2::Option Tile{
    R"({
    "names" : ["tile"],
    "description" : "Scatter only: read up to M consecutive ZMWs per seek (>= 1).",
    "type" : "integer",
    "default" : 100
})"};

// ChunkSeed is std::uint64_t, so parse as a string via ParseIntegerOrThrow<uint64_t>:
// CLIv2's "unsigned integer" type is only 32-bit and would truncate the seed range.
const CLI_v2::Option Seed{
    R"({
    "names" : ["seed"],
    "description" : "Scatter only: deterministic shuffle seed.",
    "type" : "string",
    "default" : "42"
})"};

const CLI_v2::PositionalArgument Input{
    R"({
    "name" : "input",
    "description" : "Input BAM file.",
    "type" : "file"
})"};

const CLI_v2::PositionalArgument ChunkNum{
    R"({
    "name" : "chunk",
    "description" : "1-based chunk index."
})"};

const CLI_v2::PositionalArgument Total{
    R"({
    "name" : "total",
    "description" : "Total number of chunks."
})"};

}  // namespace

CLI_v2::Interface CreateInterface()
{
    CLI_v2::Interface interface{"pbsamoa chunk",
                                "Dump a chunk of BAM records as SAM text on stdout.",
                                LibraryFormattedVersion()};
    // Chunking is process-level parallelism; this tool performs no threaded work.
    interface.DisableNumThreadsOption();
    interface.AddOptions({Mode, Tile, Seed});
    interface.AddPositionalArguments({Input, ChunkNum, Total});
    return interface;
}

int Runner(const CLI_v2::Results& results)
{
    BamRawReaderConfig readerConfig;

    const std::string mode = results[Mode];
    if (mode == "contiguous") {
        readerConfig.ChunkingMode = ChunkMode::CONTIGUOUS;
    } else if (mode == "scatter") {
        readerConfig.ChunkingMode = ChunkMode::SCATTER;
    } else {
        throw std::runtime_error{"invalid --mode: " + mode};
    }

    const std::int32_t tile{results[Tile]};
    readerConfig.ChunkTileZmws = tile;

    const std::string seedStr = results[Seed];
    readerConfig.ChunkSeed = Tools::ParseIntegerOrThrow<std::uint64_t>(seedStr, "seed");

    // IsUserProvided mirrors the original's tileOrSeedSet flag (set only when the user
    // explicitly passes --tile or --seed, not when the default is applied).
    const bool tileOrSeedSet = results[Tile].IsUserProvided() || results[Seed].IsUserProvided();

    // CLIv2 does not enforce the required positional-argument count, so guard it here.
    const std::vector<std::string>& positional{results.PositionalArguments()};
    if (positional.size() != 3) {
        throw std::runtime_error{"chunk requires exactly three arguments: <input> <chunk> <total>"};
    }

    const std::int32_t chunkNum{Tools::ParseIntegerOrThrow<std::int32_t>(positional[1], "CHUNK")};
    const std::int32_t totalChunks{
        Tools::ParseIntegerOrThrow<std::int32_t>(positional[2], "TOTAL")};

    if (tileOrSeedSet && (readerConfig.ChunkingMode != ChunkMode::SCATTER)) {
        throw std::runtime_error{"--tile/--seed require --mode scatter"};
    }

    readerConfig.ChunkNum = chunkNum;
    readerConfig.TotalChunks = totalChunks;

    BamRawReader reader{std::filesystem::path{positional[0]}, readerConfig};
    const auto& header{reader.Header()};
    std::print("{}", header.ToText());

    for (const auto& view : reader.Records()) {
        WriteViewAsSam(header, view);
    }

    return EXIT_SUCCESS;
}

}  // namespace ChunkTool
}  // namespace Samoa
}  // namespace PacBio
