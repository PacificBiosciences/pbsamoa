#include "Chunk.hpp"
#include "../CliUtils.hpp"
#include "../ParseUtils.hpp"

#include "../SamOutput.hpp"

#include <pbsamoa/io/BamRawReader.hpp>

#include <filesystem>
#include <optional>
#include <print>
#include <string_view>

#include <cstdint>
#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace ChunkTool {
namespace {

void PrintUsage()
{
    std::println(stderr,
                 "Usage: pbsamoa chunk IN.bam CHUNK TOTAL [--mode contiguous|scatter] "
                 "[--tile M] [--seed S]\n"
                 "  --mode   contiguous (default): one contiguous ZMW range per chunk\n"
                 "           scatter: seeded balanced shuffle of M-ZMW tiles, sampled\n"
                 "                    across the whole file (chaotic but deterministic)\n"
                 "  --tile M scatter only: read up to M consecutive ZMWs per seek "
                 "(default 1)\n"
                 "  --seed S scatter only: shuffle seed (default 0)\n"
                 "Running all TOTAL chunks visits every record exactly once.");
}

}  // namespace

int Runner(int argc, char** argv)
{
    BamRawReaderConfig readerConfig;
    const char* bamPath{nullptr};
    std::optional<std::int32_t> chunkNum;
    std::optional<std::int32_t> totalChunks;
    bool tileOrSeedSet{false};

    for (int i{0}; i < argc; ++i) {
        const std::string_view arg{argv[i]};

        if ((arg == "--help") || (arg == "-h")) {
            PrintUsage();
            return EXIT_SUCCESS;
        }

        if (arg == "--mode") {
            const std::string_view mode{Tools::RequireOptionValue(argc, argv, i, "--mode")};
            if (mode == "contiguous") {
                readerConfig.ChunkingMode = ChunkMode::CONTIGUOUS;
            } else if (mode == "scatter") {
                readerConfig.ChunkingMode = ChunkMode::SCATTER;
            } else {
                std::println(stderr, "Error: invalid --mode '{}' (expected contiguous|scatter)",
                             mode);
                return EXIT_FAILURE;
            }
            continue;
        }

        if (arg == "--tile") {
            readerConfig.ChunkTileZmws =
                Tools::ParseIntOption<std::int32_t>(argc, argv, i, "--tile", "tile");
            tileOrSeedSet = true;
            continue;
        }

        if (arg == "--seed") {
            readerConfig.ChunkSeed =
                Tools::ParseIntOption<std::uint64_t>(argc, argv, i, "--seed", "seed");
            tileOrSeedSet = true;
            continue;
        }

        if (!std::empty(arg) && (arg[0] == '-')) {
            std::println(stderr, "Error: unknown option '{}'", arg);
            PrintUsage();
            return EXIT_FAILURE;
        }

        if (!bamPath) {
            bamPath = argv[i];
        } else if (!chunkNum) {
            chunkNum = Tools::ParseIntegerOrThrow<std::int32_t>(arg, "CHUNK");
        } else if (!totalChunks) {
            totalChunks = Tools::ParseIntegerOrThrow<std::int32_t>(arg, "TOTAL");
        } else {
            std::println(stderr, "Error: too many positional arguments");
            PrintUsage();
            return EXIT_FAILURE;
        }
    }

    if (!bamPath || !chunkNum || !totalChunks) {
        PrintUsage();
        return EXIT_FAILURE;
    }

    if (tileOrSeedSet && (readerConfig.ChunkingMode != ChunkMode::SCATTER)) {
        std::println(stderr, "Error: --tile/--seed require --mode scatter");
        return EXIT_FAILURE;
    }

    readerConfig.ChunkNum = *chunkNum;
    readerConfig.TotalChunks = *totalChunks;

    BamRawReader reader{std::filesystem::path{bamPath}, readerConfig};
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
