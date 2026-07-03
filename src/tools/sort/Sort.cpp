#include "Sort.hpp"

#include "../CliUtils.hpp"

#include <pbsamoa/PbSamoaLibraryInfo.hpp>
#include <pbsamoa/io/BamSort.hpp>

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
namespace SortTool {
namespace {

const CLI_v2::Option Order{
    R"({
    "names" : ["order"],
    "description" : "Sort order.",
    "type" : "string",
    "choices" : ["coordinate", "queryname", "tag"],
    "default" : "coordinate"
})"};

const CLI_v2::Option Tag{
    R"({
    "names" : ["tag"],
    "description" : "Two-character tag to sort by (required with --order tag).",
    "type" : "string"
})"};

const CLI_v2::Option Memory{
    R"({
    "names" : ["memory"],
    "description" : "Total RAM budget, with an optional K/M/G suffix.",
    "type" : "string",
    "default" : "768M"
})"};

const CLI_v2::Option Threads{
    R"({
    "names" : ["threads"],
    "description" : "Worker threads; 0 = auto (min of hardware concurrency and 8).",
    "type" : "unsigned integer",
    "default" : 0
})"};

const CLI_v2::Option Compression{
    R"({
    "names" : ["compression"],
    "description" : "Output BGZF compression level in [1, 12].",
    "type" : "integer",
    "default" : 6
})"};

const CLI_v2::Option TempDir{
    R"({
    "names" : ["temp-dir"],
    "description" : "Directory for temporary run files (default: output directory).",
    "type" : "dir"
})"};

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
    CLI_v2::Interface interface{"pbsamoa sort", "Sort a BAM file (coordinate/queryname/tag).",
                                LibraryFormattedVersion()};
    // `sort` keeps its own --threads (0 = auto, capped at min(hw, 8)); the built-in
    // --num-threads resolves 0 to the raw hardware count, which would drop that cap.
    interface.DisableNumThreadsOption();
    interface.AddOptions({Order, Tag, Memory, Threads, Compression, TempDir});
    interface.AddPositionalArguments({Input, Output});
    return interface;
}

int Runner(const CLI_v2::Results& results)
{
    SortConfig config{};

    // Copy-init (not brace-init): a braced std::string from Result is ambiguous
    // because Result also converts to a char-like type for string's char-list ctor.
    const std::string order = results[Order];
    config.Order = Tools::ParseSortOrder(order);

    const std::string memory = results[Memory];
    config.MaxMemory = Tools::ParseMemory(memory);

    const std::uint32_t threads{results[Threads]};
    config.NumThreads = threads;

    const std::int32_t compression{results[Compression]};
    config.CompressionLevel = Tools::CheckCompressionLevel(compression);

    const std::string tempDir = results[TempDir];
    if (!tempDir.empty()) {
        config.TempDir = std::filesystem::path{tempDir};
    }

    // CLIv2 does not enforce the required positional-argument count, so guard it
    // here before indexing (matches the tool's prior "exactly two operands" rule).
    const std::vector<std::string>& positional{results.PositionalArguments()};
    if (positional.size() != 2) {
        throw std::runtime_error{"sort requires exactly two arguments: <input> <output>"};
    }

    const std::string tag = results[Tag];
    if (!tag.empty()) {
        config.Tag = Tools::ParseSortTag(tag);
    }
    if ((config.Order == SortOrder::TAG) && tag.empty()) {
        throw std::runtime_error{"--order tag requires --tag XX"};
    }

    config.CommandLine = std::string{"pbsamoa "} + results.InputCommandLine();

    const SortStats stats{SortBam(std::filesystem::path{positional[0]},
                                  std::filesystem::path{positional[1]}, config)};
    std::println(stderr, "pbsamoa sort: {} records, {} spilled run(s)", stats.NumRecords,
                 stats.NumRuns);
    return EXIT_SUCCESS;
}

}  // namespace SortTool
}  // namespace Samoa
}  // namespace PacBio
