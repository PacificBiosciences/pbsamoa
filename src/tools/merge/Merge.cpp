#include "Merge.hpp"

#include "../CliUtils.hpp"

#include <pbsamoa/io/BamMerge.hpp>

#include <filesystem>
#include <format>
#include <print>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <cstdint>
#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace MergeTool {

namespace {

void PrintUsage()
{
    std::println(stderr,
                 "Usage: pbsamoa merge OUT.bam IN1.bam [IN2.bam ...] [options]\n"
                 "\n"
                 "Merge already-sorted BAM files (identical @SQ required) into one.\n"
                 "\n"
                 "Options:\n"
                 "  --order ORDER     coordinate | queryname | tag   (default coordinate)\n"
                 "  --tag XX          2-char tag inputs are sorted by (required iff --order tag)\n"
                 "  --threads N       worker threads, 0/auto = min(hw,8)\n"
                 "  --compression L   output BGZF level [1,12]         (default 6)");
}

std::string BuildCommandLine(int argc, char** argv)
{
    std::string commandLine{"pbsamoa merge"};
    for (int i{0}; i < argc; ++i) {
        commandLine += ' ';
        commandLine += argv[i];
    }
    return commandLine;
}

}  // namespace

int Runner(int argc, char** argv)
{
    MergeConfig config{};
    std::vector<std::string_view> positional{};
    bool tagProvided{false};

    for (int i{0}; i < argc; ++i) {
        const std::string_view arg{argv[i]};
        if (arg == "--order") {
            const std::string_view value{Tools::RequireOptionValue(argc, argv, i, "--order")};
            if (value == "coordinate") {
                config.Order = SortOrder::COORDINATE;
            } else if (value == "queryname") {
                config.Order = SortOrder::QUERY_NAME;
            } else if (value == "tag") {
                config.Order = SortOrder::TAG;
            } else {
                throw std::runtime_error{std::format("invalid --order: {}", value)};
            }
        } else if (arg == "--tag") {
            const std::string_view value{Tools::RequireOptionValue(argc, argv, i, "--tag")};
            if (std::size(value) != 2) {
                throw std::runtime_error{
                    std::format("--tag must be exactly 2 characters: {}", value)};
            }
            config.Tag = {value[0], value[1]};
            tagProvided = true;
        } else if (arg == "--threads") {
            const std::int32_t threads{Tools::ParseIntegerOrThrow<std::int32_t>(
                Tools::RequireOptionValue(argc, argv, i, "--threads"), "--threads")};
            if (threads < 0) {
                throw std::runtime_error{"--threads must be >= 0"};
            }
            config.NumThreads = static_cast<std::size_t>(threads);
        } else if (arg == "--compression") {
            const int level{Tools::ParseIntegerOrThrow<int>(
                Tools::RequireOptionValue(argc, argv, i, "--compression"), "--compression")};
            if ((level < 1) || (level > 12)) {
                throw std::runtime_error{"--compression must be in [1, 12]"};
            }
            config.CompressionLevel = level;
        } else if (arg.starts_with("--")) {
            throw std::runtime_error{std::format("unknown option: {}", arg)};
        } else {
            positional.push_back(arg);
        }
    }

    // First positional is the output; the rest are inputs.
    if (std::size(positional) < 2) {
        PrintUsage();
        return EXIT_FAILURE;
    }
    if ((config.Order == SortOrder::TAG) && !tagProvided) {
        throw std::runtime_error{"--order tag requires --tag XX"};
    }

    const std::filesystem::path output{positional.front()};
    std::vector<std::filesystem::path> inputs{};
    inputs.reserve(std::size(positional) - 1);
    for (std::size_t i{1}; i < std::size(positional); ++i) {
        inputs.emplace_back(positional[i]);
    }

    config.CommandLine = BuildCommandLine(argc, argv);

    const MergeStats stats{MergeBam(inputs, output, config)};
    std::println(stderr, "pbsamoa merge: {} records from {} input(s)", stats.NumRecords,
                 stats.NumInputs);
    return EXIT_SUCCESS;
}

}  // namespace MergeTool
}  // namespace Samoa
}  // namespace PacBio
