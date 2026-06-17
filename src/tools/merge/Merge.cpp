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

}  // namespace

int Runner(int argc, char** argv)
{
    MergeConfig config{};
    std::vector<std::string_view> positional{};
    bool tagProvided{false};

    for (int i{0}; i < argc; ++i) {
        const std::string_view arg{argv[i]};
        if (arg == "--order") {
            config.Order = Tools::ParseSortOrder(argc, argv, i);
        } else if (arg == "--tag") {
            config.Tag = Tools::ParseSortTag(argc, argv, i);
            tagProvided = true;
        } else if (arg == "--threads") {
            config.NumThreads = Tools::ParseThreadsOption(argc, argv, i);
        } else if (arg == "--compression") {
            config.CompressionLevel = Tools::ParseCompressionLevelOption(argc, argv, i);
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

    config.CommandLine = Tools::BuildCommandLine("merge", argc, argv);

    const MergeStats stats{MergeBam(inputs, output, config)};
    std::println(stderr, "pbsamoa merge: {} records from {} input(s)", stats.NumRecords,
                 stats.NumInputs);
    return EXIT_SUCCESS;
}

}  // namespace MergeTool
}  // namespace Samoa
}  // namespace PacBio
