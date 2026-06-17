#include "Sort.hpp"

#include "../CliUtils.hpp"

#include <pbsamoa/io/BamSort.hpp>

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
namespace SortTool {

namespace {

void PrintUsage()
{
    std::println(stderr,
                 "Usage: pbsamoa sort IN.bam OUT.bam [options]\n"
                 "\n"
                 "Options:\n"
                 "  --order ORDER     coordinate | queryname | tag   (default coordinate)\n"
                 "  --tag XX          2-char tag to sort by (required iff --order tag)\n"
                 "  --memory SIZE     total RAM budget, K/M/G suffix   (default 768M)\n"
                 "  --threads N       worker threads, 0/auto = min(hw,8)\n"
                 "  --compression L   output BGZF level [1,12]         (default 6)\n"
                 "  --temp-dir DIR    directory for temporary run files");
}

/// Parse a memory size with an optional K/M/G suffix (base 1024) into bytes.
ByteLimit ParseMemory(std::string_view text)
{
    if (text.empty()) {
        throw std::runtime_error{"invalid --memory: empty value"};
    }

    std::uint64_t multiplier{1};
    std::string_view digits{text};
    switch (text.back()) {
        case 'k':
        case 'K':
            multiplier = std::uint64_t{1024};
            digits = text.substr(0, std::size(text) - 1);
            break;
        case 'm':
        case 'M':
            multiplier = std::uint64_t{1024} * 1024;
            digits = text.substr(0, std::size(text) - 1);
            break;
        case 'g':
        case 'G':
            multiplier = std::uint64_t{1024} * 1024 * 1024;
            digits = text.substr(0, std::size(text) - 1);
            break;
        default:
            break;
    }

    const std::uint64_t value{Tools::ParseIntegerOrThrow<std::uint64_t>(digits, "--memory")};
    return ByteLimit{value * multiplier};
}

}  // namespace

int Runner(int argc, char** argv)
{
    SortConfig config{};
    std::vector<std::string_view> positional{};
    bool tagProvided{false};

    for (int i{0}; i < argc; ++i) {
        const std::string_view arg{argv[i]};
        if (arg == "--order") {
            config.Order = Tools::ParseSortOrder(argc, argv, i);
        } else if (arg == "--tag") {
            config.Tag = Tools::ParseSortTag(argc, argv, i);
            tagProvided = true;
        } else if (arg == "--memory") {
            config.MaxMemory = ParseMemory(Tools::RequireOptionValue(argc, argv, i, "--memory"));
        } else if (arg == "--threads") {
            config.NumThreads = Tools::ParseThreadsOption(argc, argv, i);
        } else if (arg == "--compression") {
            config.CompressionLevel = Tools::ParseCompressionLevelOption(argc, argv, i);
        } else if (arg == "--temp-dir") {
            config.TempDir =
                std::filesystem::path{Tools::RequireOptionValue(argc, argv, i, "--temp-dir")};
        } else if (arg.starts_with("--")) {
            throw std::runtime_error{std::format("unknown option: {}", arg)};
        } else {
            positional.push_back(arg);
        }
    }

    if (std::size(positional) != 2) {
        PrintUsage();
        return EXIT_FAILURE;
    }
    if ((config.Order == SortOrder::TAG) && !tagProvided) {
        throw std::runtime_error{"--order tag requires --tag XX"};
    }

    config.CommandLine = Tools::BuildCommandLine("sort", argc, argv);

    const SortStats stats{SortBam(std::filesystem::path{positional[0]},
                                  std::filesystem::path{positional[1]}, config)};
    std::println(stderr, "pbsamoa sort: {} records, {} spilled run(s)", stats.NumRecords,
                 stats.NumRuns);
    return EXIT_SUCCESS;
}

}  // namespace SortTool
}  // namespace Samoa
}  // namespace PacBio
