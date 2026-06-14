#include "Sort.hpp"

#include "../CliUtils.hpp"

#include <pbsamoa/io/BamSort.hpp>

#include <array>
#include <filesystem>
#include <format>
#include <optional>
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

std::string BuildCommandLine(int argc, char** argv)
{
    std::string commandLine{"pbsamoa sort"};
    for (int i{0}; i < argc; ++i) {
        commandLine += ' ';
        commandLine += argv[i];
    }
    return commandLine;
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
        } else if (arg == "--memory") {
            config.MaxMemory = ParseMemory(Tools::RequireOptionValue(argc, argv, i, "--memory"));
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

    config.CommandLine = BuildCommandLine(argc, argv);

    const SortStats stats{SortBam(std::filesystem::path{positional[0]},
                                  std::filesystem::path{positional[1]}, config)};
    std::println(stderr, "pbsamoa sort: {} records, {} spilled run(s)", stats.NumRecords,
                 stats.NumRuns);
    return EXIT_SUCCESS;
}

}  // namespace SortTool
}  // namespace Samoa
}  // namespace PacBio
