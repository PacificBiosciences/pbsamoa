#include "bai-build/BaiBuild.hpp"
#include "bai-query/BaiQuery.hpp"
#include "bench/Bench.hpp"
#include "chunk/Chunk.hpp"
#include "convert/Convert.hpp"
#include "dump/Dump.hpp"
#include "zmi-build/ZmiBuild.hpp"
#include "zmi-query/ZmiQuery.hpp"

#include <pbsamoa/PbSamoaLibraryInfo.hpp>

#include <array>
#include <exception>
#include <print>
#include <ranges>
#include <string_view>

namespace {

using Runner = int (*)(int argc, char** argv);

struct CommandSpec
{
    std::string_view Name;
    Runner Run;
};

constexpr std::array COMMANDS{
    CommandSpec{"dump", PacBio::Samoa::Dump::Runner},
    CommandSpec{"chunk", PacBio::Samoa::ChunkTool::Runner},
    CommandSpec{"convert", PacBio::Samoa::Convert::Runner},
    CommandSpec{"bai-build", PacBio::Samoa::BaiBuild::Runner},
    CommandSpec{"bai-query", PacBio::Samoa::BaiQuery::Runner},
    CommandSpec{"zmi-build", PacBio::Samoa::ZmiBuild::Runner},
    CommandSpec{"zmi-query", PacBio::Samoa::ZmiQuery::Runner},
    CommandSpec{"bench", PacBio::Samoa::Bench::Runner},
};

void PrintUsage()
{
    std::print(stderr,
               "pbsamoa - SAM/BAM/BAI toolkit ({})\n"
               "\n"
               "Usage: pbsamoa <command> [args...]\n"
               "\n"
               "Commands:\n"
               "  dump       Convert BAM to SAM text on stdout\n"
               "  chunk      Dump a chunk of BAM records as SAM text\n"
               "  convert    Convert BAM/SAM to CRAM format\n"
               "  bai-build  Build BAI index for a BAM file\n"
               "  bai-query  Query BAM records by genomic region\n"
               "  zmi-build  Copy BAM and build ZMI index alongside\n"
               "  zmi-query  Query BAM records by ZMW hole number\n"
               "  bench      Benchmark pbsamoa read/write performance\n"
               "\n"
               "Examples:\n"
               "  pbsamoa dump       input.bam              Convert BAM to SAM text\n"
               "  pbsamoa convert    input.bam output.cram  Convert BAM/SAM to CRAM\n"
               "  pbsamoa bai-build  input.bam              Build BAI index\n"
               "  pbsamoa bai-query  input.bam chr1:1-1000  Region query via BAI\n"
               "  pbsamoa chunk      input.bam 1 4          Dump chunk 1 of 4 as SAM\n"
               "  pbsamoa zmi-build  input.bam output.bam   Copy BAM and build ZMI index\n"
               "  pbsamoa zmi-query  input.bam 42           Query by ZMW hole number\n"
               "  pbsamoa bench      input.bam              Run benchmarks\n",
               PacBio::Samoa::LibraryFormattedVersion());
}

}  // namespace

int main(int argc, char* argv[])
{
    using namespace PacBio::Samoa;

    try {
        if (argc < 2) {
            PrintUsage();
            return 1;
        }

        const std::string_view cmd{argv[1]};

        if (const auto command{std::ranges::find(COMMANDS, cmd, &CommandSpec::Name)};
            command != std::ranges::end(COMMANDS)) {
            return command->Run(argc - 2, argv + 2);
        }

        if ((cmd == "--help") || (cmd == "-h")) {
            PrintUsage();
            return 0;
        }

        std::println(stderr, "Unknown command: {}", argv[1]);
        PrintUsage();
        return 1;
    } catch (const std::exception& e) {
        std::println(stderr, "Error: {}", e.what());
        return 1;
    }
}
