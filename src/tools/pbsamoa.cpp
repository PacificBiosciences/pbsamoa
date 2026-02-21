#include "bai-build/BaiBuild.hpp"
#include "bai-query/BaiQuery.hpp"
#include "bench/Bench.hpp"
#include "chunk/Chunk.hpp"
#include "dump/Dump.hpp"
#include "zmi-build/ZmiBuild.hpp"
#include "zmi-query/ZmiQuery.hpp"

#include <pbsamoa/PbSamoaLibraryInfo.hpp>

#include <string_view>

#include <cstdio>

namespace {

void PrintUsage()
{
    const std::string version{PacBio::Samoa::LibraryFormattedVersion()};
    std::fprintf(stderr,
                 "pbsamoa — SAM/BAM/BAI toolkit (%s)\n"
                 "\n"
                 "Usage: pbsamoa <command> [args...]\n"
                 "\n"
                 "Commands:\n"
                 "  dump       Convert BAM to SAM text on stdout\n"
                 "  chunk      Dump a chunk of BAM records as SAM text\n"
                 "  bai-build  Build BAI index for a BAM file\n"
                 "  bai-query  Query BAM records by genomic region\n"
                 "  zmi-build  Copy BAM and build ZMI index alongside\n"
                 "  zmi-query  Query BAM records by ZMW hole number\n"
                 "  bench      Benchmark pbsamoa read/write performance\n"
                 "\n"
                 "Examples:\n"
                 "  pbsamoa dump       input.bam              Convert BAM to SAM text\n"
                 "  pbsamoa bai-build  input.bam              Build BAI index\n"
                 "  pbsamoa bai-query  input.bam chr1:1-1000  Region query via BAI\n"
                 "  pbsamoa chunk      input.bam 1 4          Dump chunk 1 of 4 as SAM\n"
                 "  pbsamoa zmi-build  input.bam output.bam   Copy BAM and build ZMI index\n"
                 "  pbsamoa zmi-query  input.bam 42           Query by ZMW hole number\n"
                 "  pbsamoa bench      input.bam              Run benchmarks\n",
                 version.c_str());
}

}  // namespace

int main(int argc, char* argv[])
{
    using namespace PacBio::Samoa;

    if (argc < 2) {
        PrintUsage();
        return 1;
    }

    const std::string_view cmd{argv[1]};

    if (cmd == "dump") {
        return Dump::Runner(argc - 2, argv + 2);
    }
    if (cmd == "chunk") {
        return ChunkTool::Runner(argc - 2, argv + 2);
    }
    if (cmd == "bai-build") {
        return BaiBuild::Runner(argc - 2, argv + 2);
    }
    if (cmd == "bai-query") {
        return BaiQuery::Runner(argc - 2, argv + 2);
    }
    if (cmd == "zmi-build") {
        return ZmiBuild::Runner(argc - 2, argv + 2);
    }
    if (cmd == "zmi-query") {
        return ZmiQuery::Runner(argc - 2, argv + 2);
    }
    if (cmd == "bench") {
        return Bench::Runner(argc - 2, argv + 2);
    }
    if ((cmd == "--help") || (cmd == "-h")) {
        PrintUsage();
        return 0;
    }

    std::fprintf(stderr, "Unknown command: %s\n", argv[1]);
    PrintUsage();
    return 1;
}
