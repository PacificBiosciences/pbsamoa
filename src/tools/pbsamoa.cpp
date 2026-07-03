#include "bai-build/BaiBuild.hpp"
#include "bai-query/BaiQuery.hpp"
#include "bench/Bench.hpp"
#include "chunk/Chunk.hpp"
#include "convert/Convert.hpp"
#include "dump/Dump.hpp"
#include "merge/Merge.hpp"
#include "sort/Sort.hpp"
#include "zmi-build/ZmiBuild.hpp"
#include "zmi-index/ZmiIndex.hpp"
#include "zmi-query/ZmiQuery.hpp"

#include <pbsamoa/PbSamoaLibraryInfo.hpp>

#include <pbcopper/cli2/CLI.h>
#include <pbcopper/cli2/MultiToolInterface.h>

namespace {

PacBio::CLI_v2::MultiToolInterface CreateMultiToolInterface()
{
    using namespace PacBio::Samoa;

    PacBio::CLI_v2::MultiToolInterface interface{"pbsamoa", "SAM/BAM/BAI toolkit",
                                                 LibraryFormattedVersion()};
    interface.AddTools({
        {"dump", Dump::CreateInterface(), &Dump::Runner},
        {"chunk", ChunkTool::CreateInterface(), &ChunkTool::Runner},
        {"sort", SortTool::CreateInterface(), &SortTool::Runner},
        {"merge", MergeTool::CreateInterface(), &MergeTool::Runner},
        {"convert", Convert::CreateInterface(), &Convert::Runner},
        {"bai-build", BaiBuild::CreateInterface(), &BaiBuild::Runner},
        {"bai-query", BaiQuery::CreateInterface(), &BaiQuery::Runner},
        {"zmi-build", ZmiBuild::CreateInterface(), &ZmiBuild::Runner},
        {"zmi-index", ZmiIndex::CreateInterface(), &ZmiIndex::Runner},
        {"zmi-query", ZmiQuery::CreateInterface(), &ZmiQuery::Runner},
        {"bench", Bench::CreateInterface(), &Bench::Runner},
    });
    interface.HelpFooter(
        "Examples:\n"
        "  pbsamoa dump       input.bam              Convert BAM to SAM text\n"
        "  pbsamoa convert    input.bam output.cram  Convert BAM/SAM to CRAM\n"
        "  pbsamoa bai-build  input.bam              Build BAI index\n"
        "  pbsamoa bai-query  input.bam chr1:1-1000  Region query via BAI\n"
        "  pbsamoa chunk      input.bam 1 4          Dump chunk 1 of 4 as SAM\n"
        "  pbsamoa sort       input.bam sorted.bam   Coordinate-sort a BAM\n"
        "  pbsamoa merge      out.bam a.bam b.bam     Merge or concatenate BAMs\n"
        "  pbsamoa zmi-build  input.bam output.bam   Copy BAM and build ZMI index\n"
        "  pbsamoa zmi-index  input.bam              Build .zmi sidecar (no rewrite)\n"
        "  pbsamoa zmi-query  input.bam 42           Query by ZMW hole number\n"
        "  pbsamoa bench      input.bam              Run benchmarks");
    return interface;
}

}  // namespace

int main(int argc, char* argv[])
{
    return PacBio::CLI_v2::Run(argc, argv, CreateMultiToolInterface());
}
