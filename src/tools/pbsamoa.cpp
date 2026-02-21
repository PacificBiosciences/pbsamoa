#include "bai-build/BaiBuild.hpp"
#include "bai-query/BaiQuery.hpp"
#include "bench/Bench.hpp"
#include "chunk/Chunk.hpp"
#include "dump/Dump.hpp"
#include "zmi-build/ZmiBuild.hpp"
#include "zmi-query/ZmiQuery.hpp"

#include <pbcopper/cli2/CLI.h>

int main(int argc, char* argv[])
{
    using namespace PacBio::Samoa;

    PacBio::CLI_v2::MultiToolInterface app{"pbsamoa", "SAM/BAM/BAI toolkit", "0.1.0"};

    app.AddTools({
        {"dump", Dump::CreateInterface(), &Dump::Runner},
        {"chunk", ChunkTool::CreateInterface(), &ChunkTool::Runner},
        {"bai-build", BaiBuild::CreateInterface(), &BaiBuild::Runner},
        {"bai-query", BaiQuery::CreateInterface(), &BaiQuery::Runner},
        {"zmi-build", ZmiBuild::CreateInterface(), &ZmiBuild::Runner},
        {"zmi-query", ZmiQuery::CreateInterface(), &ZmiQuery::Runner},
        {"bench", Bench::CreateInterface(), &Bench::Runner},
    });

    app.HelpFooter(
        "Examples:\n"
        "  pbsamoa dump       input.bam              Convert BAM to SAM text\n"
        "  pbsamoa bai-build  input.bam              Build BAI index\n"
        "  pbsamoa bai-query  input.bam chr1:1-1000  Region query via BAI\n"
        "  pbsamoa chunk      input.bam 1 4          Dump chunk 1 of 4 as SAM\n"
        "  pbsamoa zmi-build  input.bam output.bam   Copy BAM and build ZMI index\n"
        "  pbsamoa zmi-query  input.bam 42           Query by ZMW hole number\n"
        "  pbsamoa bench      input.bam              Run benchmarks\n");

    return PacBio::CLI_v2::Run(argc, argv, app);
}
