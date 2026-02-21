#include "ZmiQuery.hpp"

#include "../SamOutput.hpp"

#include <pbsamoa/core/Bgzf.hpp>
#include <pbsamoa/index/ZmwIndex.hpp>
#include <pbsamoa/io/BamRawReader.hpp>

#include <pbcopper/cli2/CLI.h>

#include <filesystem>
#include <string>
#include <vector>

#include <cstdint>
#include <cstdio>
#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace ZmiQuery {
namespace {

// clang-format off
const PacBio::CLI_v2::PositionalArgument ZmiQueryInput{
R"({
    "name" : "IN.bam",
    "description" : "Input BAM file (requires .zmi index alongside).",
    "type" : "file"
})"};

const PacBio::CLI_v2::PositionalArgument ZmiQueryZmw{
R"({
    "name" : "ZMW",
    "description" : "ZMW hole number to query.",
    "type" : "integer"
})"};
// clang-format on

}  // namespace

PacBio::CLI_v2::Interface CreateInterface()
{
    PacBio::CLI_v2::Interface iface{"zmi-query", "Query BAM records by ZMW hole number", "0.1.0"};
    iface.DisableLogFileOption();
    iface.DisableNumThreadsOption();
    iface.AddPositionalArguments({ZmiQueryInput, ZmiQueryZmw});
    return iface;
}

int Runner(const PacBio::CLI_v2::Results& results)
{
    const std::filesystem::path bamPath{results[ZmiQueryInput]};
    const std::int32_t zmw{std::stoi(std::string{results[ZmiQueryZmw]})};

    const ZmwIndex index{ZmwIndex::Open(bamPath)};
    const std::vector<std::int64_t> offsets{index.Find(zmw)};

    BamRawReader reader{bamPath};
    const auto& header{reader.Header()};

    std::fputs(header.ToText().c_str(), stdout);

    for (const std::int64_t offset : offsets) {
        reader.Seek(VirtualOffset(offset));
        const auto view{reader.ReadRecord()};
        if (view.has_value()) {
            WriteViewAsSam(header, *view);
        }
    }

    return EXIT_SUCCESS;
}

}  // namespace ZmiQuery
}  // namespace Samoa
}  // namespace PacBio
