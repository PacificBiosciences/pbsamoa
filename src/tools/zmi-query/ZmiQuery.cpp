#include "ZmiQuery.hpp"
#include "../ParseUtils.hpp"

#include "../SamOutput.hpp"

#include <pbsamoa/PbSamoaLibraryInfo.hpp>
#include <pbsamoa/core/Bgzf.hpp>
#include <pbsamoa/index/ZmwIndex.hpp>
#include <pbsamoa/io/BamRawReader.hpp>

#include <pbcopper/cli2/Interface.h>
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
namespace ZmiQuery {
namespace {

const CLI_v2::PositionalArgument Input{
    R"({
    "name" : "input",
    "description" : "Input BAM file.",
    "type" : "file"
})"};

const CLI_v2::PositionalArgument Zmw{
    R"({
    "name" : "zmw",
    "description" : "ZMW hole number to query.",
    "type" : "string"
})"};

}  // namespace

CLI_v2::Interface CreateInterface()
{
    CLI_v2::Interface interface{"pbsamoa zmi-query", "Query BAM records by ZMW hole number.",
                                LibraryFormattedVersion()};
    // zmi-query has no thread option; disable the built-in --num-threads so it
    // does not appear in --help.
    interface.DisableNumThreadsOption();
    interface.AddPositionalArguments({Input, Zmw});
    return interface;
}

int Runner(const CLI_v2::Results& results)
{
    // CLIv2 does not enforce the required positional count; guard before indexing.
    const std::vector<std::string>& positional{results.PositionalArguments()};
    if (positional.size() != 2) {
        throw std::runtime_error{"zmi-query requires exactly two arguments: <input> <zmw>"};
    }

    const std::filesystem::path bamPath{positional[0]};
    const std::int32_t zmw{Tools::ParseIntegerOrThrow<std::int32_t>(positional[1], "ZMW")};

    const ZmwIndex index{ZmwIndex::Open(bamPath)};
    const std::vector<std::int64_t> offsets{index.Find(zmw)};

    BamRawReader reader{bamPath};
    const auto& header{reader.Header()};

    std::print("{}", header.ToText());

    for (const std::int64_t offset : offsets) {
        reader.Seek(VirtualOffset(offset));
        if (const auto view{reader.ReadRecord()}; view) {
            WriteViewAsSam(header, *view);
        }
    }

    return EXIT_SUCCESS;
}

}  // namespace ZmiQuery
}  // namespace Samoa
}  // namespace PacBio
