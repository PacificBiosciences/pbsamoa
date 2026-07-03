#include "BaiQuery.hpp"

#include "../ParseUtils.hpp"

#include "../../PathUtils.hpp"
#include "../SamOutput.hpp"

#include <pbsamoa/PbSamoaLibraryInfo.hpp>
#include <pbsamoa/index/BaiIndex.hpp>
#include <pbsamoa/io/BamRawReader.hpp>

#include <pbcopper/cli2/Interface.h>
#include <pbcopper/cli2/PositionalArgument.h>
#include <pbcopper/cli2/Results.h>

#include <filesystem>
#include <format>
#include <print>
#include <stdexcept>
#include <string>
#include <vector>

#include <cstdint>
#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace BaiQuery {
namespace {

const CLI_v2::PositionalArgument Input{
    R"({
    "name" : "input",
    "description" : "Input BAM file.",
    "type" : "file"
})"};

const CLI_v2::PositionalArgument Region{
    R"({
    "name" : "region",
    "description" : "Genomic region to query (e.g. ref:start-end).",
    "type" : "string"
})"};

}  // namespace

CLI_v2::Interface CreateInterface()
{
    CLI_v2::Interface interface{"pbsamoa bai-query", "Query BAM records by genomic region.",
                                LibraryFormattedVersion()};
    // bai-query exposes no threading options; disable the built-in --num-threads so
    // it does not appear in --help (the original tool had no thread flag).
    interface.DisableNumThreadsOption();
    interface.AddPositionalArguments({Input, Region});
    return interface;
}

int Runner(const CLI_v2::Results& results)
{
    // CLIv2 does not enforce the required positional-argument count, so guard it
    // here before indexing (matches the tool's prior "exactly two operands" rule).
    const std::vector<std::string>& positional{results.PositionalArguments()};
    if (positional.size() != 2) {
        throw std::runtime_error{"bai-query requires exactly two arguments: <input> <region>"};
    }

    const std::filesystem::path bamPath{positional[0]};
    // Region is a free-form string (e.g. "ref:1-1000"); keep handler-side parsing via
    // ParseRegion rather than a CLIv2 type (CLIv2 has no region type).
    const auto region{Tools::ParseRegion(positional[1], "ref:start-end")};
    if (!region) {
        throw std::runtime_error{region.error()};
    }

    const std::filesystem::path baiPath{SidecarPath(bamPath, ".bai")};
    if (!std::filesystem::exists(baiPath)) {
        throw std::runtime_error{std::format("index file not found: {}", baiPath.string())};
    }
    const BaiIndex index{BaiIndex::FromFile(baiPath)};

    BamRawReader reader{bamPath};
    const SamHeader& header{reader.Header()};

    std::print("{}", header.ToText());

    const std::int32_t refId{header.ReferenceId(region->RefName)};
    if (refId < 0) {
        throw std::runtime_error{
            std::format("reference '{}' not found in BAM header", region->RefName)};
    }

    for (const auto& view : reader.Query(index, refId, region->Beg, region->End)) {
        WriteViewAsSam(header, view);
    }

    return EXIT_SUCCESS;
}

}  // namespace BaiQuery
}  // namespace Samoa
}  // namespace PacBio
