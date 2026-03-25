#include "BaiQuery.hpp"
#include "../ParseUtils.hpp"

#include "../../PathUtils.hpp"
#include "../SamOutput.hpp"

#include <pbsamoa/index/BaiIndex.hpp>
#include <pbsamoa/io/BamRawReader.hpp>

#include <filesystem>
#include <format>
#include <print>
#include <stdexcept>
#include <string_view>

#include <cstdint>
#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace BaiQuery {

int Runner(int argc, char** argv)
{
    if (argc < 2) {
        std::println(stderr, "Usage: pbsamoa bai-query INPUT REGION");
        return EXIT_FAILURE;
    }

    const std::filesystem::path bamPath{argv[0]};
    const std::string_view regionText{argv[1]};
    const auto region{Tools::ParseRegion(regionText, "ref:start-end")};
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
