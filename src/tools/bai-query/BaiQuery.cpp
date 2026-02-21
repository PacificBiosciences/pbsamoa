#include "BaiQuery.hpp"
#include "../ParseUtils.hpp"

#include "../SamOutput.hpp"

#include <pbsamoa/index/BaiIndex.hpp>
#include <pbsamoa/io/BamRawReader.hpp>

#include <expected>
#include <filesystem>
#include <format>
#include <print>
#include <stdexcept>
#include <string>
#include <string_view>

#include <cstdint>
#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace BaiQuery {
namespace {

struct Region
{
    std::string refName;
    std::int32_t beg;
    std::int32_t end;
};

std::expected<Region, std::string> ParseRegion(std::string_view text)
{
    constexpr std::string_view FORMAT_ERROR = "invalid region format '{}' (expected ref:start-end)";

    const std::size_t colonPos{text.find(':')};
    if (colonPos == std::string_view::npos) {
        return std::unexpected{std::format(FORMAT_ERROR, text)};
    }

    const std::string refName{text.substr(0, colonPos)};
    const std::string_view rest{text.substr(colonPos + 1)};

    const std::size_t dashPos{rest.find('-')};
    if (dashPos == std::string_view::npos) {
        return std::unexpected{std::format(FORMAT_ERROR, text)};
    }

    const auto startResult{
        Tools::ParseInteger<std::int32_t>(rest.substr(0, dashPos), "region start")};
    if (!startResult) {
        return std::unexpected{startResult.error()};
    }

    const auto endResult{Tools::ParseInteger<std::int32_t>(rest.substr(dashPos + 1), "region end")};
    if (!endResult) {
        return std::unexpected{endResult.error()};
    }

    const std::int32_t start{*startResult};
    const std::int32_t end{*endResult};
    if ((start < 1) || (end < start)) {
        return std::unexpected{"region must satisfy start >= 1 and end >= start"};
    }

    return Region{refName, start - 1, end};
}

}  // namespace

int Runner(int argc, char** argv)
{
    if (argc < 2) {
        std::println(stderr, "Usage: pbsamoa bai-query INPUT REGION");
        return EXIT_FAILURE;
    }

    const std::filesystem::path bamPath{argv[0]};
    const std::string regionStr{argv[1]};
    const auto region{ParseRegion(regionStr)};
    if (!region) {
        throw std::runtime_error{region.error()};
    }

    const std::filesystem::path baiPath{bamPath.string() + ".bai"};
    if (!std::filesystem::exists(baiPath)) {
        throw std::runtime_error{std::format("index file not found: {}", baiPath.string())};
    }
    const BaiIndex index{BaiIndex::FromFile(baiPath)};

    BamRawReader reader{bamPath};
    const SamHeader& header{reader.Header()};

    std::print("{}", header.ToText());

    const std::int32_t refId{header.ReferenceId(region->refName)};
    if (refId < 0) {
        throw std::runtime_error{
            std::format("reference '{}' not found in BAM header", region->refName)};
    }

    for (const auto& view : reader.Query(index, refId, region->beg, region->end)) {
        WriteViewAsSam(header, view);
    }

    return EXIT_SUCCESS;
}

}  // namespace BaiQuery
}  // namespace Samoa
}  // namespace PacBio
