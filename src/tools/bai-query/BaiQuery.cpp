#include "BaiQuery.hpp"

#include "../SamOutput.hpp"

#include <pbsamoa/index/BaiIndex.hpp>
#include <pbsamoa/io/BamRawReader.hpp>

#include <charconv>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

#include <cstdint>
#include <cstdio>
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

std::optional<Region> ParseRegion(std::string_view text)
{
    const std::size_t colonPos{text.find(':')};
    if (colonPos == std::string_view::npos) {
        return std::nullopt;
    }

    const std::string refName{text.substr(0, colonPos)};
    const std::string_view rest{text.substr(colonPos + 1)};

    const std::size_t dashPos{rest.find('-')};
    if (dashPos == std::string_view::npos) {
        return std::nullopt;
    }

    std::int32_t start{0};
    std::int32_t end{0};

    auto result1 = std::from_chars(std::data(rest), std::data(rest) + dashPos, start);
    if (result1.ec != std::errc{}) {
        return std::nullopt;
    }

    auto result2 =
        std::from_chars(std::data(rest) + dashPos + 1, std::data(rest) + std::size(rest), end);
    if (result2.ec != std::errc{}) {
        return std::nullopt;
    }

    // Validate 1-based coordinates
    if ((start < 1) || (end < start)) {
        return std::nullopt;
    }

    return Region{refName, start - 1, end};
}

}  // namespace

int Runner(int argc, char* argv[])
{
    if (argc < 2) {
        std::fprintf(stderr, "Usage: pbsamoa bai-query INPUT REGION\n");
        return EXIT_FAILURE;
    }

    const std::filesystem::path bamPath{argv[0]};
    const std::string regionStr{argv[1]};
    const std::optional<Region> region{ParseRegion(regionStr)};
    if (!region.has_value()) {
        throw std::runtime_error{"invalid region format '" + regionStr +
                                 "' (expected ref:start-end)"};
    }

    const std::filesystem::path baiPath{bamPath.string() + ".bai"};
    if (!std::filesystem::exists(baiPath)) {
        throw std::runtime_error{"index file not found: " + baiPath.string()};
    }
    const BaiIndex index{BaiIndex::FromFile(baiPath)};

    BamRawReader reader{bamPath};
    const SamHeader& header{reader.Header()};

    std::fputs(header.ToText().c_str(), stdout);

    const std::int32_t refId{header.ReferenceId(region->refName)};
    if (refId < 0) {
        throw std::runtime_error{"reference '" + region->refName + "' not found in BAM header"};
    }

    for (const auto& view : reader.Query(index, refId, region->beg, region->end)) {
        WriteViewAsSam(header, view);
    }

    return EXIT_SUCCESS;
}

}  // namespace BaiQuery
}  // namespace Samoa
}  // namespace PacBio
