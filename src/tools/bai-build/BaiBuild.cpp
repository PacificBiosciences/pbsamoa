#include "BaiBuild.hpp"

#include "../../PathUtils.hpp"
#include "../CliUtils.hpp"

#include <pbsamoa/index/BaiIndex.hpp>

#include <algorithm>
#include <filesystem>
#include <format>
#include <print>
#include <stdexcept>
#include <string_view>
#include <thread>
#include <vector>

#include <cstdint>
#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace BaiBuild {

namespace {

void PrintUsage()
{
    std::println(stderr,
                 "Usage: pbsamoa bai-build IN.bam [options]\n"
                 "\n"
                 "Options:\n"
                 "  --threads N       BGZF inflate workers, 0/auto = min(hw,8)   (default auto)");
}

/// Resolve the worker count: 0 means auto = min(hardware_concurrency, 8).
std::size_t ResolveThreads(std::size_t requested)
{
    if (requested != 0) {
        return requested;
    }
    const std::uint32_t hardware{std::thread::hardware_concurrency()};
    const std::size_t available{(hardware == 0) ? std::size_t{1} : hardware};
    return std::min<std::size_t>(available, 8);
}

}  // namespace

int Runner(int argc, char** argv)
{
    std::size_t requestedThreads{0};
    std::vector<std::string_view> positional{};

    for (int i{0}; i < argc; ++i) {
        const std::string_view arg{argv[i]};
        if (arg == "--threads") {
            requestedThreads = Tools::ParseThreadsOption(argc, argv, i);
        } else if (arg.starts_with("--")) {
            throw std::runtime_error{std::format("unknown option: {}", arg)};
        } else {
            positional.push_back(arg);
        }
    }

    if (std::size(positional) != 1) {
        PrintUsage();
        return EXIT_FAILURE;
    }

    const std::filesystem::path bamPath{positional[0]};
    const BaiIndex index{BaiIndex::Build(bamPath, ResolveThreads(requestedThreads))};
    const std::filesystem::path outPath{SidecarPath(bamPath, ".bai")};
    index.ToFile(outPath);
    std::println(stderr, "Index written to {}", outPath.string());
    return EXIT_SUCCESS;
}

}  // namespace BaiBuild
}  // namespace Samoa
}  // namespace PacBio
