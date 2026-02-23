#include "ZmiQuery.hpp"
#include "../ParseUtils.hpp"

#include "../SamOutput.hpp"

#include <pbsamoa/core/Bgzf.hpp>
#include <pbsamoa/index/ZmwIndex.hpp>
#include <pbsamoa/io/BamRawReader.hpp>

#include <filesystem>
#include <print>
#include <stdexcept>
#include <vector>

#include <cstdint>
#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace ZmiQuery {

int Runner(int argc, char* argv[])
{
    if (argc < 2) {
        std::println(stderr, "Usage: pbsamoa zmi-query IN.bam ZMW");
        return EXIT_FAILURE;
    }

    const std::filesystem::path bamPath{argv[0]};
    const auto zmwResult{Tools::ParseInteger<std::int32_t>(argv[1], "ZMW")};
    if (!zmwResult) {
        throw std::runtime_error{zmwResult.error()};
    }

    const ZmwIndex index{ZmwIndex::Open(bamPath)};
    const std::vector<std::int64_t> offsets{index.Find(*zmwResult)};

    BamRawReader reader{bamPath};
    const auto& header{reader.Header()};

    std::print("{}", header.ToText());

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
