#include "ZmiQuery.hpp"

#include "../SamOutput.hpp"

#include <pbsamoa/core/Bgzf.hpp>
#include <pbsamoa/index/ZmwIndex.hpp>
#include <pbsamoa/io/BamRawReader.hpp>

#include <charconv>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <cstdint>
#include <cstdio>
#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace ZmiQuery {

int Runner(int argc, char* argv[])
{
    if (argc < 2) {
        std::fprintf(stderr, "Usage: pbsamoa zmi-query IN.bam ZMW\n");
        return EXIT_FAILURE;
    }

    const std::filesystem::path bamPath{argv[0]};

    std::int32_t zmw{0};
    {
        const std::string_view zmwStr{argv[1]};
        auto r{std::from_chars(std::data(zmwStr), std::data(zmwStr) + std::size(zmwStr), zmw)};
        if (r.ec != std::errc{}) {
            throw std::runtime_error{"invalid ZMW: " + std::string{zmwStr}};
        }
    }

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
