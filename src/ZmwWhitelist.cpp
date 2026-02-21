#include <pbsamoa/index/ZmwWhitelist.hpp>

#include <algorithm>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {

ZmwWhitelist::ZmwWhitelist(std::vector<std::int32_t> zmws) : data_{std::move(zmws)} {}

ZmwWhitelist::ZmwWhitelist(std::vector<ZmwIdentity> ids) : data_{std::move(ids)} {}

std::vector<std::int64_t> ZmwWhitelist::Resolve(const ZmwIndex& index) const
{
    std::vector<std::int64_t> result;

    if (const auto* zmws = std::get_if<std::vector<std::int32_t>>(&data_); zmws != nullptr) {
        for (const std::int32_t zmw : *zmws) {
            const std::vector<std::int64_t> offsets{index.Find(zmw)};
            result.insert(std::end(result), std::begin(offsets), std::end(offsets));
        }
    } else if (const auto* ids = std::get_if<std::vector<ZmwIdentity>>(&data_); ids != nullptr) {
        result = index.Find(*ids);
    }

    // Sort and deduplicate (file order = ascending virtual offset)
    std::ranges::sort(result);
    const auto [first, last] = std::ranges::unique(result);
    result.erase(first, last);

    return result;
}

}  // namespace Samoa
}  // namespace PacBio
