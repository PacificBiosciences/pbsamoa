#include <pbsamoa/index/ZmwWhitelist.hpp>

#include <algorithm>
#include <iterator>
#include <variant>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {

namespace {

void SortAndDeduplicate(std::vector<std::int64_t>& result)
{
    std::ranges::sort(result);
    const auto [first, last] = std::ranges::unique(result);
    result.erase(first, last);
}

}  // namespace

ZmwWhitelist::ZmwWhitelist(std::vector<std::int32_t> zmws) : data_{std::move(zmws)} {}

ZmwWhitelist::ZmwWhitelist(std::vector<ZmwIdentity> ids) : data_{std::move(ids)} {}

std::vector<std::int64_t> ZmwWhitelist::Resolve(const ZmwIndex& index) const
{
    std::vector<std::int64_t> result;

    if (const std::vector<std::int32_t>* zmws{std::get_if<std::vector<std::int32_t>>(&data_)};
        zmws) {
        for (const std::int32_t zmw : *zmws) {
            std::ranges::copy(index.Find(zmw), std::back_inserter(result));
        }
    } else {
        result = index.Find(std::get<std::vector<ZmwIdentity>>(data_));
    }

    // Sort and deduplicate (file order = ascending virtual offset)
    SortAndDeduplicate(result);
    return result;
}

}  // namespace Samoa
}  // namespace PacBio
