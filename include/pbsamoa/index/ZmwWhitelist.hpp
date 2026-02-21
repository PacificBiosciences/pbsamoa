#ifndef PBSAMOA_INDEX_ZMWWHITELIST_HPP
#define PBSAMOA_INDEX_ZMWWHITELIST_HPP

#include <pbsamoa/index/ZmwIndex.hpp>

#include <variant>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief A set of ZMW hole numbers or identities to select from a BAM file.
///
/// Constructed from bare hole numbers (matches any read group) or from
/// (rgId, zmw) pairs (exact identity matching). Resolve against a ZmwIndex
/// to get sorted, deduplicated virtual offsets.
class ZmwWhitelist
{
public:
    /// \brief Construct from bare hole numbers (matches any read group).
    explicit ZmwWhitelist(std::vector<std::int32_t> zmws);

    /// \brief Construct from (rgId, zmw) pairs (exact identity matching).
    explicit ZmwWhitelist(std::vector<ZmwIdentity> ids);

    /// \brief Resolve against a ZmwIndex.
    /// \returns virtual offsets sorted in file order, deduplicated.
    std::vector<std::int64_t> Resolve(const ZmwIndex& index) const;

private:
    std::variant<std::vector<std::int32_t>, std::vector<ZmwIdentity>> data_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_INDEX_ZMWWHITELIST_HPP
