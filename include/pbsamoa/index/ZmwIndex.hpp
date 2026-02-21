#ifndef PBSAMOA_INDEX_ZMWINDEX_HPP
#define PBSAMOA_INDEX_ZMWINDEX_HPP

#include <filesystem>
#include <mutex>
#include <span>
#include <unordered_map>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief Identifies a unique ZMW as (readGroupId, holeNumber).
struct ZmwIdentity
{
    std::int32_t rgId{0};
    std::int32_t zmw{0};

    bool operator==(const ZmwIdentity&) const = default;
};

/// \brief In-memory ZMW index loaded from .zmi or .pbi files.
///
/// Stores (rgId, zmw, virtualOffset) per BAM record in file order.
/// Supports point queries, batch queries, and chunking via UniqueZmws().
class ZmwIndex
{
public:
    ZmwIndex() = default;
    ZmwIndex(const ZmwIndex&) = delete;
    ZmwIndex& operator=(const ZmwIndex&) = delete;
    ZmwIndex(ZmwIndex&& other) noexcept;
    ZmwIndex& operator=(ZmwIndex&& other) noexcept;
    ~ZmwIndex() = default;

    /// \brief Load from a .zmi file.
    static ZmwIndex FromZmi(const std::filesystem::path& path);

    /// \brief Load from a legacy .pbi file (extracts rgId, holeNumber, fileOffset).
    static ZmwIndex FromPbi(const std::filesystem::path& path);

    /// \brief Auto-detect: tries .zmi first, falls back to .pbi.
    /// \param[in] bamPath path to the BAM file (index path derived from it)
    static ZmwIndex Open(const std::filesystem::path& bamPath);

    // --- Point / batch queries ---

    /// \brief Find all virtual offsets for records matching a ZMW hole number.
    /// Returns offsets from all read groups.
    std::vector<std::int64_t> Find(std::int32_t zmw) const;

    /// \brief Find all virtual offsets for records matching (rgId, zmw).
    std::vector<std::int64_t> Find(ZmwIdentity id) const;

    /// \brief Find all virtual offsets for a set of ZMW identities.
    std::vector<std::int64_t> Find(std::span<const ZmwIdentity> ids) const;

    // --- Chunking support ---

    /// \brief All unique (rgId, zmw) pairs in file order.
    std::vector<ZmwIdentity> UniqueZmws() const;

    /// \brief Virtual offset of the first record for a given ZMW.
    /// \throws std::runtime_error if ZMW not found
    std::int64_t FirstOffset(std::int32_t zmw) const;

    /// \brief Virtual offset of the first record for a given (rgId, zmw).
    /// \throws std::runtime_error if identity not found
    std::int64_t FirstOffset(ZmwIdentity id) const;

    /// \brief Total number of records in the index.
    std::uint64_t NumRecords() const;

    /// \brief Number of unique (rgId, zmw) identities.
    std::uint64_t NumZmws() const;

private:
    /// \brief Lazily build hash indexes on first lookup.
    void BuildIndex() const;

    /// \brief Combine (rgId, zmw) into a single 64-bit key for hashing.
    static std::uint64_t IdentityKey(std::int32_t rgId, std::int32_t zmw);

    std::vector<std::int32_t> rgIds_;
    std::vector<std::int32_t> zmws_;
    std::vector<std::int64_t> offsets_;

    // Lazy hash indexes for O(1) lookups (built on first query)
    mutable std::unordered_map<std::int32_t, std::vector<std::ptrdiff_t>> zmwIndex_;
    mutable std::unordered_map<std::uint64_t, std::vector<std::ptrdiff_t>> identityIndex_;
    mutable std::uint64_t cachedNumZmws_{0};
    mutable std::once_flag indexOnce_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_INDEX_ZMWINDEX_HPP
