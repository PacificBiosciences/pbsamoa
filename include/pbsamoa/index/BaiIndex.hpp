#ifndef PBSAMOA_INDEX_BAIINDEX_HPP
#define PBSAMOA_INDEX_BAIINDEX_HPP

#include <pbsamoa/core/Bgzf.hpp>

#include <filesystem>
#include <map>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief A chunk of BAM data: a half-open range of virtual offsets.
struct Chunk
{
    VirtualOffset Begin;
    VirtualOffset End;

    /// \brief Check whether two chunks overlap in virtual offset space.
    bool Overlaps(const Chunk& other) const;
};

/// \brief Per-reference BAI index data: bins and linear index.
struct ReferenceIndex
{
    std::map<std::uint32_t, std::vector<Chunk>> bins;
    std::vector<VirtualOffset> linearIndex;
};

/// \brief BAI index: read, write, build, and query.
class BaiIndex
{
public:
    BaiIndex() = default;

    /// \brief Load a BAI index from a file.
    /// \throws std::runtime_error on bad magic, truncation, or I/O error
    static BaiIndex FromFile(const std::filesystem::path& path);

    /// \brief Write a BAI index to a file.
    /// \throws std::runtime_error on I/O error
    void ToFile(const std::filesystem::path& path) const;

    /// \brief Build a BAI index by scanning a coordinate-sorted BAM file.
    /// \throws std::runtime_error if the BAM header or record stream is not coordinate-sorted
    static BaiIndex Build(const std::filesystem::path& bamPath);

    /// \brief Query the index for chunks overlapping [beg, end) on refId.
    ///
    /// Returns an empty vector if refId is out of range or no records overlap.
    /// Chunks are sorted by Begin and merged when adjacent/overlapping.
    std::vector<Chunk> Query(std::int32_t refId, std::int32_t beg, std::int32_t end) const;

    std::int32_t NumReferences() const;
    const ReferenceIndex& Reference(std::int32_t refId) const;
    std::uint64_t MappedCount() const;
    std::uint64_t UnmappedCount() const;

private:
    std::vector<ReferenceIndex> references_;
    std::uint64_t mappedCount_{0};
    std::uint64_t unmappedCount_{0};
    std::uint64_t noCoorCount_{0};  // n_no_coor: unplaced (refId < 0) reads only
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_INDEX_BAIINDEX_HPP
