#ifndef PBSAMOA_INDEX_BAIINDEX_HPP
#define PBSAMOA_INDEX_BAIINDEX_HPP

#include <pbsamoa/core/Bgzf.hpp>

#include <filesystem>
#include <map>
#include <memory>
#include <span>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief A chunk of BAM data: a half-open range of virtual offsets.
struct Chunk
{
    VirtualOffset Begin;
    VirtualOffset End;
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
    /// \param[in] bamPath    coordinate-sorted BAM to index
    /// \param[in] numWorkers parallel BGZF-inflate workers; 0 or 1 = single-threaded.
    ///            The index is byte-identical regardless of numWorkers (deterministic).
    /// \throws std::runtime_error if the BAM header or record stream is not coordinate-sorted
    static BaiIndex Build(const std::filesystem::path& bamPath, std::size_t numWorkers = 0);

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

    friend class BaiStreamBuilder;
};

/// \brief Incrementally builds a BAI index from records observed in coordinate order.
///
/// This is the shared indexing core. BaiIndex::Build() drives it from a file scan;
/// the BAM writer's index callback drives it on the fly while merging or sorting.
/// Both paths feed identical (begin offset, record bytes) pairs, so they produce
/// byte-identical output.
///
/// A record's end virtual offset equals the begin offset of the record that follows
/// it, so each record is finalized when the next one arrives; Finalize() supplies the
/// end offset of the last record (the position just past it).
class BaiStreamBuilder
{
public:
    /// \param[in] numReferences @SQ count of the indexed BAM (sizes the reference table)
    explicit BaiStreamBuilder(std::int32_t numReferences);
    ~BaiStreamBuilder();

    BaiStreamBuilder(const BaiStreamBuilder&) = delete;
    BaiStreamBuilder& operator=(const BaiStreamBuilder&) = delete;
    BaiStreamBuilder(BaiStreamBuilder&&) noexcept;
    BaiStreamBuilder& operator=(BaiStreamBuilder&&) noexcept;

    /// \brief Observe one record, in coordinate (write) order.
    /// \param[in] recordBeginVo virtual offset of the record's block_size field
    /// \param[in] recordBytes   BAM record payload WITHOUT the 4-byte block_size prefix
    /// \throws std::runtime_error if the records are not coordinate-sorted
    void Observe(VirtualOffset recordBeginVo, std::span<const std::byte> recordBytes);

    /// \brief Finish the index. Must be called exactly once, after the last Observe().
    /// \param[in] endVo virtual offset just past the last observed record
    BaiIndex Finalize(VirtualOffset endVo);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_INDEX_BAIINDEX_HPP
