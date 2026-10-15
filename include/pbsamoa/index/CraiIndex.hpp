#ifndef PBSAMOA_INDEX_CRAIINDEX_HPP
#define PBSAMOA_INDEX_CRAIINDEX_HPP

#include <filesystem>
#include <unordered_map>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief One row from a CRAI index.
struct CraiEntry
{
    std::int32_t SequenceId{};
    std::int64_t AlignmentStart{};
    std::int64_t AlignmentSpan{};
    std::int64_t ContainerOffset{};
    std::int64_t SliceOffset{};
    std::int64_t SliceSize{};
};

/// \brief CRAI index: read and query slice-level index entries.
class CraiIndex
{
public:
    CraiIndex() = default;

    /// \brief Load a CRAI index from a gzip-compressed file.
    /// \throws std::runtime_error on parse or I/O errors
    static CraiIndex FromFile(const std::filesystem::path& path);

    /// \brief Return entries for a specific reference sequence.
    /// Includes unmapped rows when refId is -1.
    [[nodiscard]] std::vector<CraiEntry> EntriesForReference(std::int32_t refId) const;

    /// \brief The rows of the slices that overlap [beg, end) of refId (0-based, end exclusive),
    ///        one row for each slice, in file order (container offset, then slice offset).
    ///
    /// A multi-reference slice has a row for each reference: the result keeps one of them.
    /// For refId -1, the rows of the unmapped records.
    [[nodiscard]] std::vector<CraiEntry> Query(std::int32_t refId, std::int64_t beg,
                                               std::int64_t end) const;

private:
    void AddEntry(CraiEntry entry);

    std::unordered_map<std::int32_t, std::vector<CraiEntry>> entries_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_INDEX_CRAIINDEX_HPP
