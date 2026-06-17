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

private:
    void AddEntry(CraiEntry entry);

    std::unordered_map<std::int32_t, std::vector<CraiEntry>> entries_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_INDEX_CRAIINDEX_HPP
