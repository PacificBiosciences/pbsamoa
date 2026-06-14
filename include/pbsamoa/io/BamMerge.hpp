#ifndef PBSAMOA_IO_BAMMERGE_HPP
#define PBSAMOA_IO_BAMMERGE_HPP

#include <pbsamoa/io/BamSort.hpp>  // SortOrder

#include <array>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief Configuration for merging already-sorted BAM files.
struct MergeConfig
{
    /// Sort order the inputs are sorted by (and the order of the merged output).
    SortOrder Order{SortOrder::COORDINATE};

    /// Two-char tag the inputs are sorted by; required iff Order == TAG.
    std::array<char, 2> Tag{};

    /// 0 -> auto (min(hardware_concurrency, 8)); drives output BGZF compression.
    std::size_t NumThreads{0};

    /// Output BGZF compression level [1, 12].
    int CompressionLevel{6};

    /// Command line recorded in the appended @PG CL field (CLI passes argv).
    std::optional<std::string> CommandLine{};
};

/// \brief Outcome statistics from a merge.
struct MergeStats
{
    std::int64_t NumRecords{0};  ///< total records written
    std::size_t NumInputs{0};    ///< number of input files merged
};

/// \brief Merge several already-sorted BAM files into one sorted BAM file.
///
/// Streaming k-way merge of the inputs (which must already be sorted by
/// \p config.Order). All inputs must share an identical @SQ reference list; the
/// merged header unions their @RG / @PG / @CO records. Equal-key records are
/// emitted in input-file order (first input wins ties), which is deterministic.
/// Output is written atomically.
///
/// \throws std::runtime_error on I/O errors, a missing input, mismatched sort
///         order, or incompatible reference lists (mirrors BamWriter / SortBam).
MergeStats MergeBam(const std::vector<std::filesystem::path>& inputs,
                    const std::filesystem::path& output, const MergeConfig& config = {});

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_IO_BAMMERGE_HPP
