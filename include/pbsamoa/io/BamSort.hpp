#ifndef PBSAMOA_IO_BAMSORT_HPP
#define PBSAMOA_IO_BAMSORT_HPP

#include <pbsamoa/io/BamRawReader.hpp>  // ByteLimit, Literals

#include <array>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief Record ordering produced by SortBam / `pbsamoa sort`.
enum class SortOrder
{
    COORDINATE,  ///< by (refId, pos, strand); unmapped last. @HD SO:coordinate
    QUERY_NAME,  ///< natural/numeric read-name order. @HD SO:queryname
    TAG,         ///< by a 2-char auxiliary tag, coordinate fallback. @HD SO:unknown
};

/// \brief Configuration for a native external (RAM-bounded) BAM sort.
struct SortConfig
{
    SortOrder Order{SortOrder::COORDINATE};

    /// Two-char tag to sort by; required iff Order == Tag, otherwise ignored.
    std::array<char, 2> Tag{};

    /// Total in-flight record-buffer budget. Exactly one run is in flight at a
    /// time, so this is the peak run-buffer RAM (not multiplied per thread).
    /// Default 768 MiB (matches samtools' per-thread default, here as the total).
    ByteLimit MaxMemory{std::size_t{768} * 1024 * 1024};

    /// 0 -> auto (min(hardware_concurrency, 8)). Drives BGZF decompress/compress.
    std::size_t NumThreads{0};

    /// Final output BGZF compression level [1, 12]. Temp runs are forced to 1.
    int CompressionLevel{6};

    /// Directory for temporary run files. Default: directory of the output file.
    std::optional<std::filesystem::path> TempDir{};

    /// Cluster unmapped reads by their forward-strand minimiser, equivalent to
    /// `samtools sort -M -R`. Requires coordinate order.
    bool Minimise{false};

    /// Command line recorded in the appended @PG CL field (CLI passes argv).
    /// Spec lists @PG CL but the locked SortConfig had no carrier; added here so
    /// the library can populate it without re-deriving the invocation.
    std::optional<std::string> CommandLine{};

    /// Write the sorted output straight to \p output instead of via a temporary
    /// file and atomic rename. Required for destinations that cannot be renamed
    /// onto — `/dev/stdout`, `/dev/fd/*`, or a FIFO — which a filesystem probe
    /// alone cannot distinguish (a redirected `/dev/stdout` stats as the regular
    /// file it points to). The caller, which knows the destination is a stream,
    /// sets this. Spilled runs then fall back to the current directory when no
    /// TempDir is given. Existing non-regular outputs are also auto-detected.
    bool DirectToOutput{false};
};

/// \brief Outcome statistics from a sort.
struct SortStats
{
    std::int64_t NumRecords{0};  ///< total records written
    std::size_t NumRuns{0};      ///< spilled run files (0 => sorted fully in memory)
};

/// \brief Sort a single BAM file into a single sorted BAM file.
///
/// External merge sort bounded by \p config.MaxMemory: builds RAM-bounded runs of
/// RawRecords, stable-sorts and spills each, then k-way merges to the output in
/// bounded fan-in passes so the open run files stay within the process file-descriptor
/// limit regardless of how many runs the memory budget produces.
/// Output is validly sorted (correct order, unmapped-last, deterministic stable
/// tie-breaking); record content matches `samtools sort`, BGZF bytes may differ.
///
/// \throws std::runtime_error on I/O or invariant errors (mirrors BamWriter).
SortStats SortBam(const std::filesystem::path& input, const std::filesystem::path& output,
                  const SortConfig& config = {});

namespace detail {

/// \brief Natural/numeric string comparison replicating htslib `strnum_cmp`.
///
/// Compares run-by-run: digit runs compare by numeric value (leading zeros
/// ignored, longer number larger, leading-zero count as final tiebreak); other
/// characters compare by byte. Returns <0, 0, or >0. Exposed for unit testing.
int StrNumCmp(std::string_view a, std::string_view b);

}  // namespace detail

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_IO_BAMSORT_HPP
