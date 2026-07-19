#include <pbsamoa/io/BamMerge.hpp>
#include <pbsamoa/io/BamSort.hpp>

#include "BinaryUtils.hpp"
#include "WriterUtils.hpp"

#include <pbsamoa/core/Bgzf.hpp>
#include <pbsamoa/core/Endian.hpp>
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/core/Tags.hpp>
#include <pbsamoa/index/BaiIndex.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamWriter.hpp>
#include <pbsamoa/io/MergeReadAhead.hpp>

#include <pbsamoa/PbSamoaLibraryInfo.hpp>

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <format>
#include <fstream>
#include <ios>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <thread>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

#include <cstddef>
#include <cstdint>
#include <ctime>

namespace PacBio {
namespace Samoa {

namespace detail {

namespace {

constexpr bool IsAsciiDigit(char c) { return (c >= '0') && (c <= '9'); }

}  // namespace

int StrNumCmp(std::string_view a, std::string_view b)
{
    const std::size_t na{std::size(a)};
    const std::size_t nb{std::size(b)};
    std::size_t ia{0};
    std::size_t ib{0};

    while ((ia < na) && (ib < nb)) {
        if (IsAsciiDigit(a[ia]) && IsAsciiDigit(b[ib])) {
            // Skip leading zeros in each number independently.
            while ((ia < na) && (a[ia] == '0')) {
                ++ia;
            }
            while ((ib < nb) && (b[ib] == '0')) {
                ++ib;
            }
            // Advance over the common, equal-valued digit prefix.
            while ((ia < na) && (ib < nb) && IsAsciiDigit(a[ia]) && IsAsciiDigit(b[ib]) &&
                   (a[ia] == b[ib])) {
                ++ia;
                ++ib;
            }
            const bool da{(ia < na) && IsAsciiDigit(a[ia])};
            const bool db{(ib < nb) && IsAsciiDigit(b[ib])};
            if (da && db) {
                // Both diverge mid-digit-run: the number with more remaining
                // digits is the larger value; equal-length runs compare by digit.
                std::size_t k{0};
                while (((ia + k) < na) && ((ib + k) < nb) && IsAsciiDigit(a[ia + k]) &&
                       IsAsciiDigit(b[ib + k])) {
                    ++k;
                }
                const bool moreA{((ia + k) < na) && IsAsciiDigit(a[ia + k])};
                const bool moreB{((ib + k) < nb) && IsAsciiDigit(b[ib + k])};
                if (moreA != moreB) {
                    return moreA ? 1 : -1;
                }
                return (static_cast<unsigned char>(a[ia]) < static_cast<unsigned char>(b[ib])) ? -1
                                                                                               : 1;
            }
            if (da != db) {
                return da ? 1 : -1;  // trailing digits => larger number
            }
            // Both numbers ended with equal value; fewer consumed characters
            // (i.e. fewer leading zeros) sorts first.
            if (ia != ib) {
                return (ia < ib) ? -1 : 1;
            }
            // Identical so far — fall through and continue with the remainder.
        } else {
            if (a[ia] != b[ib]) {
                return (static_cast<unsigned char>(a[ia]) < static_cast<unsigned char>(b[ib])) ? -1
                                                                                               : 1;
            }
            ++ia;
            ++ib;
        }
    }

    if (ia < na) {
        return 1;
    }
    if (ib < nb) {
        return -1;
    }
    return 0;
}

}  // namespace detail

namespace {

/// CPU time consumed by the calling thread so far, in nanoseconds. The merge runs
/// its heap on one consumer thread; measuring that thread's own CPU (vs inferring
/// it from wall minus waits) separates real merge compute from time the thread sat
/// off-core because the compress pool oversubscribed the machine. POSIX
/// CLOCK_THREAD_CPUTIME_ID is available on both Linux and macOS; on the unexpected
/// failure path it returns 0, yielding a harmless zero delta.
std::uint64_t ThreadCpuNanos()
{
    timespec ts{};
    if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &ts) != 0) {
        return 0;
    }
    return (static_cast<std::uint64_t>(ts.tv_sec) * 1'000'000'000ULL) +
           static_cast<std::uint64_t>(ts.tv_nsec);
}

/// Per-sort comparison context: shared by run-build, spill, and merge.
struct SortContext
{
    SortOrder Order{SortOrder::COORDINATE};
    TagKey Tag{};
    std::int32_t NumReferences{0};
    bool Minimise{false};
};

/// A record paired with its precomputed comparison keys. Computed once per
/// record (never per comparison): `CoordKey` always; `TagSortValue` only for tag
/// order (engaged == tag present; nullopt == tag missing, which sorts first).
struct Item
{
    RawRecord Record;
    std::uint64_t CoordKey{0};
    std::uint64_t MinHash{0};
    std::int32_t MinHashPosition{0};
    std::optional<TagValue> TagSortValue{};
};

std::pair<std::uint64_t, std::int32_t> MinimiserKey(const RawRecord& record)
{
    constexpr std::uint32_t KMER{20};
    constexpr std::int32_t KMER_OFFSET{19};
    constexpr std::uint64_t MASK{(std::uint64_t{1} << (2 * KMER)) - 1};
    constexpr std::uint64_t XOR{0xdead7878beef7878ULL};

    const SequenceView sequence{record.Seq()};
    const auto base = [](const char nucleotide) {
        switch (nucleotide) {
            case 'C':
                return std::uint64_t{1};
            case 'G':
                return std::uint64_t{2};
            case 'T':
                return std::uint64_t{3};
            default:
                return std::uint64_t{0};
        }
    };

    std::uint64_t hash{0};
    std::uint64_t minimum{std::numeric_limits<std::uint64_t>::max()};
    std::int32_t minimumPosition{0};
    std::uint32_t i{0};
    for (; (i < (KMER - 1)) && (i < sequence.Size()); ++i) {
        hash = (hash << 2U) | base(sequence[i]);
    }
    for (; i < sequence.Size(); ++i) {
        hash = (hash << 2U) | base(sequence[i]);
        const std::uint64_t candidate{(hash ^ XOR) & MASK};
        if (candidate < minimum) {
            minimum = candidate;
            // i indexes into the sequence, whose length is BAM l_qseq (std::int32_t),
            // so it never exceeds INT32_MAX.
            minimumPosition = static_cast<std::int32_t>(i);
        }
    }

    minimum += std::uint64_t{1} << 30U;
    minimumPosition -= KMER_OFFSET;
    return {minimum, minimumPosition <= 65535 ? 65535 - minimumPosition : 0};
}

/// samtools coordinate key: (refAdj << 32) | ((pos+1) << 1) | rev, with the
/// unmapped reference id mapped to NumReferences so unmapped records sort last
/// and forward precedes reverse at equal position.
std::uint64_t CoordinateSortKey(const RawRecord& record, std::int32_t numReferences)
{
    const std::int32_t refId{record.RefId()};
    const std::uint64_t refAdj{(refId < 0) ? static_cast<std::uint64_t>(numReferences)
                                           : static_cast<std::uint64_t>(refId)};
    const std::int64_t pos{record.Pos()};  // 0-based; -1 for unmapped/unplaced
    const std::uint64_t posPart{static_cast<std::uint64_t>(pos + 1)};
    const std::uint64_t rev{record.IsReverseStrand() ? 1U : 0U};
    return (refAdj << 32U) | (posPart << 1U) | rev;
}

Item MakeItem(RawRecord record, const SortContext& ctx)
{
    const std::uint64_t coordKey{CoordinateSortKey(record, ctx.NumReferences)};
    std::uint64_t minHash{0};
    std::int32_t minHashPosition{0};
    if (ctx.Minimise && (record.RefId() < 0)) {
        std::tie(minHash, minHashPosition) = MinimiserKey(record);
    }
    std::optional<TagValue> tagValue{};
    if (ctx.Order == SortOrder::TAG) {
        const TagMap tags{record.ParseTags()};
        if (const TagValue* value{tags.Get(ctx.Tag)}; value != nullptr) {
            tagValue = *value;
        }
    }
    return Item{
        .Record = std::move(record),
        .CoordKey = coordKey,
        .MinHash = minHash,
        .MinHashPosition = minHashPosition,
        .TagSortValue = std::move(tagValue),
    };
}

bool TagAsNumber(const TagValue& value, double& out)
{
    if (const auto* p{std::get_if<std::int64_t>(&value)}) {
        out = static_cast<double>(*p);
        return true;
    }
    if (const auto* p{std::get_if<float>(&value)}) {
        out = static_cast<double>(*p);
        return true;
    }
    return false;
}

bool TagAsText(const TagValue& value, std::string& out)
{
    if (const auto* p{std::get_if<char>(&value)}) {
        out.assign(1, *p);
        return true;
    }
    if (const auto* p{std::get_if<std::string>(&value)}) {
        out = *p;
        return true;
    }
    return false;
}

/// Type-normalized tag comparison: numbers numerically, char/string lexically.
/// Unsupported types (H, B) or mixed numeric/textual pairs return 0, deferring
/// to the coordinate fallback in RecordCompare.
int CompareTagValues(const TagValue& a, const TagValue& b)
{
    double na{0.0};
    double nb{0.0};
    if (TagAsNumber(a, na) && TagAsNumber(b, nb)) {
        if (std::holds_alternative<std::int64_t>(a) && std::holds_alternative<std::int64_t>(b)) {
            const std::int64_t ia{std::get<std::int64_t>(a)};
            const std::int64_t ib{std::get<std::int64_t>(b)};
            return (ia < ib) ? -1 : (ia > ib) ? 1 : 0;
        }
        return (na < nb) ? -1 : (na > nb) ? 1 : 0;
    }

    std::string sa{};
    std::string sb{};
    if (TagAsText(a, sa) && TagAsText(b, sb)) {
        const int c{sa.compare(sb)};
        return (c < 0) ? -1 : (c > 0) ? 1 : 0;
    }

    return 0;
}

/// Three-way record comparison under the active order. Reads precomputed keys;
/// reads Name()/Flag() live (alloc-free) for query-name order.
int RecordCompare(const SortContext& ctx, const Item& a, const Item& b)
{
    if (ctx.Minimise && (a.Record.RefId() < 0) && (b.Record.RefId() < 0)) {
        if (a.MinHash != b.MinHash) {
            return a.MinHash < b.MinHash ? -1 : 1;
        }
        if (a.MinHashPosition != b.MinHashPosition) {
            return a.MinHashPosition > b.MinHashPosition ? -1 : 1;
        }
    }
    switch (ctx.Order) {
        case SortOrder::COORDINATE:
            return (a.CoordKey < b.CoordKey) ? -1 : (a.CoordKey > b.CoordKey) ? 1 : 0;
        case SortOrder::QUERY_NAME: {
            if (const int c{detail::StrNumCmp(a.Record.Name(), b.Record.Name())}; c != 0) {
                return c;
            }
            const std::int32_t fa{a.Record.Flag() & 0xC0};
            const std::int32_t fb{b.Record.Flag() & 0xC0};
            return (fa < fb) ? -1 : (fa > fb) ? 1 : 0;
        }
        case SortOrder::TAG: {
            const bool ha{a.TagSortValue.has_value()};
            const bool hb{b.TagSortValue.has_value()};
            if (ha != hb) {
                return ha ? 1 : -1;  // missing tag sorts first
            }
            if (ha && hb) {
                if (const int c{CompareTagValues(*a.TagSortValue, *b.TagSortValue)}; c != 0) {
                    return c;
                }
            }
            return (a.CoordKey < b.CoordKey)   ? -1
                   : (a.CoordKey > b.CoordKey) ? 1
                                               : 0;  // fallback on equal/both-missing
        }
    }
    return 0;
}

constexpr std::size_t ITEM_SLACK{64};

/// Approximate in-RAM footprint of one buffered record (raw bytes + the separate
/// CIGAR copy RawRecord keeps + per-Item bookkeeping slack).
std::size_t Footprint(const Item& item)
{
    return std::size(item.Record.RawData()) +
           (item.Record.CigarOps().size() * sizeof(std::uint32_t)) + sizeof(Item) + ITEM_SLACK;
}

std::size_t ResolveThreads(std::size_t requested)
{
    if (requested != 0) {
        return requested;
    }
    const std::uint32_t hardware{std::thread::hardware_concurrency()};
    const std::size_t available{(hardware == 0) ? std::size_t{1} : hardware};
    return std::min<std::size_t>(available, 8);
}

/// Removes registered temp files on scope exit — success or exception.
class TempFileCleanup
{
public:
    TempFileCleanup() = default;
    TempFileCleanup(const TempFileCleanup&) = delete;
    TempFileCleanup& operator=(const TempFileCleanup&) = delete;
    TempFileCleanup(TempFileCleanup&&) = delete;
    TempFileCleanup& operator=(TempFileCleanup&&) = delete;

    void Add(std::filesystem::path path) { paths_.push_back(std::move(path)); }

    ~TempFileCleanup()
    {
        for (const std::filesystem::path& path : paths_) {
            std::error_code ec{};
            std::filesystem::remove(path, ec);
        }
    }

private:
    std::vector<std::filesystem::path> paths_;
};

/// A single merge input: either a spilled run file or the in-memory final run.
struct MergeSource
{
    struct MemoryState
    {
        std::vector<Item> Items;
        std::size_t Pos{0};
    };

    std::variant<std::unique_ptr<BamRawReader>, MemoryState> Source;

    std::optional<Item> Next(const SortContext& ctx)
    {
        return std::visit(
            [&ctx](auto& s) -> std::optional<Item> {
                if constexpr (std::is_same_v<std::decay_t<decltype(s)>,
                                             std::unique_ptr<BamRawReader>>) {
                    std::optional<RawRecord> record{s->ReadRecord()};
                    if (!record) {
                        return std::nullopt;
                    }
                    return MakeItem(std::move(*record), ctx);
                } else {
                    if (s.Pos < std::size(s.Items)) {
                        return std::move(s.Items[s.Pos++]);
                    }
                    return std::nullopt;
                }
            },
            Source);
    }
};

std::string SortOrderText(SortOrder order)
{
    switch (order) {
        case SortOrder::COORDINATE:
            return "coordinate";
        case SortOrder::QUERY_NAME:
            return "queryname";
        case SortOrder::TAG:
            return "unknown";
    }
    std::unreachable();
}

/// Append a chained @PG describing this tool invocation (PP = current last @PG).
void AppendProgramRecord(SamHeader& header, std::string id,
                         const std::optional<std::string>& commandLine)
{
    ProgramRecord program{std::move(id)};
    program.SetTag("PN", "pbsamoa");
    program.SetTag("VN", std::string{GetLibraryInfo().Release});
    if (commandLine) {
        program.SetTag("CL", *commandLine);
    }
    if (const std::span<const ProgramRecord> existing{header.ProgramRecords()}; !existing.empty()) {
        program.SetTag("PP", std::string{existing.back().Id()});
    }
    header.AddProgramRecord(std::move(program));
}

/// Upper bound on run files (== open file descriptors) merged in a single pass.
/// Keeps the simultaneously-open run count well under a conservative process
/// open-file limit (RLIMIT_NOFILE, e.g. macOS's default soft limit of 256); see
/// MergeRuns for the multi-pass fallback when there are more runs than this.
constexpr std::size_t MERGE_MAX_FANIN{64};

/// Heap-merge already-open \p sources into \p writer (left open for the caller to
/// close). Returns the number of records written.
///
/// A max-heap whose "max" is the smallest record: invert the record order so the
/// smallest record floats to the top; ties break on lower source index (== earlier
/// input run), reconstructing input order for equal keys.
std::int64_t MergeSourcesIntoWriter(std::vector<MergeSource>& sources, BamWriter& writer,
                                    const SortContext& ctx)
{
    struct HeapEntry
    {
        Item Value;
        std::size_t Source{0};
    };

    const auto heapLess = [&ctx](const HeapEntry& a, const HeapEntry& b) {
        const int c{RecordCompare(ctx, a.Value, b.Value)};
        if (c != 0) {
            return c > 0;
        }
        return a.Source > b.Source;
    };

    std::vector<HeapEntry> heap{};
    heap.reserve(std::size(sources));
    for (std::size_t i{0}; i < std::size(sources); ++i) {
        if (std::optional<Item> first{sources[i].Next(ctx)}) {
            heap.push_back(HeapEntry{.Value = std::move(*first), .Source = i});
        }
    }
    std::ranges::make_heap(heap, heapLess);

    std::int64_t recordsWritten{0};
    while (!heap.empty()) {
        std::ranges::pop_heap(heap, heapLess);
        const HeapEntry entry{std::move(heap.back())};
        heap.pop_back();

        writer.Write(entry.Value.Record);
        ++recordsWritten;

        if (std::optional<Item> next{sources[entry.Source].Next(ctx)}) {
            heap.push_back(HeapEntry{.Value = std::move(*next), .Source = entry.Source});
            std::ranges::push_heap(heap, heapLess);
        }
    }
    return recordsWritten;
}

/// Open each run file as a synchronous (single-fd) reader. Sources are read
/// synchronously (BgzfWorkers = 0): a background decompression pipeline per source
/// does not help, since the heap consumes one record at a time from a single source
/// while the other readers' pipeline threads would busy-spin on their SPSC queues,
/// burning CPU without improving wall time. One fd per source also keeps the
/// open-file count predictable for fan-in budgeting.
std::vector<MergeSource> OpenRunSources(std::span<const std::filesystem::path> runPaths)
{
    std::vector<MergeSource> sources{};
    sources.reserve(std::size(runPaths));
    for (const std::filesystem::path& path : runPaths) {
        sources.push_back(MergeSource{
            std::make_unique<BamRawReader>(path, BamRawReaderConfig{.BgzfWorkers = 0}),
        });
    }
    return sources;
}

/// k-way merge of sorted runs (spilled run files plus an optional already-sorted
/// in-memory run) into \p outPath. Returns the number of records written.
///
/// Each open run holds a file descriptor, so opening every run at once fails once
/// the run count approaches the process open-file limit (RLIMIT_NOFILE) — a tiny
/// --memory budget over a large input can produce thousands of runs. At most
/// MERGE_MAX_FANIN runs are open at any time: when there are more, runs are merged
/// in fixed-size batches into intermediate run files across as many passes as needed
/// (a balanced k-way merge tree) until a final batch fits. Tie-break stability is
/// preserved across passes because batches are contiguous in run order and merged
/// lowest-index-first at every level; the in-memory tail (newest records) is merged
/// last so equal keys keep input order.
std::int64_t MergeRuns(const std::vector<std::filesystem::path>& runPaths,
                       std::vector<Item> memoryRun, const std::filesystem::path& outPath,
                       const SamHeader& outHeader, const BamWriterConfig& outWriterConfig,
                       const SortContext& ctx)
{
    const bool hasMemory{!memoryRun.empty()};

    // Single pass: every run file plus the in-memory tail fit under the fan-in
    // budget, so merge straight to the output and keep the tail in RAM.
    if ((std::size(runPaths) + (hasMemory ? std::size_t{1} : std::size_t{0})) <= MERGE_MAX_FANIN) {
        std::vector<MergeSource> sources{OpenRunSources(runPaths)};
        if (hasMemory) {
            sources.push_back(MergeSource{MergeSource::MemoryState{std::move(memoryRun)}});
        }
        BamWriter writer{outPath, outHeader, outWriterConfig};
        const std::int64_t recordsWritten{MergeSourcesIntoWriter(sources, writer, ctx)};
        writer.Close();
        return recordsWritten;
    }

    // Multi-pass: too many runs to open at once. Intermediate run files are written
    // with cheap compression (they are re-read immediately) and cleaned up here; the
    // caller owns runPaths, which are never removed by this function.
    TempFileCleanup intermediateCleanup{};
    BamWriterConfig passWriterConfig{outWriterConfig};
    passWriterConfig.BgzfConfig.CompressionLevel = 1;

    // Spill the in-memory tail to a file so every merge input is uniform, appended
    // last to preserve its newest-records tie-break position.
    std::vector<std::filesystem::path> level{runPaths};
    if (hasMemory) {
        std::filesystem::path memPath{detail::TemporaryWritePath(outPath)};
        memPath += ".memrun.bam";
        intermediateCleanup.Add(memPath);
        BamWriter writer{memPath, outHeader, passWriterConfig};
        for (const Item& item : memoryRun) {
            writer.Write(item.Record);
        }
        writer.Close();
        level.push_back(std::move(memPath));
    }

    // Reduce the run count by merging fixed-size batches until a single batch remains.
    while (std::size(level) > MERGE_MAX_FANIN) {
        std::vector<std::filesystem::path> next{};
        next.reserve((std::size(level) + MERGE_MAX_FANIN - 1) / MERGE_MAX_FANIN);
        for (std::size_t start{0}; start < std::size(level); start += MERGE_MAX_FANIN) {
            const std::size_t count{std::min(MERGE_MAX_FANIN, std::size(level) - start)};
            std::filesystem::path interPath{detail::TemporaryWritePath(outPath)};
            interPath += ".merge.bam";
            intermediateCleanup.Add(interPath);

            std::vector<MergeSource> sources{
                OpenRunSources(std::span{level}.subspan(start, count))};
            BamWriter writer{interPath, outHeader, passWriterConfig};
            MergeSourcesIntoWriter(sources, writer, ctx);
            writer.Close();

            next.push_back(std::move(interPath));
        }
        level = std::move(next);
    }

    std::vector<MergeSource> sources{OpenRunSources(level)};
    BamWriter writer{outPath, outHeader, outWriterConfig};
    const std::int64_t recordsWritten{MergeSourcesIntoWriter(sources, writer, ctx)};
    writer.Close();
    return recordsWritten;
}

// --- merge-endpoint helpers ---

/// True iff both headers carry an identical @SQ list (same names, lengths, order).
/// Merging requires this so per-record refIds stay valid without remapping.
bool SameReferences(const SamHeader& a, const SamHeader& b)
{
    const std::span<const ReferenceSequence> ra{a.ReferenceSequences()};
    const std::span<const ReferenceSequence> rb{b.ReferenceSequences()};
    if (std::size(ra) != std::size(rb)) {
        return false;
    }
    for (std::size_t i{0}; i < std::size(ra); ++i) {
        if ((ra[i].Name() != rb[i].Name()) || (ra[i].Length() != rb[i].Length())) {
            return false;
        }
    }
    return true;
}

/// Throws unless \p header advertises the sort order the merge expects. Tag order
/// cannot be validated from @HD (inputs carry SO:unknown), so it is accepted.
void ValidateMergeSortOrder(std::string_view actualOrder, SortOrder order,
                            const std::filesystem::path& path)
{
    std::string_view required{};
    if (order == SortOrder::COORDINATE) {
        required = "coordinate";
    } else if (order == SortOrder::QUERY_NAME) {
        required = "queryname";
    } else {
        return;
    }
    if (actualOrder != required) {
        throw std::runtime_error{
            std::format("MergeBam: input {} has @HD SO:{}, expected SO:{} for this merge order",
                        path.string(), actualOrder, required)};
    }
}

void UnionReadGroups(SamHeader& dst, const SamHeader& src)
{
    for (const ReadGroup& group : src.ReadGroups()) {
        const bool present{std::ranges::any_of(dst.ReadGroups(), [&](const ReadGroup& existing) {
            return existing.Id() == group.Id();
        })};
        if (!present) {
            dst.AddReadGroup(group);
        }
    }
}

void UnionProgramRecords(SamHeader& dst, const SamHeader& src)
{
    for (const ProgramRecord& program : src.ProgramRecords()) {
        const bool present{std::ranges::any_of(
            dst.ProgramRecords(),
            [&](const ProgramRecord& existing) { return existing.Id() == program.Id(); })};
        if (!present) {
            dst.AddProgramRecord(program);
        }
    }
}

void UnionComments(SamHeader& dst, const SamHeader& src)
{
    for (const std::string& comment : src.Comments()) {
        const bool present{std::ranges::any_of(
            dst.Comments(), [&](const std::string& existing) { return existing == comment; })};
        if (!present) {
            dst.AddComment(comment);
        }
    }
}

// --- concat (BGZF block passthrough) ---

constexpr std::size_t CONCAT_MAX_BLOCK{65536};                   // largest possible BGZF block
constexpr std::size_t CONCAT_COPY_BUFFER{std::size_t{1} << 20};  // verbatim-copy chunk
constexpr std::streamoff EOF_MARKER_LEN{static_cast<std::streamoff>(std::size(BGZF_EOF_MARKER))};

void WriteAll(std::ofstream& out, std::span<const std::byte> bytes)
{
    out.write(reinterpret_cast<const char*>(std::data(bytes)),
              static_cast<std::streamsize>(std::size(bytes)));
    if (!out) {
        throw std::runtime_error{"ConcatBam: output write failed"};
    }
}

/// Emit \p data as one or more framed BGZF blocks, chunked to the max payload.
void WriteAsBgzfBlocks(std::ofstream& out, std::span<const std::byte> data, const int level)
{
    for (std::size_t offset{0}; offset < std::size(data); offset += BGZF_MAX_UNCOMPRESSED_BLOCK) {
        const std::size_t len{std::min(BGZF_MAX_UNCOMPRESSED_BLOCK, std::size(data) - offset)};
        const std::vector<std::byte> frame{CompressBgzfBlock(data.subspan(offset, len), level)};
        WriteAll(out, frame);
    }
}

/// Byte-level concatenation of \p inputs via BGZF block passthrough: writes
/// \p mergedHeader, then for each input drops its header and copies the compressed
/// record blocks verbatim, decompressing only the single block that straddles the
/// header/record boundary. Identical @SQ across inputs (enforced by the caller)
/// keeps every record's refId valid without remapping. Records are not counted
/// (returns -1) — counting would require decompressing every block, defeating the
/// passthrough.
std::int64_t ConcatBam(const std::vector<std::filesystem::path>& inputs,
                       const std::filesystem::path& outPath, const SamHeader& mergedHeader,
                       const int level)
{
    std::ofstream out{outPath, std::ios::binary};
    if (!out) {
        throw std::runtime_error{std::format("ConcatBam: cannot open output {}", outPath.string())};
    }

    const std::vector<std::byte> headerBytes{mergedHeader.ToBamHeaderBlock()};
    WriteAsBgzfBlocks(out, headerBytes, level);

    for (const std::filesystem::path& input : inputs) {
        std::ifstream in{input, std::ios::binary};
        if (!in) {
            throw std::runtime_error{
                std::format("ConcatBam: cannot open input {}", input.string())};
        }
        in.seekg(0, std::ios::end);
        const std::streamoff fileSize{static_cast<std::streamoff>(in.tellg())};
        in.seekg(0, std::ios::beg);

        // Walk BGZF blocks from the start, decompressing each, until the BAM header
        // ends. The block whose cumulative decompressed size first reaches the header
        // length is the boundary block: its record tail is recompressed, then every
        // remaining compressed block is copied verbatim.
        std::vector<std::byte> headerPrefix{};
        std::size_t cumUncompressed{0};
        std::streamoff verbatimStart{-1};

        while (true) {
            const std::streamoff blockStart{static_cast<std::streamoff>(in.tellg())};
            std::array<std::byte, 18> blockHeader{};
            in.read(reinterpret_cast<char*>(std::data(blockHeader)), 18);
            const std::streamsize headerRead{in.gcount()};
            if (headerRead == 0) {
                break;  // clean EOF: no record region
            }
            if (headerRead < 18) {
                throw std::runtime_error{
                    std::format("ConcatBam: truncated BGZF header in {}", input.string())};
            }
            const std::optional<BgzfBlockInfo> info{ParseBgzfBlockHeader(blockHeader)};
            if (!info) {
                throw std::runtime_error{
                    std::format("ConcatBam: malformed BGZF block in {}", input.string())};
            }

            std::vector<std::byte> block(info->blockSize);
            std::ranges::copy(blockHeader, std::begin(block));
            const std::streamsize bodyBytes{static_cast<std::streamsize>(info->blockSize) - 18};
            in.read(reinterpret_cast<char*>(std::data(block) + 18), bodyBytes);
            if (in.gcount() < bodyBytes) {
                throw std::runtime_error{
                    std::format("ConcatBam: truncated BGZF block in {}", input.string())};
            }
            if (IsBgzfEofMarker(block)) {
                break;  // header-only input: no records to copy
            }

            std::vector<std::byte> decompressed(CONCAT_MAX_BLOCK);
            const std::optional<std::size_t> produced{
                DecompressBgzfBlock(block, *info, decompressed)};
            if (!produced) {
                throw std::runtime_error{
                    std::format("ConcatBam: BGZF decompression failed in {}", input.string())};
            }
            decompressed.resize(*produced);

            const std::size_t cumBefore{cumUncompressed};
            headerPrefix.insert(std::end(headerPrefix), std::begin(decompressed),
                                std::end(decompressed));
            cumUncompressed += *produced;

            const std::size_t headerLen{ComputeHeaderSize(headerPrefix)};
            if (headerLen == 0) {
                continue;  // header continues into the next block
            }

            const std::size_t recordStart{headerLen - cumBefore};
            if (recordStart < *produced) {
                WriteAsBgzfBlocks(out, std::span{decompressed}.subspan(recordStart), level);
            }
            verbatimStart = blockStart + static_cast<std::streamoff>(info->blockSize);
            break;
        }

        if (verbatimStart < 0) {
            continue;  // header-only or empty: nothing to copy
        }

        // Exclude a trailing BGZF EOF marker from the verbatim copy.
        std::streamoff eofStart{fileSize};
        if ((fileSize - verbatimStart) >= EOF_MARKER_LEN) {
            std::array<std::byte, 28> trailer{};
            in.clear();
            in.seekg(fileSize - EOF_MARKER_LEN, std::ios::beg);
            in.read(reinterpret_cast<char*>(std::data(trailer)),
                    static_cast<std::streamsize>(std::size(trailer)));
            if (IsBgzfEofMarker(trailer)) {
                eofStart = fileSize - EOF_MARKER_LEN;
            }
        }

        in.clear();
        in.seekg(verbatimStart, std::ios::beg);
        std::vector<char> buffer(CONCAT_COPY_BUFFER);
        std::streamoff remaining{eofStart - verbatimStart};
        while (remaining > 0) {
            const std::streamsize chunk{static_cast<std::streamsize>(std::min<std::streamoff>(
                remaining, static_cast<std::streamoff>(std::size(buffer))))};
            in.read(std::data(buffer), chunk);
            if (in.gcount() != chunk) {
                throw std::runtime_error{
                    std::format("ConcatBam: truncated record region in {}", input.string())};
            }
            out.write(std::data(buffer), chunk);
            if (!out) {
                throw std::runtime_error{"ConcatBam: output write failed"};
            }
            remaining -= chunk;
        }
    }

    WriteAll(out, BGZF_EOF_MARKER);
    out.flush();
    if (!out) {
        throw std::runtime_error{"ConcatBam: output write failed"};
    }
    return -1;
}

// --- parallel read-ahead k-way merge ---

/// Coordinate sort key computed directly from a record's raw bytes (refID@0,
/// pos@4, FLAG@14), mirroring CoordinateSortKey without materializing a RawRecord.
std::uint64_t CoordinateSortKeyFromBytes(std::span<const std::byte> bytes,
                                         std::int32_t numReferences)
{
    const std::int32_t refId{ReadI32LE(std::data(bytes))};
    const std::uint64_t refAdj{(refId < 0) ? static_cast<std::uint64_t>(numReferences)
                                           : static_cast<std::uint64_t>(refId)};
    const std::int64_t pos{ReadI32LE(std::data(bytes) + 4)};
    const std::uint64_t posPart{static_cast<std::uint64_t>(pos + 1)};
    const std::uint64_t rev{((ReadU16LE(std::data(bytes) + 14) & 0x10U) != 0U) ? 1U : 0U};
    return (refAdj << 32U) | (posPart << 1U) | rev;
}

/// Read name (string_view) from a record's raw bytes: l_read_name@8 (counts the
/// trailing NUL), name@32.
std::string_view NameFromBytes(std::span<const std::byte> bytes)
{
    const std::uint8_t nameLen{static_cast<std::uint8_t>(bytes[8])};
    return std::string_view{reinterpret_cast<const char*>(std::data(bytes) + 32),
                            static_cast<std::size_t>(nameLen - 1)};
}

/// One record in the parallel-merge heap: a view into a live read-ahead buffer plus
/// the precomputed keys (matching Item, but never owning the record bytes).
struct MergeView
{
    std::span<const std::byte> Bytes;
    std::uint64_t CoordKey{0};
    std::uint64_t MinHash{0};
    std::int32_t MinHashPosition{0};
    std::optional<TagValue> TagSortValue{};
};

MergeView MakeMergeView(std::span<const std::byte> bytes, const SortContext& ctx)
{
    std::optional<RawRecord> record{};
    if (ctx.Minimise || (ctx.Order == SortOrder::TAG)) {
        record.emplace(bytes);
    }
    std::uint64_t minHash{0};
    std::int32_t minHashPosition{0};
    if (ctx.Minimise && (record->RefId() < 0)) {
        std::tie(minHash, minHashPosition) = MinimiserKey(*record);
    }
    std::optional<TagValue> tagValue{};
    if (ctx.Order == SortOrder::TAG) {
        // Tag order is rare; materialize a RawRecord only to parse the tag map.
        const TagMap tags{record->ParseTags()};
        if (const TagValue* value{tags.Get(ctx.Tag)}; value != nullptr) {
            tagValue = *value;
        }
    }
    return MergeView{
        .Bytes = bytes,
        .CoordKey = CoordinateSortKeyFromBytes(bytes, ctx.NumReferences),
        .MinHash = minHash,
        .MinHashPosition = minHashPosition,
        .TagSortValue = std::move(tagValue),
    };
}

/// Three-way comparison of two MergeViews under the active order; mirrors
/// RecordCompare but reads name/flag live from the record bytes.
int MergeViewCompare(const SortContext& ctx, const MergeView& a, const MergeView& b)
{
    if (ctx.Minimise && (ReadI32LE(std::data(a.Bytes)) < 0) &&
        (ReadI32LE(std::data(b.Bytes)) < 0)) {
        if (a.MinHash != b.MinHash) {
            return a.MinHash < b.MinHash ? -1 : 1;
        }
        if (a.MinHashPosition != b.MinHashPosition) {
            return a.MinHashPosition > b.MinHashPosition ? -1 : 1;
        }
    }
    switch (ctx.Order) {
        case SortOrder::COORDINATE:
            return (a.CoordKey < b.CoordKey) ? -1 : (a.CoordKey > b.CoordKey) ? 1 : 0;
        case SortOrder::QUERY_NAME: {
            if (const int c{detail::StrNumCmp(NameFromBytes(a.Bytes), NameFromBytes(b.Bytes))};
                c != 0) {
                return c;
            }
            const std::int32_t fa{ReadU16LE(std::data(a.Bytes) + 14) & 0xC0};
            const std::int32_t fb{ReadU16LE(std::data(b.Bytes) + 14) & 0xC0};
            return (fa < fb) ? -1 : (fa > fb) ? 1 : 0;
        }
        case SortOrder::TAG: {
            const bool ha{a.TagSortValue.has_value()};
            const bool hb{b.TagSortValue.has_value()};
            if (ha != hb) {
                return ha ? 1 : -1;  // missing tag sorts first
            }
            if (ha && hb) {
                if (const int c{CompareTagValues(*a.TagSortValue, *b.TagSortValue)}; c != 0) {
                    return c;
                }
            }
            return (a.CoordKey < b.CoordKey) ? -1 : (a.CoordKey > b.CoordKey) ? 1 : 0;
        }
    }
    return 0;
}

/// Streaming k-way merge over the parallel read-ahead: decompression runs on the
/// read-ahead's shared pool, the heap consumes record views (no per-record copy),
/// and the writer compresses in parallel. Equal-key records keep input order via
/// the lower-source-index tie-break, so the output is deterministic regardless of
/// thread count or budget. Returns the number of records written.
std::int64_t MergeReadAheadIntoWriter(MergeReadAhead& readAhead, BamWriter& writer,
                                      const SortContext& ctx)
{
    struct HeapEntry
    {
        MergeView View;
        std::size_t Source{0};
    };

    const auto heapLess = [&ctx](const HeapEntry& a, const HeapEntry& b) {
        const int c{MergeViewCompare(ctx, a.View, b.View)};
        if (c != 0) {
            return c > 0;
        }
        return a.Source > b.Source;
    };

    std::vector<HeapEntry> heap{};
    heap.reserve(readAhead.NumSources());
    for (std::size_t i{0}; i < readAhead.NumSources(); ++i) {
        if (const std::optional<std::span<const std::byte>> bytes{readAhead.Next(i)}) {
            heap.push_back(HeapEntry{.View = MakeMergeView(*bytes, ctx), .Source = i});
        }
    }
    std::ranges::make_heap(heap, heapLess);

    std::int64_t recordsWritten{0};
    while (!heap.empty()) {
        std::ranges::pop_heap(heap, heapLess);
        const std::size_t source{heap.back().Source};
        // The writer copies the bytes synchronously, so the source's next Next()
        // call may safely recycle the buffer this span views into.
        writer.Write(heap.back().View.Bytes);
        ++recordsWritten;

        if (const std::optional<std::span<const std::byte>> bytes{readAhead.Next(source)}) {
            heap.back() = HeapEntry{.View = MakeMergeView(*bytes, ctx), .Source = source};
            std::ranges::push_heap(heap, heapLess);
        } else {
            heap.pop_back();
        }
    }
    return recordsWritten;
}

// --- disjoint-chain verbatim passthrough ---

/// Outcome of probing a coordinate merge for the disjoint-chain fast path.
struct DisjointChainProbe
{
    bool Disjoint{false};
    std::vector<std::size_t> Order{};  ///< input indices, min-sorted, empties excluded
    std::int64_t TotalRecords{0};      ///< exact record count (valid iff Disjoint)
};

/// Probe whether \p inputs form a totally-ordered disjoint coordinate chain, so the
/// merge can emit via verbatim BGZF block passthrough (ConcatBam) instead of the heap
/// merge. Reads each input's first record for its minimum key, sorts indices by
/// (minKey, index), then iterates each input in that order — verifying records are
/// non-decreasing and strictly precede the next input's minimum, while counting.
///
/// The full decompress only completes when the chain is actually disjoint; an
/// overlapping or internally-unsorted input aborts the scan at the first crossing or
/// decreasing record, so the heap fallback pays almost nothing for the probe.
/// Header-only inputs contribute no records and are dropped from the order.
DisjointChainProbe ProbeDisjointChain(const std::vector<std::filesystem::path>& inputs,
                                      std::int32_t numReferences, std::size_t decodeWorkers)
{
    struct Entry
    {
        std::size_t Index;
        std::uint64_t MinKey;
    };

    std::vector<Entry> entries{};
    entries.reserve(std::size(inputs));
    for (std::size_t i{0}; i < std::size(inputs); ++i) {
        BamRawReader reader{inputs[i], BamRawReaderConfig{.BgzfWorkers = 0}};
        auto range{reader.Records()};
        const auto it{range.begin()};
        if (it == range.end()) {
            continue;  // header-only input: no records, no constraint
        }
        entries.push_back(
            Entry{.Index = i, .MinKey = CoordinateSortKeyFromBytes(it->RawData(), numReferences)});
    }

    std::ranges::sort(entries, std::less{},
                      [](const Entry& e) { return std::tuple{e.MinKey, e.Index}; });

    std::int64_t total{0};
    for (std::size_t k{0}; k < std::size(entries); ++k) {
        const bool hasNext{(k + 1) < std::size(entries)};
        const std::uint64_t nextMin{hasNext ? entries[k + 1].MinKey : std::uint64_t{0}};

        BamRawReader reader{inputs[entries[k].Index],
                            BamRawReaderConfig{.BgzfWorkers = decodeWorkers}};
        std::uint64_t prevKey{0};
        bool firstRecord{true};
        for (const RawRecord& record : reader.Records()) {
            const std::uint64_t key{CoordinateSortKeyFromBytes(record.RawData(), numReferences)};
            if (!firstRecord && (key < prevKey)) {
                return DisjointChainProbe{};  // not internally sorted -> heap fallback
            }
            if (hasNext && (key >= nextMin)) {
                return DisjointChainProbe{};  // overlaps the next input -> heap fallback
            }
            prevKey = key;
            firstRecord = false;
            ++total;
        }
    }

    DisjointChainProbe probe{.Disjoint = true, .TotalRecords = total};
    probe.Order.reserve(std::size(entries));
    for (const Entry& entry : entries) {
        probe.Order.push_back(entry.Index);
    }
    return probe;
}

// --- merge-mode resolution ---

enum class MergeMode
{
    SORTED,
    CONCAT,
};

struct MergePlan
{
    MergeMode Mode{MergeMode::SORTED};
    SortOrder Order{SortOrder::COORDINATE};  // meaningful iff Mode == SORTED
};

bool IsUnsortedOrder(std::string_view sortOrder)
{
    return sortOrder.empty() || (sortOrder == "unsorted") || (sortOrder == "unknown");
}

/// Decide how to merge. An explicit \p config.Order forces a validated sorted merge
/// of that order; \p config.Concat forces concat; otherwise the mode is auto-detected
/// from the inputs' @HD SO (all-coordinate / all-queryname / all-unsorted), and mixed
/// orders are rejected.
MergePlan ResolveMergePlan(const MergeConfig& config, const std::vector<std::string>& sortOrders,
                           const std::vector<std::filesystem::path>& inputs)
{
    if (config.Concat) {
        return MergePlan{.Mode = MergeMode::CONCAT, .Order = SortOrder::COORDINATE};
    }
    if (config.Order.has_value()) {
        for (std::size_t i{0}; i < std::size(inputs); ++i) {
            ValidateMergeSortOrder(sortOrders[i], *config.Order, inputs[i]);
        }
        return MergePlan{.Mode = MergeMode::SORTED, .Order = *config.Order};
    }

    if (std::ranges::all_of(sortOrders, [](std::string_view so) { return so == "coordinate"; })) {
        return MergePlan{.Mode = MergeMode::SORTED, .Order = SortOrder::COORDINATE};
    }
    if (std::ranges::all_of(sortOrders, [](std::string_view so) { return so == "queryname"; })) {
        return MergePlan{.Mode = MergeMode::SORTED, .Order = SortOrder::QUERY_NAME};
    }
    if (std::ranges::all_of(sortOrders, IsUnsortedOrder)) {
        return MergePlan{.Mode = MergeMode::CONCAT, .Order = SortOrder::COORDINATE};
    }

    std::string orders{};
    for (std::size_t i{0}; i < std::size(sortOrders); ++i) {
        orders += std::format("{}{}", (i == 0) ? "" : ", ",
                              sortOrders[i].empty() ? "unknown" : sortOrders[i]);
    }
    throw std::runtime_error{std::format(
        "MergeBam: inputs have mixed sort orders ({}); pass --order or --concat", orders)};
}

}  // namespace

SortStats SortBam(const std::filesystem::path& input, const std::filesystem::path& output,
                  const SortConfig& config)
{
    if (!std::filesystem::exists(input)) {
        throw std::runtime_error{
            std::format("SortBam: input file does not exist: {}", input.string())};
    }
    if (config.Minimise && (config.Order != SortOrder::COORDINATE)) {
        throw std::runtime_error{"SortBam: minimiser clustering requires coordinate order"};
    }

    const std::size_t numThreads{ResolveThreads(config.NumThreads)};
    const int finalLevel{std::clamp(config.CompressionLevel, 1, 12)};
    const std::size_t maxMemory{config.MaxMemory.Value()};

    // A stream output (e.g. /dev/stdout or a FIFO) cannot be written through a
    // temp file and atomic rename; the final BamWriter streams straight to it
    // instead, and spilled runs fall back to the current directory when no TempDir
    // is set. The caller flags this explicitly (a redirected /dev/stdout stats as
    // the regular file it points to, so a probe cannot detect it); an existing
    // non-regular output is also auto-detected as a convenience.
    std::error_code outEc{};
    const bool streamOutput{config.DirectToOutput ||
                            (std::filesystem::exists(output, outEc) &&
                             !std::filesystem::is_regular_file(output, outEc))};

    TempFileCleanup cleanup{};
    const std::filesystem::path runBase{
        config.TempDir ? (*config.TempDir / output.filename())
        : streamOutput
            ? detail::TemporaryWritePath(std::filesystem::current_path() / output.filename())
            : detail::TemporaryWritePath(output)};

    SortContext ctx{
        .Order = config.Order,
        .Tag = TagKey{config.Tag[0], config.Tag[1]},
        .NumReferences = 0,
        .Minimise = config.Minimise,
    };
    const auto less = [&ctx](const Item& a, const Item& b) { return RecordCompare(ctx, a, b) < 0; };

    const BamWriterConfig runWriterConfig{
        .BgzfConfig = {.CompressionLevel = 1, .BgzfWorkers = numThreads},
        .UseTempFile = false,
    };

    std::vector<Item> run{};
    std::vector<std::filesystem::path> runPaths{};
    std::int64_t numRecords{0};
    SamHeader outHeader{};

    {
        BamRawReader reader{input, BamRawReaderConfig{.BgzfWorkers = numThreads}};
        const SamHeader& inputHeader{reader.Header()};
        ctx.NumReferences = inputHeader.NumReferences();

        SamHeader runHeader{inputHeader};
        runHeader.SetSortOrder("unsorted");

        std::size_t runBytes{0};
        bool hasMapped{false};
        const auto spill = [&]() {
            std::ranges::stable_sort(run, less);
            std::filesystem::path path{runBase};
            path += std::format(".run{}.bam", std::size(runPaths));
            cleanup.Add(path);
            BamWriter writer{path, runHeader, runWriterConfig};
            for (const Item& item : run) {
                writer.Write(item.Record);
            }
            writer.Close();
            runPaths.push_back(std::move(path));
            run.clear();
            runBytes = 0;
        };

        while (std::optional<RawRecord> record{reader.ReadRecord()}) {
            hasMapped = hasMapped || (record->RefId() >= 0);
            Item item{MakeItem(std::move(*record), ctx)};
            const std::size_t footprint{Footprint(item)};
            if (!run.empty() && ((runBytes + footprint) > maxMemory)) {
                spill();
            }
            runBytes += footprint;
            run.push_back(std::move(item));
            ++numRecords;
        }

        outHeader = inputHeader;
        if (config.Minimise) {
            outHeader.SetSortOrder(hasMapped ? "coordinate" : "unsorted");
            outHeader.SetSubSort(hasMapped ? "coordinate:minhash" : "unsorted:minhash");
        } else {
            outHeader.SetSortOrder(SortOrderText(config.Order));
            outHeader.SetSubSort("");  // drop any stale sub-sort; this order has none
        }
        outHeader.SetGroupOrder("");  // drop any stale group-order
        AppendProgramRecord(outHeader, "pbsamoa.sort", config.CommandLine);
    }  // input reader closed here; its BGZF worker threads are released

    // The remaining `run` is the final (unspilled) run for both paths below.
    std::ranges::stable_sort(run, less);

    // Stream output is written in place; a regular file is written to a temp path
    // and atomically renamed so a partial sort never clobbers the destination.
    const std::filesystem::path writePath{streamOutput ? output
                                                       : detail::TemporaryWritePath(output)};
    if (!streamOutput) {
        cleanup.Add(writePath);
    }
    const BamWriterConfig outWriterConfig{
        .BgzfConfig = {.CompressionLevel = finalLevel, .BgzfWorkers = numThreads},
        .UseTempFile = false,
    };

    if (runPaths.empty()) {
        // Fast path: whole input fit the budget — one run, no merge.
        BamWriter writer{writePath, outHeader, outWriterConfig};
        for (const Item& item : run) {
            writer.Write(item.Record);
        }
        writer.Close();
    } else {
        MergeRuns(runPaths, std::move(run), writePath, outHeader, outWriterConfig, ctx);
    }

    if (!streamOutput) {
        detail::AtomicRename(writePath, output, "sort");
    }

    return SortStats{
        .NumRecords = numRecords,
        .NumRuns = std::size(runPaths),
    };
}

MergeStats MergeBam(const std::vector<std::filesystem::path>& inputs,
                    const std::filesystem::path& output, const MergeConfig& config)
{
    if (inputs.empty()) {
        throw std::runtime_error{"MergeBam: no input files given"};
    }
    for (const std::filesystem::path& input : inputs) {
        if (!std::filesystem::exists(input)) {
            throw std::runtime_error{
                std::format("MergeBam: input file does not exist: {}", input.string())};
        }
    }

    const auto wallStart{std::chrono::steady_clock::now()};
    const auto nsToSeconds{[](std::uint64_t ns) { return static_cast<double>(ns) * 1e-9; }};

    const std::size_t numThreads{ResolveThreads(config.NumThreads)};
    // Decode and compress pools are sized independently; each falls back to the
    // shared resolved thread count when its knob is left at 0 (auto).
    const std::size_t decodeThreads{(config.DecodeThreads != 0) ? config.DecodeThreads
                                                                : numThreads};
    const std::size_t compressThreads{(config.CompressThreads != 0) ? config.CompressThreads
                                                                    : numThreads};
    const int finalLevel{std::clamp(config.CompressionLevel, 1, 12)};

    // Build the merged header from the first input, requiring every input to share
    // an identical @SQ list; union @RG / @PG / @CO. Collect each input's @HD SO so
    // the merge mode (sorted vs concat) can be resolved.
    SamHeader merged{};
    std::vector<std::string> sortOrders{};
    sortOrders.reserve(std::size(inputs));
    {
        const BamRawReader reader{inputs.front(), BamRawReaderConfig{.BgzfWorkers = 0}};
        merged = reader.Header();
        sortOrders.emplace_back(reader.Header().SortOrder());
    }
    for (std::size_t i{1}; i < std::size(inputs); ++i) {
        const BamRawReader reader{inputs[i], BamRawReaderConfig{.BgzfWorkers = 0}};
        const SamHeader& header{reader.Header()};
        if (!SameReferences(merged, header)) {
            throw std::runtime_error{
                std::format("MergeBam: input {} has @SQ references incompatible with {} "
                            "(identical reference lists are required)",
                            inputs[i].string(), inputs.front().string())};
        }
        UnionReadGroups(merged, header);
        UnionProgramRecords(merged, header);
        UnionComments(merged, header);
        sortOrders.emplace_back(header.SortOrder());
    }

    const MergePlan plan{ResolveMergePlan(config, sortOrders, inputs)};

    // BAI can only index a coordinate-sorted, record-decoded output. Reject the
    // combinations that cannot produce one up front; the disjoint-chain passthrough
    // (also undecodable) is caught below, once its probe has run.
    if (config.BaiOutput) {
        if (plan.Mode == MergeMode::CONCAT) {
            throw std::runtime_error{
                "MergeBam: --bai cannot index --concat output (it is unsorted)"};
        }
        if (plan.Order != SortOrder::COORDINATE) {
            throw std::runtime_error{
                std::format("MergeBam: --bai requires a coordinate-sorted merge, not SO:{}",
                            SortOrderText(plan.Order))};
        }
    }

    merged.SetGroupOrder("");
    TempFileCleanup cleanup{};
    const std::filesystem::path finalTmp{detail::TemporaryWritePath(output)};
    cleanup.Add(finalTmp);
    std::optional<std::filesystem::path> baiTmp{};
    if (config.BaiOutput) {
        baiTmp = detail::TemporaryWritePath(*config.BaiOutput);
        cleanup.Add(*baiTmp);
    }

    std::int64_t numRecords{0};
    bool passthrough{false};
    MergeRuntimeStats runtime{};
    runtime.DecodeThreads = decodeThreads;
    runtime.CompressThreads = compressThreads;
    if (plan.Mode == MergeMode::CONCAT) {
        // Concat preserves each input's records in order but interleaves files, so the
        // output is not globally sorted.
        merged.SetSortOrder("unsorted");
        AppendProgramRecord(merged, "pbsamoa.merge", config.CommandLine);
        numRecords = ConcatBam(inputs, finalTmp, merged, finalLevel);
    } else {
        merged.SetSortOrder(SortOrderText(plan.Order));
        AppendProgramRecord(merged, "pbsamoa.merge", config.CommandLine);

        // Coordinate inputs occupying strictly disjoint coordinate ranges need no
        // interleaving: emit them via verbatim BGZF block passthrough (no
        // recompression), reordered by minimum coordinate. Overlapping, internally
        // unsorted, or non-coordinate merges fall back to the heap merge.
        DisjointChainProbe probe{};
        if (plan.Order == SortOrder::COORDINATE) {
            probe = ProbeDisjointChain(inputs, merged.NumReferences(), decodeThreads);
        }

        if (probe.Disjoint) {
            if (config.BaiOutput) {
                throw std::runtime_error{
                    "MergeBam: --bai is not supported on the disjoint-chain passthrough (records "
                    "are copied verbatim, not decoded); omit --bai and run bai-build on the "
                    "output"};
            }
            std::vector<std::filesystem::path> ordered{};
            ordered.reserve(std::size(probe.Order));
            for (const std::size_t idx : probe.Order) {
                ordered.push_back(inputs[idx]);
            }
            ConcatBam(ordered, finalTmp, merged, finalLevel);
            numRecords = probe.TotalRecords;
            passthrough = true;
        } else {
            const SortContext ctx{
                .Order = plan.Order,
                .Tag = TagKey{config.Tag[0], config.Tag[1]},
                .NumReferences = merged.NumReferences(),
            };
            BamWriterConfig outWriterConfig{
                .BgzfConfig = {.CompressionLevel = finalLevel, .BgzfWorkers = compressThreads},
                .UseTempFile = false,
            };
            if (config.WriterQueueCapacity != 0) {
                outWriterConfig.BgzfConfig.InputQueueCapacity = config.WriterQueueCapacity;
            }

            // Parallel k-way merge: every input is decompressed ahead on a shared pool
            // and framed into per-source batches, while the heap consumes zero-copy
            // record views and the writer compresses in parallel.
            MergeReadAhead readAhead{inputs, decodeThreads, config.ReadAheadMemory,
                                     config.BatchBytes};

            // On-the-fly BAI: the writer's index callback feeds each record's begin
            // offset and bytes to the builder from the IO thread; Finalize() runs after
            // Close() joins that thread, so the builder is never touched concurrently.
            // baiBuilder outlives writer (declared first), so the callback's pointer
            // stays valid for the writer's whole life, including its flush in Close().
            std::optional<BaiStreamBuilder> baiBuilder{};
            if (config.BaiOutput) {
                baiBuilder.emplace(merged.NumReferences());
            }
            BamWriter writer{[&] {
                if (baiBuilder) {
                    BaiStreamBuilder* const builderPtr{&*baiBuilder};
                    return BamWriter{
                        finalTmp, merged, outWriterConfig,
                        [builderPtr](std::int64_t vo, std::span<const std::byte> bytes) {
                            builderPtr->Observe(VirtualOffset{static_cast<std::uint64_t>(vo)},
                                                bytes);
                        }};
                }
                return BamWriter{finalTmp, merged, outWriterConfig};
            }()};

            // Measure the consumer thread's real CPU across the heap-merge loop only
            // (not Close(), which just joins the flush). Reported as MergeCpuSeconds so
            // the verdict can tell genuine merge-CPU limits from core starvation.
            const std::uint64_t mergeCpuStart{ThreadCpuNanos()};
            numRecords = MergeReadAheadIntoWriter(readAhead, writer, ctx);
            runtime.MergeCpuSeconds = nsToSeconds(ThreadCpuNanos() - mergeCpuStart);
            writer.Close();

            // The last record's end offset is the position just past the final byte
            // written (EndVirtualOffset), available once Close() has flushed everything.
            if (baiBuilder) {
                baiBuilder->Finalize(writer.EndVirtualOffset()).ToFile(*baiTmp);
            }

            // Snapshot bottleneck counters after the pipeline has fully drained.
            runtime.HeapMerge = true;
            const MergeReadAheadStats raStats{readAhead.Stats()};
            const WriterMetrics writerMetrics{writer.GetMetrics()};
            runtime.InputWaitSeconds = nsToSeconds(raStats.ConsumerInputWaitNs);
            runtime.BudgetWaitSeconds = nsToSeconds(raStats.ProducerBudgetWaitNs);
            runtime.InputIoReadSeconds = nsToSeconds(raStats.IoReadNs);
            runtime.InputDecompressSeconds = nsToSeconds(raStats.DecompressNs);
            runtime.PeakInFlightBytes = static_cast<std::int64_t>(raStats.PeakInFlightBytes);
            runtime.BudgetBytes = static_cast<std::int64_t>(raStats.BudgetBytes);
            runtime.OutputWaitSeconds = nsToSeconds(writerMetrics.Bgzf.CallerStallNs);
            runtime.OutputCompressSeconds = nsToSeconds(writerMetrics.Bgzf.CompressNs);
            runtime.OutputWriteSeconds = nsToSeconds(writerMetrics.Bgzf.IoWriteNs);
            runtime.WriterStalls = static_cast<std::int64_t>(writerMetrics.Bgzf.WriterStalls);
        }
    }

    detail::AtomicRename(finalTmp, output, "merge");
    if (baiTmp) {
        detail::AtomicRename(*baiTmp, *config.BaiOutput, "merge BAI");
    }

    runtime.WallSeconds =
        nsToSeconds(static_cast<std::uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                                                   std::chrono::steady_clock::now() - wallStart)
                                                   .count()));
    std::int64_t bytesIn{0};
    for (const std::filesystem::path& input : inputs) {
        std::error_code sizeError{};
        const std::uintmax_t size{std::filesystem::file_size(input, sizeError)};
        if (!sizeError) {
            bytesIn += static_cast<std::int64_t>(size);
        }
    }
    runtime.BytesIn = bytesIn;
    {
        std::error_code sizeError{};
        const std::uintmax_t size{std::filesystem::file_size(output, sizeError)};
        runtime.BytesOut = sizeError ? 0 : static_cast<std::int64_t>(size);
    }

    return MergeStats{
        .NumRecords = numRecords,
        .NumInputs = std::size(inputs),
        .Passthrough = passthrough,
        .Runtime = runtime,
    };
}

}  // namespace Samoa
}  // namespace PacBio
