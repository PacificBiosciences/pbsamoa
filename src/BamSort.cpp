#include <pbsamoa/io/BamMerge.hpp>
#include <pbsamoa/io/BamSort.hpp>

#include "WriterUtils.hpp"

#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/core/Tags.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamWriter.hpp>

#include <pbsamoa/PbSamoaLibraryInfo.hpp>

#include <algorithm>
#include <filesystem>
#include <format>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <thread>
#include <variant>
#include <vector>

#include <cstddef>
#include <cstdint>

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

/// Per-sort comparison context: shared by run-build, spill, and merge.
struct SortContext
{
    SortOrder Order{SortOrder::COORDINATE};
    TagKey Tag{};
    std::int32_t NumReferences{0};
};

/// A record paired with its precomputed comparison keys. Computed once per
/// record (never per comparison): `CoordKey` always; `TagSortValue` only for tag
/// order (engaged == tag present; nullopt == tag missing, which sorts first).
struct Item
{
    RawRecord Record;
    std::uint64_t CoordKey{0};
    std::optional<TagValue> TagSortValue{};
};

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
        .TagSortValue = std::move(tagValue),
    };
}

int CompareU64(std::uint64_t a, std::uint64_t b) { return (a < b) ? -1 : (a > b) ? 1 : 0; }

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
int RecordCompare(SortOrder order, const Item& a, const Item& b)
{
    switch (order) {
        case SortOrder::COORDINATE:
            return CompareU64(a.CoordKey, b.CoordKey);
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
            return CompareU64(a.CoordKey, b.CoordKey);  // fallback on equal/both-missing
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
    std::unique_ptr<BamRawReader> Reader;  // engaged => file-backed
    std::vector<Item> Memory;              // used iff Reader == nullptr
    std::size_t MemoryPos{0};

    std::optional<Item> Next(const SortContext& ctx)
    {
        if (Reader) {
            std::optional<RawRecord> record{Reader->ReadRecord()};
            if (!record) {
                return std::nullopt;
            }
            return MakeItem(std::move(*record), ctx);
        }
        if (MemoryPos < std::size(Memory)) {
            return std::move(Memory[MemoryPos++]);
        }
        return std::nullopt;
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
    return "unknown";
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

SamHeader BuildOutputHeader(const SamHeader& input, const SortConfig& config)
{
    SamHeader header{input};
    header.SetSortOrder(SortOrderText(config.Order));
    header.SetGroupOrder("");  // drop any stale group-order
    AppendProgramRecord(header, "pbsamoa.sort", config.CommandLine);
    return header;
}

/// k-way merge of sorted sources (spilled run files plus an optional already-sorted
/// in-memory run) into \p outPath. Returns the number of records written.
std::int64_t MergeRuns(const std::vector<std::filesystem::path>& runPaths,
                       std::vector<Item> memoryRun, const std::filesystem::path& outPath,
                       const SamHeader& outHeader, const BamWriterConfig& outWriterConfig,
                       const SortContext& ctx)
{
    std::vector<MergeSource> sources{};
    sources.reserve(std::size(runPaths) + 1);
    for (const std::filesystem::path& path : runPaths) {
        MergeSource source{};
        source.Reader = std::make_unique<BamRawReader>(path, BamRawReaderConfig{.BgzfWorkers = 0});
        sources.push_back(std::move(source));
    }
    {
        MergeSource source{};
        source.Memory = std::move(memoryRun);
        sources.push_back(std::move(source));
    }

    struct HeapEntry
    {
        Item Value;
        std::size_t Source{0};
    };

    // A max-heap whose "max" is the smallest record: invert the record order so
    // the smallest record floats to the top; ties break on lower source index
    // (== earlier input run), reconstructing input order for equal keys.
    const auto heapLess = [&ctx](const HeapEntry& a, const HeapEntry& b) {
        const int c{RecordCompare(ctx.Order, a.Value, b.Value)};
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
    BamWriter writer{outPath, outHeader, outWriterConfig};
    while (!heap.empty()) {
        std::ranges::pop_heap(heap, heapLess);
        HeapEntry entry{std::move(heap.back())};
        heap.pop_back();

        writer.Write(entry.Value.Record);
        ++recordsWritten;

        if (std::optional<Item> next{sources[entry.Source].Next(ctx)}) {
            heap.push_back(HeapEntry{.Value = std::move(*next), .Source = entry.Source});
            std::ranges::push_heap(heap, heapLess);
        }
    }
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
void ValidateMergeSortOrder(const SamHeader& header, SortOrder order,
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
    if (header.SortOrder() != required) {
        throw std::runtime_error{
            std::format("MergeBam: input {} has @HD SO:{}, expected SO:{} for this merge order",
                        path.string(), header.SortOrder(), required)};
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

}  // namespace

SortStats SortBam(const std::filesystem::path& input, const std::filesystem::path& output,
                  const SortConfig& config)
{
    if (!std::filesystem::exists(input)) {
        throw std::runtime_error{
            std::format("SortBam: input file does not exist: {}", input.string())};
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
    };
    const auto less = [&ctx](const Item& a, const Item& b) {
        return RecordCompare(ctx.Order, a, b) < 0;
    };

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
            Item item{MakeItem(std::move(*record), ctx)};
            const std::size_t footprint{Footprint(item)};
            if (!run.empty() && ((runBytes + footprint) > maxMemory)) {
                spill();
            }
            runBytes += footprint;
            run.push_back(std::move(item));
            ++numRecords;
        }

        outHeader = BuildOutputHeader(inputHeader, config);
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

    const std::size_t numThreads{ResolveThreads(config.NumThreads)};
    const int finalLevel{std::clamp(config.CompressionLevel, 1, 12)};

    // Build the merged header from the first input, requiring every input to share
    // an identical @SQ list and the expected sort order; union @RG / @PG / @CO.
    SamHeader merged{};
    {
        BamRawReader reader{inputs.front(), BamRawReaderConfig{.BgzfWorkers = 0}};
        merged = reader.Header();
    }
    ValidateMergeSortOrder(merged, config.Order, inputs.front());

    for (std::size_t i{1}; i < std::size(inputs); ++i) {
        BamRawReader reader{inputs[i], BamRawReaderConfig{.BgzfWorkers = 0}};
        const SamHeader& header{reader.Header()};
        if (!SameReferences(merged, header)) {
            throw std::runtime_error{
                std::format("MergeBam: input {} has @SQ references incompatible with {} "
                            "(identical reference lists are required)",
                            inputs[i].string(), inputs.front().string())};
        }
        ValidateMergeSortOrder(header, config.Order, inputs[i]);
        UnionReadGroups(merged, header);
        UnionProgramRecords(merged, header);
        UnionComments(merged, header);
    }

    const SortContext ctx{
        .Order = config.Order,
        .Tag = TagKey{config.Tag[0], config.Tag[1]},
        .NumReferences = merged.NumReferences(),
    };

    merged.SetSortOrder(SortOrderText(config.Order));
    merged.SetGroupOrder("");
    AppendProgramRecord(merged, "pbsamoa.merge", config.CommandLine);

    TempFileCleanup cleanup{};
    const std::filesystem::path finalTmp{detail::TemporaryWritePath(output)};
    cleanup.Add(finalTmp);
    const BamWriterConfig outWriterConfig{
        .BgzfConfig = {.CompressionLevel = finalLevel, .BgzfWorkers = numThreads},
        .UseTempFile = false,
    };

    const std::int64_t numRecords{MergeRuns(inputs, {}, finalTmp, merged, outWriterConfig, ctx)};

    detail::AtomicRename(finalTmp, output, "merge");

    return MergeStats{
        .NumRecords = numRecords,
        .NumInputs = std::size(inputs),
    };
}

}  // namespace Samoa
}  // namespace PacBio
