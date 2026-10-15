#include "Depth.hpp"

#include "../../ParallelUtils.hpp"
#include "../../PathUtils.hpp"
#include "../ParseUtils.hpp"
#include "../StatusStream.hpp"

#include <pbsamoa/PbSamoaLibraryInfo.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/index/BaiIndex.hpp>
#include <pbsamoa/io/BamRawReader.hpp>

#include <pbcopper/cli2/Interface.h>
#include <pbcopper/cli2/Option.h>
#include <pbcopper/cli2/PositionalArgument.h>
#include <pbcopper/cli2/Results.h>
#include <pbcopper/parallel/ThreadPool.h>

#include <algorithm>
#include <atomic>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <iterator>
#include <memory>
#include <optional>
#include <print>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <exception>

namespace PacBio {
namespace Samoa {
namespace Depth {
namespace {

constexpr std::int32_t MAX_THREADS{64};
// The default --excl-flags of samtools depth: UNMAP, SECONDARY, QCFAIL, DUP.
constexpr std::uint16_t EXCLUDED_FLAGS{0x704};
constexpr std::string_view HEADER_LINE{"#contig\tstart\tend\tmean_depth\n"};

const CLI_v2::Option BinSize{
    R"({
    "names" : ["bin-size"],
    "description" : "Size of the bins in bp. The last bin of a contig is shorter.",
    "type" : "integer",
    "default" : 1000
})"};

const CLI_v2::Option Threads{
    R"({
    "names" : ["threads"],
    "description" : "Threads for BGZF reading and for the depth sums; -1 = auto (min of hardware concurrency and 8). Each thread keeps 8 bytes per bin.",
    "type" : "integer",
    "default" : -1
})"};

const CLI_v2::Option StatusJson{
    R"({
    "names" : ["status-json"],
    "description" : "Write a live status stream to this file: one JSON object per line (run info, phases, progress each second, final).",
    "type" : "string",
    "default" : ""
})"};

const CLI_v2::PositionalArgument Input{
    R"({
    "name" : "input",
    "description" : "Aligned input BAM file. It does not have to be sorted.",
    "type" : "file"
})"};

const CLI_v2::PositionalArgument Output{
    R"({
    "name" : "output",
    "description" : "Output TSV file: contig, start (0-based), end, mean depth of each bin.",
    "type" : "file"
})"};

/// The bins of all contigs in one array, in header order. The bins of contig i are
/// [Offsets[i], Offsets[i + 1]).
struct BinLayout
{
    std::int64_t BinSize{0};
    std::vector<std::string> Names;
    std::vector<std::int64_t> Lengths;
    std::vector<std::int64_t> Offsets;

    [[nodiscard]] std::int64_t NumBins() const { return Offsets.back(); }
};

BinLayout MakeLayout(const SamHeader& header, const std::int64_t binSize)
{
    BinLayout layout{.BinSize = binSize, .Offsets = {0}};
    for (const ReferenceSequence& sequence : header.ReferenceSequences()) {
        const std::int64_t length{sequence.Length()};
        layout.Names.emplace_back(sequence.Name());
        layout.Lengths.push_back(length);
        layout.Offsets.push_back(layout.Offsets.back() + ((length + binSize - 1) / binSize));
    }
    return layout;
}

struct DepthCounts
{
    std::int64_t Records{0};
    // records that pass the flag filter and add to the depth
    std::int64_t CountedRecords{0};
    // the sum of the depth over all positions
    std::int64_t DepthBases{0};
};

/// The sums of one worker thread. Each worker owns one, so the workers do not share
/// memory that they write.
struct alignas(detail::PIPELINE_CACHE_LINE_SIZE) WorkerSums
{
    std::vector<std::int64_t> Sums;
    DepthCounts Counts;
};

/// Add 1 to each position of [begin, end) of a contig whose bins start at \p sums.
void AddSpan(const std::span<std::int64_t> sums, const std::int64_t binSize, std::int64_t begin,
             const std::int64_t end)
{
    while (begin < end) {
        const std::int64_t bin{begin / binSize};
        const std::int64_t binEnd{std::min(end, (bin + 1) * binSize)};
        sums[bin] += binEnd - begin;
        begin = binEnd;
    }
}

/// Add the positions that a record covers, as `samtools depth -J`: M, =, X and D count,
/// N does not.
void AddRecord(const RawRecordView& view, const BinLayout& layout, WorkerSums& worker)
{
    ++worker.Counts.Records;
    const std::int32_t refId{view.RefId()};
    if (((view.Flag() & EXCLUDED_FLAGS) != 0) || (refId < 0)) {
        return;
    }
    if (refId >= std::ssize(layout.Lengths)) {
        throw std::runtime_error{
            std::format("[pbsamoa] record {} has reference id {}, but the header has {} references",
                        view.Name(), refId, std::size(layout.Lengths))};
    }
    ++worker.Counts.CountedRecords;
    const std::int64_t length{layout.Lengths[refId]};
    const std::span<std::int64_t> sums{
        std::span<std::int64_t>{worker.Sums}.subspan(layout.Offsets[refId])};
    std::int64_t pos{view.Pos()};
    for (const CigarOp op : view.CigarOps()) {
        const std::int64_t opLength{op.Length()};
        switch (op.Type()) {
            using enum CigarOpType;
            case M:
            case D:
            case EQ:
            case X: {
                // a record that runs past the end of its contig counts up to the end
                const std::int64_t end{std::min(pos + opLength, length)};
                worker.Counts.DepthBases += std::max(end - pos, std::int64_t{0});
                AddSpan(sums, layout.BinSize, pos, end);
                pos += opLength;
                break;
            }
            case N:
                pos += opLength;
                break;
            case I:
            case S:
            case H:
            case P:
                break;
        }
    }
}

/// Read the BAM once. Worker threads walk the CIGARs into their own sums; the sums add up
/// at the end. Adds the records of each batch to \p done for the status stream.
WorkerSums SumDepth(BamRawReader& reader, const BinLayout& layout, const std::size_t numThreads,
                    std::atomic<std::int64_t>& done)
{
    std::vector<WorkerSums> workers(numThreads);
    for (WorkerSums& worker : workers) {
        worker.Sums.resize(layout.NumBins());
    }
    {
        Parallel::ThreadPool<> pool{
            Parallel::ThreadPool<>::Config{.NumThreads = numThreads, .QueueMultiplier = 2}};
        while (std::optional<RawRecordBatch> nextBatch{reader.ReadBatch()}) {
            const auto batch{std::make_shared<RawRecordBatch>(std::move(*nextBatch))};
            const std::size_t recordCount{batch->RecordCount()};
            for (std::size_t first{0}; first < recordCount;
                 first += detail::MAX_PARALLEL_CHUNK_SIZE) {
                const std::size_t last{
                    std::min(recordCount, first + detail::MAX_PARALLEL_CHUNK_SIZE)};
                pool.Submit([&workers, &layout, &done, batch, first,
                             last](const Parallel::ThreadIndex thread) {
                    WorkerSums& worker{workers[thread.Value()]};
                    for (std::size_t i{first}; i < last; ++i) {
                        AddRecord(batch->View(i), layout, worker);
                    }
                    const std::int64_t count = last - first;
                    done.fetch_add(count, std::memory_order::relaxed);
                });
            }
        }
        pool.Finalize();
    }
    WorkerSums total{std::move(workers.front())};
    for (const WorkerSums& worker : std::span{workers}.subspan(1)) {
        std::ranges::transform(total.Sums, worker.Sums, std::begin(total.Sums), std::plus{});
        total.Counts.Records += worker.Counts.Records;
        total.Counts.CountedRecords += worker.Counts.CountedRecords;
        total.Counts.DepthBases += worker.Counts.DepthBases;
    }
    return total;
}

/// One line for each bin. The mean prints with the shortest text that reads back as the
/// same double, so mean * (end - start) gives the integer sum again.
void WriteBins(const std::filesystem::path& path, const BinLayout& layout,
               const std::span<const std::int64_t> sums)
{
    std::ofstream out{path, std::ios::binary};
    out.write(std::data(HEADER_LINE), std::ssize(HEADER_LINE));
    std::string buffer;
    for (std::size_t contig{0}; contig < std::size(layout.Names); ++contig) {
        buffer.clear();
        const std::string_view name{layout.Names[contig]};
        const std::int64_t length{layout.Lengths[contig]};
        const std::int64_t offset{layout.Offsets[contig]};
        for (std::int64_t bin{0}; bin < (layout.Offsets[contig + 1] - offset); ++bin) {
            const std::int64_t start{bin * layout.BinSize};
            const std::int64_t end{std::min(start + layout.BinSize, length)};
            const double sum = sums[offset + bin];
            const double span = end - start;
            std::format_to(std::back_inserter(buffer), "{}\t{}\t{}\t{}\n", name, start, end,
                           sum / span);
        }
        out.write(std::data(buffer), std::ssize(buffer));
    }
    out.close();
    if (!out) {
        throw std::runtime_error{std::format("[pbsamoa] could not write {}", path.string())};
    }
}

/// The counts as the metrics of the status stream.
JSON::Json CountsJson(const DepthCounts& counts, const BinLayout& layout)
{
    return {{"records", counts.Records},
            {"counted_records", counts.CountedRecords},
            {"depth_bases", counts.DepthBases},
            {"bins", layout.NumBins()},
            {"bin_size", layout.BinSize}};
}

/// The records of the input from its BAI (the mapped and unmapped counts of all references and
/// the unplaced reads), as the total of the progress. None without a BAI, with a BAI older than
/// the BAM, or with a BAI without counts. The total only feeds the status stream.
std::optional<std::int64_t> RecordTotal(const std::filesystem::path& input)
{
    const std::filesystem::path bai{SidecarPath(input, ".bai")};
    std::error_code ec;
    if (!std::filesystem::exists(bai, ec) ||
        (std::filesystem::last_write_time(bai, ec) < std::filesystem::last_write_time(input, ec)) ||
        ec) {
        return std::nullopt;
    }
    try {
        const BaiIndex index{BaiIndex::FromFile(bai)};
        const std::int64_t total = index.MappedCount() + index.UnmappedCount();
        return total > 0 ? std::optional{total} : std::nullopt;
    } catch (const std::exception&) {
        return std::nullopt;
    }
}

struct Settings
{
    std::filesystem::path Input;
    std::filesystem::path Output;
    std::int64_t BinSize{0};
    std::size_t NumThreads{0};
    std::filesystem::path StatusJson;
};

Settings ParseSettings(const CLI_v2::Results& results)
{
    const std::vector<std::string>& positional{results.PositionalArguments()};
    if (std::size(positional) != 2) {
        throw std::runtime_error{"depth requires exactly two arguments: <input.bam> <output.tsv>"};
    }
    Settings settings{.Input = positional[0], .Output = positional[1]};
    const std::int64_t binSize{results[BinSize]};
    if (binSize < 1) {
        throw std::runtime_error{std::format("--bin-size must be at least 1: {}", binSize)};
    }
    settings.BinSize = binSize;
    const std::int32_t requestedThreads{results[Threads]};
    settings.NumThreads =
        std::ranges::max(Tools::ResolveNumWorkers(requestedThreads, MAX_THREADS), std::size_t{1});
    const std::string statusJson = results[StatusJson];
    settings.StatusJson = statusJson;
    return settings;
}

}  // namespace

CLI_v2::Interface CreateInterface()
{
    CLI_v2::Interface interface{
        "pbsamoa depth",
        "Mean read depth in fixed bins of each contig, as samtools depth -J: deletions count, "
        "unmapped, secondary, QC-fail and duplicate records do not.",
        LibraryFormattedVersion()};
    // `depth` keeps its own --threads (-1 = auto, capped at min(hw, 8)), as `subsample`.
    interface.DisableNumThreadsOption();
    interface.AddOptions({BinSize, Threads, StatusJson});
    interface.AddPositionalArguments({Input, Output});
    return interface;
}

int Runner(const CLI_v2::Results& results)
{
    const Settings settings{ParseSettings(results)};
    const std::string commandLine{"pbsamoa " + results.InputCommandLine()};
    // outside the try: the catch writes the final line of a failure
    Tools::StatusStream status;
    status.Open(settings.StatusJson,
                Tools::MakeRunInfo("pbsamoa", "depth", LibraryFormattedVersion(), commandLine,
                                   settings.NumThreads, {"depth"}));
    try {
        BamRawReader reader{settings.Input, BamRawReaderConfig{.BgzfWorkers = settings.NumThreads}};
        const BinLayout layout{MakeLayout(reader.Header(), settings.BinSize)};
        std::atomic<std::int64_t> done{0};
        const std::optional<std::int64_t> indexRecords{RecordTotal(settings.Input)};
        Tools::StatusStream::Phase phase{&status, "depth", &done, "reads", indexRecords};
        const WorkerSums total{SumDepth(reader, layout, settings.NumThreads, done)};
        try {
            WriteBins(settings.Output, layout, total.Sums);
        } catch (...) {
            std::error_code ignored;
            std::filesystem::remove(settings.Output, ignored);
            throw;
        }
        // `=`, not braces: a braced JSON::Json makes a one-element array
        JSON::Json metrics = CountsJson(total.Counts, layout);
        if (indexRecords) {
            metrics["index_records"] = *indexRecords;
        }
        phase.End(metrics);
        std::println(stderr,
                     "pbsamoa depth: {} of {} records counted, {} depth bases in {} bins of {} bp",
                     total.Counts.CountedRecords, total.Counts.Records, total.Counts.DepthBases,
                     layout.NumBins(), layout.BinSize);
        status.Final(
            "FINISHED", {}, std::move(metrics),
            JSON::Json::array({{{"path", std::filesystem::absolute(settings.Output).string()},
                                {"kind", "tsv"}}}));
    } catch (const std::exception& e) {
        status.Final("FAILED", e.what());
        throw;
    }
    return EXIT_SUCCESS;
}

}  // namespace Depth
}  // namespace Samoa
}  // namespace PacBio
