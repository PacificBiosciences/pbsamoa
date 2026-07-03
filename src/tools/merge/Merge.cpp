#include "Merge.hpp"

#include "../../PathUtils.hpp"
#include "../CliUtils.hpp"
#include "../MetricUtils.hpp"

#include <pbsamoa/PbSamoaLibraryInfo.hpp>
#include <pbsamoa/io/BamMerge.hpp>

#include <pbcopper/cli2/Interface.h>
#include <pbcopper/cli2/Option.h>
#include <pbcopper/cli2/PositionalArgument.h>
#include <pbcopper/cli2/Results.h>

#include <algorithm>
#include <filesystem>
#include <format>
#include <print>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <cstddef>
#include <cstdint>
#include <cstdlib>

#include <sys/resource.h>

namespace PacBio {
namespace Samoa {
namespace MergeTool {

namespace {

// CLIv2 option constants ---------------------------------------------------

const CLI_v2::Option Order{
    R"({
    "names" : ["order"],
    "description" : "Sort order: coordinate | queryname | tag (default: auto-detect from inputs).",
    "type" : "string"
})"};
// Note: choices omitted — an empty default is not valid inside pbcopper's
// choices set. Runner-side ParseSortOrder validates the value instead.

const CLI_v2::Option Concat{
    R"({
    "names" : ["concat"],
    "description" : "Byte-concatenate inputs (BGZF passthrough), output SO:unsorted."
})"};

const CLI_v2::Option Tag{
    R"({
    "names" : ["tag"],
    "description" : "Two-character tag inputs are sorted by (required with --order tag).",
    "type" : "string"
})"};

const CLI_v2::Option Threads{
    R"({
    "names" : ["threads"],
    "description" : "CPU worker pools (decode + compress); 0/auto = min(hw, 8).",
    "type" : "unsigned integer",
    "default" : 0
})"};

const CLI_v2::Option DecodeThreads{
    R"({
    "names" : ["decode-threads"],
    "description" : "Input-decompression pool; 0 = inherit --threads.",
    "type" : "unsigned integer",
    "default" : 0
})"};

const CLI_v2::Option CompressThreads{
    R"({
    "names" : ["compress-threads"],
    "description" : "Output-compression pool; 0 = inherit --threads.",
    "type" : "unsigned integer",
    "default" : 0
})"};

const CLI_v2::Option Memory{
    R"({
    "names" : ["memory"],
    "description" : "Input read-ahead budget; K/M/G suffix accepted.",
    "type" : "string",
    "default" : "768M"
})"};

const CLI_v2::Option BatchBytes{
    R"({
    "names" : ["batch-bytes"],
    "description" : "Per-input handoff batch target; K/M/G suffix accepted (default 256K).",
    "type" : "string",
    "default" : "256K"
})"};

const CLI_v2::Option WriterQueue{
    R"({
    "names" : ["writer-queue"],
    "description" : "Output writer queue depth in blocks (>= 1, default 256).",
    "type" : "integer",
    "default" : 256
})"};

const CLI_v2::Option Compression{
    R"({
    "names" : ["compression"],
    "description" : "Output BGZF compression level in [1, 12].",
    "type" : "integer",
    "default" : 6
})"};

const CLI_v2::Option Bai{
    R"({
    "names" : ["bai"],
    "description" : "Write OUT.bam.bai during the merge (coordinate heap-merge path only; not --concat, queryname, tag, or disjoint-chain passthrough)."
})"};

const CLI_v2::PositionalArgument Output{
    R"({
    "name" : "output",
    "description" : "Output BAM file.",
    "type" : "file"
})"};

const CLI_v2::PositionalArgument Input{
    R"({
    "name" : "input",
    "description" : "Input BAM file(s).",
    "type" : "file"
})"};

// Runtime profiling helpers ------------------------------------------------

double TimevalSeconds(const timeval& tv)
{
    return static_cast<double>(tv.tv_sec) + (static_cast<double>(tv.tv_usec) * 1e-6);
}

/// Process CPU time (user + system, all threads) consumed between two samples.
double ProcessCpuSeconds(const rusage& before, const rusage& after)
{
    return (TimevalSeconds(after.ru_utime) - TimevalSeconds(before.ru_utime)) +
           (TimevalSeconds(after.ru_stime) - TimevalSeconds(before.ru_stime));
}

/// Process peak resident set, normalized to bytes (ru_maxrss is bytes on macOS,
/// kibibytes on Linux). This is a process-wide high-water mark, not a delta.
std::int64_t PeakRssBytes(const rusage& usage)
{
#if defined(__APPLE__)
    return static_cast<std::int64_t>(usage.ru_maxrss);
#else
    return static_cast<std::int64_t>(usage.ru_maxrss) * 1024;
#endif
}

/// Print the runtime profile and a bottleneck verdict to stderr. The k-way merge
/// funnels every record through one consumer thread, so its wall time splits into
/// input wait (Next), output wait (writer queue), measured merge CPU (a per-thread
/// clock), and a stall/other remainder. That remainder — wall minus the two waits
/// minus measured CPU — is the key signal: a large remainder with measured CPU low
/// means the consumer is starved of cores, not compute-bound. The verdict weighs
/// these against the compress/decode wall-equivalents (pool CPU / pool threads) and
/// thread oversubscription, so a compression-CPU floor is named as such rather than
/// mislabeled "serial merge".
void PrintRuntimeReport(const MergeStats& stats, double cpuSeconds, std::int64_t peakRssBytes)
{
    const MergeRuntimeStats& rt{stats.Runtime};
    const double wall{rt.WallSeconds};
    const double cpuPerWall{(wall > 0.0) ? (cpuSeconds / wall) : 0.0};
    const double inMiB{Tools::ToMiB(rt.BytesIn)};
    const double outMiB{Tools::ToMiB(rt.BytesOut)};

    std::println(stderr, "  runtime: wall {:.2f}s, peak RSS {:.0f} MB, cpu {:.1f}s ({:.1f}x wall)",
                 wall, Tools::ToMiB(peakRssBytes), cpuSeconds, cpuPerWall);
    if (stats.NumRecords >= 0) {
        std::println(stderr, "    records {} ({:.0f}/s)", stats.NumRecords,
                     Tools::RateOrZero(stats.NumRecords, wall));
    }
    std::println(stderr, "    input  {:.0f} MB ({:.0f} MB/s)   output {:.0f} MB ({:.0f} MB/s)",
                 inMiB, Tools::RateOrZero(inMiB, wall), outMiB, Tools::RateOrZero(outMiB, wall));

    if (!rt.HeapMerge) {
        std::println(stderr, "  pipeline: BGZF block passthrough (no decode/compress pipeline)");
        std::println(stderr, "  verdict: I/O-COPY bound (limited by raw read + write bandwidth)");
        return;
    }

    const double fIn{(wall > 0.0) ? (rt.InputWaitSeconds / wall) : 0.0};
    const double fOut{(wall > 0.0) ? (rt.OutputWaitSeconds / wall) : 0.0};
    const double mergeCpu{rt.MergeCpuSeconds};
    const double fMergeCpu{(wall > 0.0) ? (mergeCpu / wall) : 0.0};

    // Time the single consumer thread was neither computing nor in a tracked wait:
    // CPU-starvation (the compress pool oversubscribing cores) plus the merge's
    // setup/flush outside the timed loop. Clamped because the output-queue spin is a
    // busy yield, so its wall (output wait) overlaps a sliver of measured merge CPU.
    const double other{std::max(0.0, wall - rt.InputWaitSeconds - rt.OutputWaitSeconds - mergeCpu)};
    const double fOther{(wall > 0.0) ? (other / wall) : 0.0};

    // Pool wall-equivalents: aggregate pool CPU spread across the pool's own threads,
    // i.e. the minimum wall each stage needs if its threads each held a core. When a
    // stage's wall-equivalent approaches total wall it is the throughput floor,
    // independent of queue depth — nothing downstream can finish sooner.
    const double compressWallEq{(rt.CompressThreads > 0) ? (rt.OutputCompressSeconds /
                                                            static_cast<double>(rt.CompressThreads))
                                                         : rt.OutputCompressSeconds};
    const double decodeWallEq{
        (rt.DecodeThreads > 0) ? (rt.InputDecompressSeconds / static_cast<double>(rt.DecodeThreads))
                               : rt.InputDecompressSeconds};
    const double fCompress{(wall > 0.0) ? (compressWallEq / wall) : 0.0};

    std::println(stderr, "  pipeline [heap k-way merge; consumer thread = serializer]");
    std::println(stderr,
                 "    input wait  {:.2f}s ({:.0f}%)   io-read {:.0f}ms  decompress {:.0f}ms "
                 "({:.1f}s wall-equiv / {} thr)",
                 rt.InputWaitSeconds, fIn * 100.0, rt.InputIoReadSeconds * 1e3,
                 rt.InputDecompressSeconds * 1e3, decodeWallEq, rt.DecodeThreads);
    std::println(stderr,
                 "    output wait {:.2f}s ({:.0f}%)   compress {:.0f}ms ({:.1f}s wall-equiv / {} "
                 "thr)  write {:.0f}ms  writer-stalls {}",
                 rt.OutputWaitSeconds, fOut * 100.0, rt.OutputCompressSeconds * 1e3, compressWallEq,
                 rt.CompressThreads, rt.OutputWriteSeconds * 1e3, rt.WriterStalls);
    std::println(stderr, "    merge cpu   {:.2f}s ({:.0f}%)   heap + key compare (measured)",
                 mergeCpu, fMergeCpu * 100.0);
    std::println(stderr, "    stall/other {:.2f}s ({:.0f}%)   CPU contention + setup/flush", other,
                 fOther * 100.0);
    std::println(stderr, "    read-ahead  peak {:.0f}/{:.0f} MB budget, budget wait {:.2f}s",
                 Tools::ToMiB(rt.PeakInFlightBytes), Tools::ToMiB(rt.BudgetBytes),
                 rt.BudgetWaitSeconds);

    if (wall < 0.05) {
        std::println(stderr, "  verdict: run too short to profile reliably");
        return;
    }

    // A run is memory-bound when the read-ahead budget is pinned at its cap and
    // producers spend real time blocked on it: more budget would widen read-ahead.
    const bool budgetPinned{(rt.BudgetBytes > 0) &&
                            (rt.PeakInFlightBytes >= ((rt.BudgetBytes * 9) / 10)) &&
                            ((rt.BudgetWaitSeconds / wall) >= 0.15)};

    // The consumer, packer, and IO writer each need a core on top of the two pools,
    // so a decode+compress thread sum near or above hardware concurrency starves the
    // serial consumer — surfacing as stall/other, not merge CPU. hardware_concurrency
    // counts logical cores; 0 means unknown (skip the check).
    const unsigned hwCores{std::thread::hardware_concurrency()};
    const std::size_t poolThreads{rt.DecodeThreads + rt.CompressThreads};
    const bool oversubscribed{(hwCores > 0) &&
                              ((poolThreads + 3) > static_cast<std::size_t>(hwCores))};

    std::string verdict{};
    std::string hint{};
    if (fCompress >= 0.7) {
        // Compress alone needs ~all the wall on its threads: the throughput floor,
        // even when queues flow and the consumer never blocks on output.
        verdict = "OUTPUT-bound (compression CPU)";
        if (oversubscribed) {
            hint = std::format(
                "compress needs ~{:.0f}s on {} threads (~wall); the {} pool threads exceed {} "
                "cores, so more --compress-threads won't help — lower --compression or cut "
                "--decode-threads",
                compressWallEq, rt.CompressThreads, poolThreads, hwCores);
        } else {
            hint = std::format(
                "compress needs ~{:.0f}s on {} threads (~wall); lower --compression or raise "
                "--compress-threads (now {}) if cores are free",
                compressWallEq, rt.CompressThreads, rt.CompressThreads);
        }
    } else if (budgetPinned && (fIn >= fOut)) {
        verdict = "MEMORY-bound (read-ahead budget capped)";
        hint =
            std::format("raise --memory (budget {:.0f} MB pinned)", Tools::ToMiB(rt.BudgetBytes));
    } else if ((fIn >= 0.34) && (fIn >= fOut) && (fIn >= fMergeCpu) && (fIn >= fOther)) {
        if (rt.InputDecompressSeconds > rt.InputIoReadSeconds) {
            verdict = "INPUT-bound (decode CPU)";
            hint = std::format("raise --decode-threads (now {})", rt.DecodeThreads);
        } else {
            verdict = "INPUT-bound (disk read)";
            hint = "faster input storage, or stage inputs on local disk";
        }
    } else if ((fOut >= 0.34) && (fOut >= fMergeCpu) && (fOut >= fOther)) {
        if (rt.OutputCompressSeconds >= rt.OutputWriteSeconds) {
            verdict = "OUTPUT-bound (compression CPU)";
            hint = std::format("lower --compression or raise --compress-threads (now {})",
                               rt.CompressThreads);
        } else {
            verdict = "OUTPUT-bound (disk write)";
            hint = "faster output storage, or lower --compression";
        }
    } else if (oversubscribed && (fOther >= 0.34) && (fOther >= fMergeCpu)) {
        verdict = "CONTENTION-bound (CPU oversubscription)";
        hint = std::format(
            "the serial consumer is starved: {} pool threads (+3) vs {} cores; lower "
            "--decode-threads and/or --compress-threads",
            poolThreads, hwCores);
    } else if ((fMergeCpu >= 0.34) && (fMergeCpu >= fOther)) {
        verdict = "CPU-bound (serial k-way merge)";
        hint = "consumer compute dominates; fewer inputs or a cheaper sort key";
    } else {
        verdict = "BALANCED (no single limiter)";
        hint = "stages overlap well; scale them together for more throughput";
    }

    std::println(stderr, "  verdict: {} -> {}", verdict, hint);
}

}  // namespace

CLI_v2::Interface CreateInterface()
{
    CLI_v2::Interface interface{
        "pbsamoa merge",
        "Merge BAM files into one (sorted merge or BGZF passthrough concatenation).",
        LibraryFormattedVersion()};
    // `merge` keeps its own thread flags with 0=auto semantics; the built-in
    // --num-threads resolves 0 to the raw hardware count, which would drop each
    // flag's min(hw,8) cap.
    interface.DisableNumThreadsOption();
    interface.AddOptions({Order, Concat, Tag, Threads, DecodeThreads, CompressThreads, Memory,
                          BatchBytes, WriterQueue, Compression, Bai});
    interface.AddPositionalArguments({Output, Input});
    return interface;
}

int Runner(const CLI_v2::Results& results)
{
    MergeConfig config{};

    // --order (optional; empty string = auto-detect from inputs' @HD SO:)
    const std::string order = results[Order];
    if (!order.empty()) {
        config.Order = Tools::ParseSortOrder(order);
    }

    // --concat
    const bool concat = results[Concat];
    config.Concat = concat;

    // --tag (optional; required when --order tag is set)
    const std::string tag = results[Tag];
    if (!tag.empty()) {
        config.Tag = Tools::ParseSortTag(tag);
    }
    if ((config.Order == SortOrder::TAG) && tag.empty()) {
        throw std::runtime_error{"--order tag requires --tag XX"};
    }

    // Thread knobs: read as uint32_t (no exact size_t converter in Results), assign
    // to size_t fields.
    const std::uint32_t threads{results[Threads]};
    config.NumThreads = threads;

    const std::uint32_t decodeThreads{results[DecodeThreads]};
    config.DecodeThreads = decodeThreads;

    const std::uint32_t compressThreads{results[CompressThreads]};
    config.CompressThreads = compressThreads;

    // --memory
    const std::string memory = results[Memory];
    config.ReadAheadMemory = Tools::ParseMemory(memory);

    // --batch-bytes (default "256K" = library's internal default; validate >= 1)
    const std::string batchBytes = results[BatchBytes];
    config.BatchBytes = Tools::ParseMemory(batchBytes, "--batch-bytes").Value();
    if (config.BatchBytes == 0) {
        throw std::runtime_error{"--batch-bytes must be >= 1"};
    }

    // --writer-queue (default 256; must be >= 1)
    const std::int32_t writerQueue{results[WriterQueue]};
    if (writerQueue < 1) {
        throw std::runtime_error{"--writer-queue must be >= 1"};
    }
    config.WriterQueueCapacity = static_cast<std::size_t>(writerQueue);

    // --compression
    const std::int32_t compression{results[Compression]};
    config.CompressionLevel = Tools::CheckCompressionLevel(compression);

    // CLIv2 does not enforce the required positional count; guard before indexing.
    // merge needs at least: <output> <input> (one output + one or more inputs).
    const std::vector<std::string>& positional{results.PositionalArguments()};
    if (positional.size() < 2) {
        throw std::runtime_error{
            "merge requires at least two arguments: <output> <input> [<input> ...]"};
    }

    const std::filesystem::path output{positional[0]};

    // --bai
    const bool bai = results[Bai];
    if (bai) {
        config.BaiOutput = SidecarPath(output, ".bai");
    }

    std::vector<std::filesystem::path> inputs{};
    inputs.reserve(positional.size() - 1);
    for (std::size_t i{1}; i < positional.size(); ++i) {
        inputs.emplace_back(positional[i]);
    }

    config.CommandLine = std::string{"pbsamoa "} + results.InputCommandLine();

    rusage usageBefore{};
    getrusage(RUSAGE_SELF, &usageBefore);

    const MergeStats stats{MergeBam(inputs, output, config)};

    rusage usageAfter{};
    getrusage(RUSAGE_SELF, &usageAfter);

    if (stats.NumRecords < 0) {
        std::println(stderr, "pbsamoa merge: concatenated {} input(s) (records not counted)",
                     stats.NumInputs);
    } else {
        std::println(stderr, "pbsamoa merge: {} records from {} input(s){}", stats.NumRecords,
                     stats.NumInputs, stats.Passthrough ? " (passthrough)" : "");
    }

    PrintRuntimeReport(stats, ProcessCpuSeconds(usageBefore, usageAfter), PeakRssBytes(usageAfter));
    return EXIT_SUCCESS;
}

}  // namespace MergeTool
}  // namespace Samoa
}  // namespace PacBio
