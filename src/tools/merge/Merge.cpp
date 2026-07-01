#include "Merge.hpp"

#include "../../PathUtils.hpp"
#include "../CliUtils.hpp"
#include "../MetricUtils.hpp"

#include <pbsamoa/io/BamMerge.hpp>

#include <algorithm>
#include <filesystem>
#include <format>
#include <print>
#include <stdexcept>
#include <string>
#include <string_view>
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

void PrintUsage()
{
    std::println(stderr,
                 "Usage: pbsamoa merge OUT.bam IN1.bam [IN2.bam ...] [options]\n"
                 "\n"
                 "Merge BAM files (identical @SQ required) into one. With no --order/--concat\n"
                 "the mode is auto-detected from the inputs' @HD SO: all-coordinate or\n"
                 "all-queryname does that sorted merge, all-unsorted/unknown concatenates,\n"
                 "a mix errors.\n"
                 "\n"
                 "Options:\n"
                 "  --order ORDER     coordinate | queryname | tag   (default: auto-detect)\n"
                 "  --concat          byte-concatenate inputs (BGZF passthrough), SO:unsorted\n"
                 "  --tag XX          2-char tag inputs are sorted by (required iff --order tag)\n"
                 "  --threads N       CPU worker pools (decode + compress), 0/auto = min(hw,8)\n"
                 "  --decode-threads N   input-decompression pool, 0 = inherit --threads\n"
                 "  --compress-threads N output-compression pool,  0 = inherit --threads\n"
                 "  --memory SIZE     input read-ahead budget, K/M/G suffix (default 768M)\n"
                 "  --batch-bytes SIZE   per-input handoff batch target, K/M/G (default 256K)\n"
                 "  --writer-queue N  output writer queue depth in blocks (default 256)\n"
                 "  --compression L   output BGZF level [1,12]         (default 6)\n"
                 "  --bai             write OUT.bam.bai during the merge (coordinate merges\n"
                 "                    that take the heap path only; not --concat, queryname,\n"
                 "                    tag, or the disjoint-chain passthrough)\n"
                 "\n"
                 "Thread/memory knobs (except --compression) apply to the sorted-merge path\n"
                 "only; --concat is single-threaded passthrough.");
}

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

int Runner(int argc, char** argv)
{
    MergeConfig config{};
    std::vector<std::string_view> positional{};
    bool tagProvided{false};
    bool baiRequested{false};

    for (int i{0}; i < argc; ++i) {
        const std::string_view arg{argv[i]};
        if (arg == "--order") {
            config.Order = Tools::ParseSortOrder(argc, argv, i);
        } else if (arg == "--concat") {
            config.Concat = true;
        } else if (arg == "--tag") {
            config.Tag = Tools::ParseSortTag(argc, argv, i);
            tagProvided = true;
        } else if (arg == "--threads") {
            config.NumThreads = Tools::ParseThreadsOption(argc, argv, i);
        } else if (arg == "--decode-threads") {
            config.DecodeThreads = Tools::ParseThreadsOption(argc, argv, i, "--decode-threads");
        } else if (arg == "--compress-threads") {
            config.CompressThreads = Tools::ParseThreadsOption(argc, argv, i, "--compress-threads");
        } else if (arg == "--memory") {
            config.ReadAheadMemory = Tools::ParseMemoryOption(argc, argv, i);
        } else if (arg == "--batch-bytes") {
            config.BatchBytes = Tools::ParseMemoryOption(argc, argv, i, "--batch-bytes").Value();
            if (config.BatchBytes == 0) {
                throw std::runtime_error{"--batch-bytes must be >= 1"};
            }
        } else if (arg == "--writer-queue") {
            config.WriterQueueCapacity =
                Tools::ParsePositiveCountOption(argc, argv, i, "--writer-queue");
        } else if (arg == "--compression") {
            config.CompressionLevel = Tools::ParseCompressionLevelOption(argc, argv, i);
        } else if (arg == "--bai") {
            baiRequested = true;
        } else if (arg.starts_with("--")) {
            throw std::runtime_error{std::format("unknown option: {}", arg)};
        } else {
            positional.push_back(arg);
        }
    }

    // First positional is the output; the rest are inputs.
    if (std::size(positional) < 2) {
        PrintUsage();
        return EXIT_FAILURE;
    }
    if ((config.Order == SortOrder::TAG) && !tagProvided) {
        throw std::runtime_error{"--order tag requires --tag XX"};
    }

    const std::filesystem::path output{positional.front()};
    if (baiRequested) {
        config.BaiOutput = SidecarPath(output, ".bai");
    }
    std::vector<std::filesystem::path> inputs{};
    inputs.reserve(std::size(positional) - 1);
    for (std::size_t i{1}; i < std::size(positional); ++i) {
        inputs.emplace_back(positional[i]);
    }

    config.CommandLine = Tools::BuildCommandLine("merge", argc, argv);

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
