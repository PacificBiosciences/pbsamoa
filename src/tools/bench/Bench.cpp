#include "Bench.hpp"

#include "../../PathUtils.hpp"
#include "../MetricUtils.hpp"
#include "../ParseUtils.hpp"

#include <pbsamoa/PbSamoaLibraryInfo.hpp>
#include <pbsamoa/core/Metrics.hpp>
#include <pbsamoa/index/BaiIndex.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamRecordReader.hpp>
#include <pbsamoa/io/BamWriter.hpp>

#include <pbcopper/cli2/Interface.h>
#include <pbcopper/cli2/Option.h>
#include <pbcopper/cli2/PositionalArgument.h>
#include <pbcopper/cli2/Results.h>

#include <filesystem>
#include <print>
#include <stdexcept>
#include <string>
#include <vector>

#include <cstdint>
#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace Bench {
namespace {

// --bgzf-threads is long-only: `-j` belongs to the built-in --num-threads (disabled
// below via DisableNumThreadsOption()), whose auto-resolution differs from bench's
// -1 = auto (min of hardware concurrency and 10).
const CLI_v2::Option BgzfThreads{
    R"({
    "names" : ["bgzf-threads"],
    "description" : "BGZF decompression worker threads; -1 = auto (min of hardware concurrency and 10).",
    "type" : "integer",
    "default" : -1
})"};

const CLI_v2::Option DecodeThreads{
    R"({
    "names" : ["decode-threads"],
    "description" : "BAM record decode worker threads; 0 = serial.",
    "type" : "integer",
    "default" : 4
})"};

const CLI_v2::PositionalArgument Input{
    R"({
    "name" : "input",
    "description" : "Input BAM file.",
    "type" : "file"
})"};

struct ReadCounts
{
    std::size_t records;
    std::size_t bytes;
};

double ElapsedSecs(std::chrono::steady_clock::time_point start)
{
    const auto elapsed{std::chrono::steady_clock::now() - start};
    return std::chrono::duration<double>(elapsed).count();
}

ReadCounts CountRecordsAndBytes(BamRawReader& reader)
{
    ReadCounts counts{};
    for (const auto& view : reader.Records()) {
        ++counts.records;
        counts.bytes += std::size(view.RawData());
    }
    return counts;
}

void BenchSequentialRead(const std::filesystem::path& path)
{
    const auto start{std::chrono::steady_clock::now()};
    BamRawReader reader{path};
    const ReadCounts counts{CountRecordsAndBytes(reader)};

    const double secs{ElapsedSecs(start)};
    const double inputMiB{Tools::ToMiB(counts.bytes)};
    const double mbPerSec{inputMiB / secs};

    std::println(
        "sequential_read: {} records, {:.1f} MB, {:.3f} sec, {:.1f} "
        "MB/s, {:.0f} records/s",
        counts.records, inputMiB, secs, mbPerSec, Tools::ToRate(counts.records, secs));
}

void BenchBatchRead(const std::filesystem::path& path, ByteLimit limit)
{
    const auto start{std::chrono::steady_clock::now()};
    BamRawReader reader{path};
    std::size_t records{0};
    std::size_t batches{0};

    while (const auto batch = reader.ReadBatch(limit)) {
        records += batch->RecordCount();
        ++batches;
    }

    const double secs{ElapsedSecs(start)};

    std::println(
        "batch_read({} KiB): {} records, {} batches, {:.3f} sec, {:.0f} "
        "records/s",
        limit.Value() / 1024, records, batches, secs, Tools::ToRate(records, secs));
}

void BenchWrite(const std::filesystem::path& srcPath)
{
    const auto tmpPath{std::filesystem::temp_directory_path() / "pbsamoa_bench_write.bam"};

    BamRawReader srcReader{srcPath};

    const auto batch{srcReader.ReadBatch(ByteLimit{256U * 1024U * 1024U})};
    if (!batch) {
        std::println("write: no records to benchmark");
        return;
    }
    const std::size_t numRecords{batch->RecordCount()};

    const auto start{std::chrono::steady_clock::now()};
    {
        BamWriter writer{tmpPath, srcReader.Header()};
        for (std::size_t i{0}; i < numRecords; ++i) {
            writer.Write(batch->RecordData(i));
        }
    }
    const double secs{ElapsedSecs(start)};

    const auto fileSize{std::filesystem::file_size(tmpPath)};

    const double outputMiB{Tools::ToMiB(fileSize)};

    std::println(
        "write: {} records, {:.1f} MB output, {:.3f} sec, {:.1f} MB/s, "
        "{:.0f} records/s",
        numRecords, outputMiB, secs, outputMiB / secs, Tools::ToRate(numRecords, secs));

    std::filesystem::remove(tmpPath);
}

void PrintBgzfMetrics(const BgzfMetrics& m)
{
    const double mbRead{Tools::ToMiB(m.BytesRead)};
    const double mbDecomp{Tools::ToMiB(m.BytesDecompressed)};
    const double ioMs{Tools::ToMs(m.IoReadNs)};
    const double decompMs{Tools::ToMs(m.DecompressNs)};
    const double parseMs{Tools::ToMs(m.RecordParseNs)};

    std::println("    bgzf: {:.1f} MB compressed, {:.1f} MB decompressed, {} blocks", mbRead,
                 mbDecomp, m.BlocksRead);
    std::println("    timing: io {:.1f}ms, decomp {:.1f}ms, parse {:.1f}ms", ioMs, decompMs,
                 parseMs);
    std::println("    pool peak: queue {}, workers {}, results {}", m.Pool.PeakQueueDepth,
                 m.Pool.PeakActiveWorkers, m.Pool.PeakResultQueueDepth);
    std::println("    stalls: io={}, consumer={}, reader={}", m.IoStalls, m.ConsumerStalls,
                 m.ReaderStalls);
}

void PrintDecodeMetrics(const DecodeMetrics& m)
{
    const double decodeMs{Tools::ToMs(m.DecodeNs)};
    const double batchReadMs{Tools::ToMs(m.BatchReadNs)};

    std::println("    decode: {} batches, {} records", m.BatchesDecoded, m.RecordsDecoded);
    std::println("    timing: decode {:.1f}ms, batch_read {:.1f}ms", decodeMs, batchReadMs);
    std::println("    pool peak: queue {}, workers {}", m.Pool.PeakQueueDepth,
                 m.Pool.PeakActiveWorkers);
    std::println("    stalls: producer={}, consumer={}", m.ProducerStalls, m.ConsumerStalls);
}

void BenchPipelineRead(const std::filesystem::path& path, std::size_t bgzfThreads)
{
    const auto start{std::chrono::steady_clock::now()};
    BamRawReader reader{path, BamRawReaderConfig{.BgzfWorkers = bgzfThreads}};
    const ReadCounts counts{CountRecordsAndBytes(reader)};

    const double secs{ElapsedSecs(start)};
    const double inputMiB{Tools::ToMiB(counts.bytes)};
    const double mbPerSec{inputMiB / secs};

    std::println(
        "pipeline_raw_read(bgzf={}): {} records, {:.1f} MB, {:.3f} sec, "
        "{:.1f} MB/s, {:.0f} "
        "records/s",
        bgzfThreads, counts.records, inputMiB, secs, mbPerSec, Tools::ToRate(counts.records, secs));

    PrintBgzfMetrics(reader.GetMetrics());
}

void BenchRecordReader(const std::filesystem::path& path, std::size_t bgzfThreads,
                       std::size_t decodeThreads)
{
    const auto start{std::chrono::steady_clock::now()};
    BamRecordReader reader{path, BamRecordReaderConfig{
                                     .RawReaderConfig = {.BgzfWorkers = bgzfThreads},
                                     .DecodeWorkers = decodeThreads,
                                 }};
    std::size_t records{0};

    for (const auto& rec : reader.Records()) {
        ++records;
        (void)rec.Name();
    }

    const double secs{ElapsedSecs(start)};

    std::println(
        "record_reader(bgzf={}, decode={}): {} records, {:.3f} sec, "
        "{:.0f} records/s",
        bgzfThreads, decodeThreads, records, secs, Tools::ToRate(records, secs));

    const ReaderMetrics m{reader.GetMetrics()};
    if (m.ParallelBgzf) {
        PrintBgzfMetrics(m.Bgzf);
    }
    if (m.ParallelDecode) {
        PrintDecodeMetrics(m.Decode);
    }
}

void BenchRegionQuery(const std::filesystem::path& bamPath)
{
    const std::filesystem::path baiPath{SidecarPath(bamPath, ".bai")};
    if (!std::filesystem::exists(baiPath)) {
        std::println("region_query: skipped (no .bai)");
        return;
    }

    const auto index{BaiIndex::FromFile(baiPath)};
    BamRawReader reader{bamPath};
    const auto& header{reader.Header()};

    if (header.NumReferences() == 0) {
        std::println("region_query: skipped (no references)");
        return;
    }

    const std::int32_t refLen{header.ReferenceSequences()[0].Length()};
    const std::int32_t queryLen{std::ranges::min(refLen, std::int32_t{100000})};

    const auto start{std::chrono::steady_clock::now()};
    std::size_t records{0};
    const std::int32_t step{std::ranges::max(queryLen / 10, std::int32_t{1})};

    for (std::int32_t beg{0}; beg < refLen; beg += step) {
        const std::int32_t end{std::ranges::min(beg + queryLen, refLen)};
        for (const auto& view : reader.Query(index, 0, beg, end)) {
            ++records;
            (void)view.Name();
        }
    }

    const double secs{ElapsedSecs(start)};

    std::println(
        "region_query: {} records from sliding window, {:.3f} sec, "
        "{:.0f} records/s",
        records, secs, Tools::ToRate(records, secs));
}

}  // namespace

CLI_v2::Interface CreateInterface()
{
    CLI_v2::Interface interface{"pbsamoa bench",
                                "Benchmark pbsamoa read/write performance on a BAM file.",
                                LibraryFormattedVersion()};
    // bench keeps its own -j/--bgzf-threads (default -1 = auto, capped at 10);
    // the built-in --num-threads/-j resolves 0 to raw hardware count and would
    // shadow -j, so disable it.
    interface.DisableNumThreadsOption();
    interface.AddOptions({BgzfThreads, DecodeThreads});
    interface.AddPositionalArguments({Input});
    return interface;
}

int Runner(const CLI_v2::Results& results)
{
    // Read into exact-width ints (never std::size_t — ambiguous conversion).
    const std::int32_t bgzfOpt{results[BgzfThreads]};
    const std::int32_t decodeOpt{results[DecodeThreads]};

    // CLIv2 does not enforce the required positional-argument count.
    const std::vector<std::string>& pos{results.PositionalArguments()};
    if (pos.size() != 1) {
        throw std::runtime_error{"bench requires exactly one argument: INPUT"};
    }

    const std::filesystem::path path{pos[0]};

    const std::size_t bgzfWorkers{Tools::ResolveNumWorkers(bgzfOpt, /*explicitCap=*/10)};
    std::size_t decodeWorkers{4U};
    if (decodeOpt >= 0) {
        decodeWorkers = static_cast<std::size_t>(decodeOpt);
    }

    std::println("=== pbsamoa benchmark: {} ===", path.filename().string());
    std::println("    bgzf-threads={}, decode-threads={}\n", bgzfWorkers, decodeWorkers);

    std::println("--- BamRawReader (sync) ---");
    BenchSequentialRead(path);
    BenchBatchRead(path, ByteLimit{64U * 1024U * 1024U});

    if (bgzfWorkers > 0) {
        std::println("\n--- BamRawReader (pipeline, bgzf={}) ---", bgzfWorkers);
        BenchPipelineRead(path, bgzfWorkers);
    }

    std::println("\n--- BamRecordReader (serial) ---");
    BenchRecordReader(path, 0, 0);

    if ((bgzfWorkers > 0) || (decodeWorkers > 0)) {
        std::println("\n--- BamRecordReader (bgzf={}, decode={}) ---\n", bgzfWorkers,
                     decodeWorkers);
        BenchRecordReader(path, bgzfWorkers, decodeWorkers);
    }

    std::println("\n--- Write ---");
    BenchWrite(path);
    std::println("\n--- Region Query ---");
    BenchRegionQuery(path);

    std::println("\n=== done ===");

    return EXIT_SUCCESS;
}

}  // namespace Bench
}  // namespace Samoa
}  // namespace PacBio
