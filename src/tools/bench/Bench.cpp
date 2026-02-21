#include "Bench.hpp"

#include <pbsamoa/core/Metrics.hpp>
#include <pbsamoa/index/BaiIndex.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamRecordReader.hpp>
#include <pbsamoa/io/BamWriter.hpp>

#include <pbcopper/cli2/CLI.h>
#include <pbcopper/cli2/internal/BuiltinOptions.h>
#include <pbcopper/utility/Stopwatch.h>

#include <algorithm>
#include <filesystem>
#include <span>

#include <cstdint>
#include <cstdio>
#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace Bench {
namespace {

double ElapsedSecs(const PacBio::Utility::Stopwatch& timer)
{
    return timer.ElapsedNanoseconds() / 1.0e9;
}

void BenchSequentialRead(const std::filesystem::path& path)
{
    const PacBio::Utility::Stopwatch timer;
    BamRawReader reader{path};
    std::size_t records{0};
    std::size_t bytes{0};

    for (const auto& view : reader.Records()) {
        ++records;
        bytes += std::size(view.RawData());
    }

    const double secs{ElapsedSecs(timer)};
    const double mbPerSec{bytes / (1024.0 * 1024.0) / secs};

    std::fprintf(stdout,
                 "sequential_read: %zu records, %.1f MB, %.3f sec, "
                 "%.1f MB/s, %.0f records/s\n",
                 records, bytes / (1024.0 * 1024.0), secs, mbPerSec, 1.0 * records / secs);
}

void BenchBatchRead(const std::filesystem::path& path, ByteLimit limit)
{
    const PacBio::Utility::Stopwatch timer;
    BamRawReader reader{path};
    std::size_t records{0};
    std::size_t batches{0};

    while (const auto batch = reader.ReadBatch(limit)) {
        records += batch->RecordCount();
        ++batches;
    }

    const double secs{ElapsedSecs(timer)};

    std::fprintf(stdout,
                 "batch_read(%zu KiB): %zu records, %zu batches, "
                 "%.3f sec, %.0f records/s\n",
                 limit.Value() / 1024, records, batches, secs, 1.0 * records / secs);
}

void BenchWrite(const std::filesystem::path& srcPath)
{
    const auto tmpPath{std::filesystem::temp_directory_path() / "pbsamoa_bench_write.bam"};

    BamRawReader srcReader{srcPath};

    const auto batch{srcReader.ReadBatch(ByteLimit{256U * 1024U * 1024U})};
    if (!batch.has_value()) {
        std::fprintf(stdout, "write: no records to benchmark\n");
        return;
    }
    const std::size_t numRecords{batch->RecordCount()};

    const PacBio::Utility::Stopwatch timer;
    {
        BamWriter writer{tmpPath, srcReader.Header()};
        for (std::size_t i{0}; i < numRecords; ++i) {
            writer.Write(batch->RecordData(i));
        }
    }
    const double secs{ElapsedSecs(timer)};

    const auto fileSize{std::filesystem::file_size(tmpPath)};

    std::fprintf(stdout,
                 "write: %zu records, %.1f MB output, %.3f sec, "
                 "%.1f MB/s, %.0f records/s\n",
                 numRecords, fileSize / (1024.0 * 1024.0), secs,
                 fileSize / (1024.0 * 1024.0) / secs, 1.0 * numRecords / secs);

    std::filesystem::remove(tmpPath);
}

void PrintBgzfMetrics(const BgzfMetrics& m)
{
    const double mbRead{static_cast<double>(m.BytesRead) / (1024.0 * 1024.0)};
    const double mbDecomp{static_cast<double>(m.BytesDecompressed) / (1024.0 * 1024.0)};
    const double ioMs{static_cast<double>(m.IoReadNs) / 1e6};
    const double decompMs{static_cast<double>(m.DecompressNs) / 1e6};
    const double parseMs{static_cast<double>(m.RecordParseNs) / 1e6};

    std::fprintf(stdout,
                 "    bgzf: %.1f MB compressed, %.1f MB decompressed, "
                 "%llu blocks\n",
                 mbRead, mbDecomp, static_cast<unsigned long long>(m.BlocksRead));
    std::fprintf(stdout, "    timing: io %.1fms, decomp %.1fms, parse %.1fms\n", ioMs, decompMs,
                 parseMs);
    std::fprintf(stdout, "    pool peak: queue %zu, workers %zu, results %zu\n",
                 m.PoolPeakQueueDepth, m.PoolPeakActiveWorkers, m.PoolPeakResultQueueDepth);
    std::fprintf(stdout, "    stalls: io=%llu, consumer=%llu, reader=%llu\n",
                 static_cast<unsigned long long>(m.IoStalls),
                 static_cast<unsigned long long>(m.ConsumerStalls),
                 static_cast<unsigned long long>(m.ReaderStalls));
}

void PrintDecodeMetrics(const DecodeMetrics& m)
{
    const double decodeMs{static_cast<double>(m.DecodeNs) / 1e6};
    const double batchReadMs{static_cast<double>(m.BatchReadNs) / 1e6};

    std::fprintf(stdout, "    decode: %llu batches, %llu records\n",
                 static_cast<unsigned long long>(m.BatchesDecoded),
                 static_cast<unsigned long long>(m.RecordsDecoded));
    std::fprintf(stdout, "    timing: decode %.1fms, batch_read %.1fms\n", decodeMs, batchReadMs);
    std::fprintf(stdout, "    pool peak: queue %zu, workers %zu\n", m.PoolPeakQueueDepth,
                 m.PoolPeakActiveWorkers);
    std::fprintf(stdout, "    stalls: producer=%llu, consumer=%llu\n",
                 static_cast<unsigned long long>(m.ProducerStalls),
                 static_cast<unsigned long long>(m.ConsumerStalls));
}

void BenchPipelineRead(const std::filesystem::path& path, std::size_t bgzfThreads)
{
    const PacBio::Utility::Stopwatch timer;
    BamRawReader reader{path, BamRawReaderConfig{.BgzfWorkers = bgzfThreads}};
    std::size_t records{0};
    std::size_t bytes{0};

    for (const auto& view : reader.Records()) {
        ++records;
        bytes += std::size(view.RawData());
    }

    const double secs{ElapsedSecs(timer)};
    const double mbPerSec{bytes / (1024.0 * 1024.0) / secs};

    std::fprintf(stdout,
                 "pipeline_raw_read(bgzf=%zu): %zu records, %.1f MB, %.3f sec, "
                 "%.1f MB/s, %.0f records/s\n",
                 bgzfThreads, records, bytes / (1024.0 * 1024.0), secs, mbPerSec,
                 1.0 * records / secs);

    PrintBgzfMetrics(reader.GetMetrics());
}

void BenchRecordReader(const std::filesystem::path& path, std::size_t bgzfThreads,
                       std::size_t decodeThreads)
{
    const PacBio::Utility::Stopwatch timer;
    BamRecordReader reader{path, BamRecordReaderConfig{
                                     .ViewConfig = {.BgzfWorkers = bgzfThreads},
                                     .DecodeWorkers = decodeThreads,
                                 }};
    std::size_t records{0};

    for (const auto& rec : reader.Records()) {
        ++records;
        (void)rec.Name();
    }

    const double secs{ElapsedSecs(timer)};

    std::fprintf(stdout,
                 "record_reader(bgzf=%zu, decode=%zu): %zu records, %.3f sec, "
                 "%.0f records/s\n",
                 bgzfThreads, decodeThreads, records, secs, 1.0 * records / secs);

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
    const auto baiPath{std::filesystem::path{bamPath.string() + ".bai"}};
    if (!std::filesystem::exists(baiPath)) {
        std::fprintf(stdout, "region_query: skipped (no .bai)\n");
        return;
    }

    const auto index{BaiIndex::FromFile(baiPath)};
    BamRawReader reader{bamPath};
    const auto& header{reader.Header()};

    if (header.NumReferences() == 0) {
        std::fprintf(stdout, "region_query: skipped (no references)\n");
        return;
    }

    const std::int32_t refLen{header.ReferenceSequences()[0].Length()};
    const std::int32_t queryLen{std::ranges::min(refLen, std::int32_t{100000})};

    const PacBio::Utility::Stopwatch timer;
    std::size_t records{0};
    const std::int32_t step{std::ranges::max(queryLen / 10, std::int32_t{1})};

    for (std::int32_t beg{0}; beg < refLen; beg += step) {
        const std::int32_t end{std::ranges::min(beg + queryLen, refLen)};
        for (const auto& view : reader.Query(index, 0, beg, end)) {
            ++records;
            (void)view.Name();
        }
    }

    const double secs{ElapsedSecs(timer)};

    std::fprintf(stdout,
                 "region_query: %zu records from sliding window, "
                 "%.3f sec, %.0f records/s\n",
                 records, secs, 1.0 * records / secs);
}

std::size_t ResolveNumWorkers(const PacBio::CLI_v2::Results& results)
{
    constexpr std::int32_t MAX_AUTO_WORKERS{8};
    constexpr std::int32_t MAX_PARSE_WORKERS{10};

    const std::int32_t requestedNumThreads{results[PacBio::CLI_v2::Builtin::NumThreads]};
    const std::int32_t normalizedNumThreads{results.NumThreads()};

    if (requestedNumThreads == 0) {
        return static_cast<std::size_t>(std::ranges::min(normalizedNumThreads, MAX_AUTO_WORKERS));
    }
    return static_cast<std::size_t>(std::ranges::min(normalizedNumThreads, MAX_PARSE_WORKERS));
}

// clang-format off
const PacBio::CLI_v2::PositionalArgument BenchInput{
R"({
    "name"        : "input",
    "description" : "Input BAM file.",
    "type"        : "file"
})"};

const PacBio::CLI_v2::Option BgzfThreads{
R"({
    "names" : ["bgzf-threads"],
    "description" : "Number of BGZF decompression threads. 0 = synchronous.",
    "type" : "integer",
    "default" : -1
})"};

const PacBio::CLI_v2::Option DecodeThreads{
R"({
    "names" : ["decode-threads"],
    "description" : "Number of BamRecord decode threads for record_reader bench. 0 = serial.",
    "type" : "integer",
    "default" : 4
})"};
// clang-format on

}  // namespace

PacBio::CLI_v2::Interface CreateInterface()
{
    PacBio::CLI_v2::Interface iface{"bench", "Benchmark pbsamoa read/write performance", "0.1.0"};
    iface.DisableLogFileOption();
    iface.AddPositionalArgument(BenchInput);
    iface.AddOptionGroup("Threading", {BgzfThreads, DecodeThreads});
    return iface;
}

int Runner(const PacBio::CLI_v2::Results& results)
{
    const std::filesystem::path path{results[BenchInput]};

    // Resolve BGZF threads: --bgzf-threads overrides -j
    const std::int32_t bgzfOpt{results[BgzfThreads]};
    const std::size_t bgzfWorkers{(bgzfOpt >= 0) ? static_cast<std::size_t>(bgzfOpt)
                                                 : ResolveNumWorkers(results)};

    const std::int32_t decodeOpt{results[DecodeThreads]};
    const std::size_t decodeWorkers{(decodeOpt >= 0) ? static_cast<std::size_t>(decodeOpt) : 4};

    std::fprintf(stdout, "=== pbsamoa benchmark: %s ===\n", path.filename().c_str());
    std::fprintf(stdout, "    bgzf-threads=%zu, decode-threads=%zu\n\n", bgzfWorkers,
                 decodeWorkers);

    // --- BamRawReader benchmarks ---
    std::fprintf(stdout, "--- BamRawReader (sync) ---\n");
    BenchSequentialRead(path);
    // BenchBatchRead(path, ByteLimit{64U * 1024U});
    // BenchBatchRead(path, ByteLimit{1024U * 1024U});
    BenchBatchRead(path, ByteLimit{64U * 1024U * 1024U});

    if (bgzfWorkers > 0) {
        std::fprintf(stdout, "\n--- BamRawReader (pipeline, bgzf=%zu) ---\n", bgzfWorkers);
        BenchPipelineRead(path, bgzfWorkers);
    }

    // --- BamRecordReader benchmarks ---
    std::fprintf(stdout, "\n--- BamRecordReader (serial) ---\n");
    BenchRecordReader(path, 0, 0);

    if ((bgzfWorkers > 0) || (decodeWorkers > 0)) {
        std::fprintf(stdout, "\n--- BamRecordReader (bgzf=%zu, decode=%zu) ---\n", bgzfWorkers,
                     decodeWorkers);
        BenchRecordReader(path, bgzfWorkers, decodeWorkers);
    }

    // --- Write + region query ---
    std::fprintf(stdout, "\n--- Write ---\n");
    BenchWrite(path);
    std::fprintf(stdout, "\n--- Region Query ---\n");
    BenchRegionQuery(path);

    std::fprintf(stdout, "\n=== done ===\n");

    return EXIT_SUCCESS;
}

}  // namespace Bench
}  // namespace Samoa
}  // namespace PacBio
