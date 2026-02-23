#include "Bench.hpp"
#include "../ParseUtils.hpp"

#include <pbsamoa/core/Metrics.hpp>
#include <pbsamoa/index/BaiIndex.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamRecordReader.hpp>
#include <pbsamoa/io/BamWriter.hpp>
#include <print>

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>

#include <cstdint>
#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace Bench {
namespace {

double ElapsedSecs(std::chrono::steady_clock::time_point start)
{
    const auto elapsed{std::chrono::steady_clock::now() - start};
    return std::chrono::duration<double>(elapsed).count();
}

void BenchSequentialRead(const std::filesystem::path& path)
{
    const auto start{std::chrono::steady_clock::now()};
    BamRawReader reader{path};
    std::size_t records{0};
    std::size_t bytes{0};

    for (const auto& view : reader.Records()) {
        ++records;
        bytes += std::size(view.RawData());
    }

    const double secs{ElapsedSecs(start)};
    const double mbPerSec{bytes / (1024.0 * 1024.0) / secs};

    std::println(
        "sequential_read: {} records, {:.1f} MB, {:.3f} sec, {:.1f} MB/s, {:.0f} records/s",
        records, bytes / (1024.0 * 1024.0), secs, mbPerSec, 1.0 * records / secs);
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

    std::println("batch_read({} KiB): {} records, {} batches, {:.3f} sec, {:.0f} records/s",
                 limit.Value() / 1024, records, batches, secs, 1.0 * records / secs);
}

void BenchWrite(const std::filesystem::path& srcPath)
{
    const auto tmpPath{std::filesystem::temp_directory_path() / "pbsamoa_bench_write.bam"};

    BamRawReader srcReader{srcPath};

    const auto batch{srcReader.ReadBatch(ByteLimit{256U * 1024U * 1024U})};
    if (!batch.has_value()) {
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

    std::println("write: {} records, {:.1f} MB output, {:.3f} sec, {:.1f} MB/s, {:.0f} records/s",
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

    std::println("    bgzf: {:.1f} MB compressed, {:.1f} MB decompressed, {} blocks", mbRead,
                 mbDecomp, m.BlocksRead);
    std::println("    timing: io {:.1f}ms, decomp {:.1f}ms, parse {:.1f}ms", ioMs, decompMs,
                 parseMs);
    std::println("    pool peak: queue {}, workers {}, results {}", m.PoolPeakQueueDepth,
                 m.PoolPeakActiveWorkers, m.PoolPeakResultQueueDepth);
    std::println("    stalls: io={}, consumer={}, reader={}", m.IoStalls, m.ConsumerStalls,
                 m.ReaderStalls);
}

void PrintDecodeMetrics(const DecodeMetrics& m)
{
    const double decodeMs{static_cast<double>(m.DecodeNs) / 1e6};
    const double batchReadMs{static_cast<double>(m.BatchReadNs) / 1e6};

    std::println("    decode: {} batches, {} records", m.BatchesDecoded, m.RecordsDecoded);
    std::println("    timing: decode {:.1f}ms, batch_read {:.1f}ms", decodeMs, batchReadMs);
    std::println("    pool peak: queue {}, workers {}", m.PoolPeakQueueDepth,
                 m.PoolPeakActiveWorkers);
    std::println("    stalls: producer={}, consumer={}", m.ProducerStalls, m.ConsumerStalls);
}

void BenchPipelineRead(const std::filesystem::path& path, std::size_t bgzfThreads)
{
    const auto start{std::chrono::steady_clock::now()};
    BamRawReader reader{path, BamRawReaderConfig{.BgzfWorkers = bgzfThreads}};
    std::size_t records{0};
    std::size_t bytes{0};

    for (const auto& view : reader.Records()) {
        ++records;
        bytes += std::size(view.RawData());
    }

    const double secs{ElapsedSecs(start)};
    const double mbPerSec{bytes / (1024.0 * 1024.0) / secs};

    std::println(
        "pipeline_raw_read(bgzf={}): {} records, {:.1f} MB, {:.3f} sec, {:.1f} MB/s, {:.0f} "
        "records/s",
        bgzfThreads, records, bytes / (1024.0 * 1024.0), secs, mbPerSec, 1.0 * records / secs);

    PrintBgzfMetrics(reader.GetMetrics());
}

void BenchRecordReader(const std::filesystem::path& path, std::size_t bgzfThreads,
                       std::size_t decodeThreads)
{
    const auto start{std::chrono::steady_clock::now()};
    BamRecordReader reader{path, BamRecordReaderConfig{
                                     .ViewConfig = {.BgzfWorkers = bgzfThreads},
                                     .DecodeWorkers = decodeThreads,
                                 }};
    std::size_t records{0};

    for (const auto& rec : reader.Records()) {
        ++records;
        (void)rec.Name();
    }

    const double secs{ElapsedSecs(start)};

    std::println("record_reader(bgzf={}, decode={}): {} records, {:.3f} sec, {:.0f} records/s",
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

    std::println("region_query: {} records from sliding window, {:.3f} sec, {:.0f} records/s",
                 records, secs, 1.0 * records / secs);
}

std::size_t ResolveNumWorkers(std::int32_t requested)
{
    constexpr std::int32_t MAX_AUTO_WORKERS{8};
    constexpr std::int32_t MAX_PARSE_WORKERS{10};

    const std::int32_t hwThreads{
        static_cast<std::int32_t>(std::ranges::max(std::thread::hardware_concurrency(), 1U))};

    if (requested < 0) {
        return static_cast<std::size_t>(std::ranges::min(hwThreads, MAX_AUTO_WORKERS));
    }
    return static_cast<std::size_t>(std::ranges::min(requested, MAX_PARSE_WORKERS));
}

}  // namespace

int Runner(int argc, char* argv[])
{
    // Parse options
    std::int32_t bgzfOpt{-1};
    std::int32_t decodeOpt{4};
    const char* inputFile{nullptr};

    for (int i{0}; i < argc; ++i) {
        const std::string_view arg{argv[i]};
        if ((arg == "--bgzf-threads") && (i + 1 < argc)) {
            const auto parsed{Tools::ParseInteger<std::int32_t>(argv[++i], "bgzf-threads")};
            if (!parsed) {
                throw std::runtime_error{parsed.error()};
            }
            bgzfOpt = *parsed;
        } else if ((arg == "--decode-threads") && (i + 1 < argc)) {
            const auto parsed{Tools::ParseInteger<std::int32_t>(argv[++i], "decode-threads")};
            if (!parsed) {
                throw std::runtime_error{parsed.error()};
            }
            decodeOpt = *parsed;
        } else if ((arg == "-j") && (i + 1 < argc)) {
            const auto parsed{Tools::ParseInteger<std::int32_t>(argv[++i], "bgzf-threads")};
            if (!parsed) {
                throw std::runtime_error{parsed.error()};
            }
            bgzfOpt = *parsed;
        } else if (arg[0] != '-') {
            inputFile = argv[i];
        }
    }

    if (inputFile == nullptr) {
        std::println(stderr, "Usage: pbsamoa bench [--bgzf-threads N] [--decode-threads N] INPUT");
        return EXIT_FAILURE;
    }

    const std::filesystem::path path{inputFile};

    const std::size_t bgzfWorkers{ResolveNumWorkers(bgzfOpt)};
    const std::size_t decodeWorkers{(decodeOpt >= 0) ? static_cast<std::size_t>(decodeOpt) : 4};

    std::println("=== pbsamoa benchmark: {} ===", path.filename().string());
    std::println("    bgzf-threads={}, decode-threads={}\n", bgzfWorkers, decodeWorkers);

    // --- BamRawReader benchmarks ---
    std::println("--- BamRawReader (sync) ---");
    BenchSequentialRead(path);
    BenchBatchRead(path, ByteLimit{64U * 1024U * 1024U});

    if (bgzfWorkers > 0) {
        std::println("\n--- BamRawReader (pipeline, bgzf={}) ---", bgzfWorkers);
        BenchPipelineRead(path, bgzfWorkers);
    }

    // --- BamRecordReader benchmarks ---
    std::println("\n--- BamRecordReader (serial) ---");
    BenchRecordReader(path, 0, 0);

    if ((bgzfWorkers > 0) || (decodeWorkers > 0)) {
        std::println("\n--- BamRecordReader (bgzf={}, decode={}) ---\n", bgzfWorkers,
                     decodeWorkers);
        BenchRecordReader(path, bgzfWorkers, decodeWorkers);
    }

    // --- Write + region query ---
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
