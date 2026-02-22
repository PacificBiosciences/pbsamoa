#include "Dump.hpp"

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/Metrics.hpp>
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/Sequence.hpp>
#include <pbsamoa/core/Tags.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/SamReader.hpp>

#include <parallel/ThreadPool.h>

#include <algorithm>
#include <array>
#include <charconv>
#include <chrono>
#include <exception>
#include <filesystem>
#include <future>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <thread>

#include <cstdint>
#include <cstdio>
#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace Dump {
namespace {

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

void AppendInt(std::string& out, std::int64_t v)
{
    std::array<char, 24> buf{};
    const auto [ptr, ec]{std::to_chars(std::data(buf), std::data(buf) + std::size(buf), v)};
    out.append(std::data(buf), ptr);
}

void FormatRecord(const RawRecord& view, const SamHeader& header, std::string& buf)
{
    buf.clear();

    buf.append(view.Name());
    buf += '\t';

    AppendInt(buf, view.Flag());
    buf += '\t';

    const std::int32_t refId{view.RefId()};
    if (refId < 0) {
        buf += '*';
    } else {
        buf.append(header.ReferenceName(refId));
    }
    buf += '\t';

    const std::int32_t pos{view.Pos()};
    AppendInt(buf, (pos >= 0) ? (pos + 1) : 0);
    buf += '\t';

    AppendInt(buf, view.MapQ());
    buf += '\t';

    WriteCigarTo(view.CigarOps(), buf);
    buf += '\t';

    const std::int32_t nextRefId{view.NextRefId()};
    if (nextRefId < 0) {
        buf += '*';
    } else if (nextRefId == refId) {
        buf += '=';
    } else {
        buf.append(header.ReferenceName(nextRefId));
    }
    buf += '\t';

    const std::int32_t nextPos{view.NextPos()};
    AppendInt(buf, (nextPos >= 0) ? (nextPos + 1) : 0);
    buf += '\t';

    AppendInt(buf, view.Tlen());
    buf += '\t';

    if (view.SeqLength() == 0) {
        buf += '*';
    } else {
        view.Seq().WriteTo(buf);
    }
    buf += '\t';

    const std::span<const std::uint8_t> qual{view.Qual()};
    const bool qualUnavailable{std::empty(qual) ||
                               std::ranges::all_of(qual, [](std::uint8_t q) { return q == 0xFF; })};
    if (qualUnavailable) {
        buf += '*';
    } else {
        const std::size_t startPos{std::size(buf)};
        buf.resize(startPos + std::size(qual));
        for (std::size_t qi{0}; qi < std::size(qual); ++qi) {
            buf[startPos + qi] = static_cast<char>(qual[qi] + 33);
        }
    }

    SerializeRawTagsToSam(view.AuxData(), buf);
    buf += '\n';
}

void FormatRecord(const BamRecord& record, const SamHeader& header, std::string& buf)
{
    buf.clear();

    buf.append(record.Name());
    buf += '\t';

    AppendInt(buf, record.Flag());
    buf += '\t';

    const std::int32_t refId{record.RefId()};
    if (refId < 0) {
        buf += '*';
    } else {
        buf.append(header.ReferenceName(refId));
    }
    buf += '\t';

    const std::int32_t pos{record.Pos()};
    AppendInt(buf, (pos >= 0) ? (pos + 1) : 0);
    buf += '\t';

    AppendInt(buf, record.MapQ());
    buf += '\t';

    WriteCigarTo(record.Cigar(), buf);
    buf += '\t';

    const std::int32_t nextRefId{record.NextRefId()};
    if (nextRefId < 0) {
        buf += '*';
    } else if (nextRefId == refId) {
        buf += '=';
    } else {
        buf.append(header.ReferenceName(nextRefId));
    }
    buf += '\t';

    const std::int32_t nextPos{record.NextPos()};
    AppendInt(buf, (nextPos >= 0) ? (nextPos + 1) : 0);
    buf += '\t';

    AppendInt(buf, record.Tlen());
    buf += '\t';

    if (std::empty(record.Sequence())) {
        buf += '*';
    } else {
        buf.append(record.Sequence());
    }
    buf += '\t';

    const std::span<const std::uint8_t> qual{record.Qualities()};
    const bool qualUnavailable{std::empty(qual) ||
                               std::ranges::all_of(qual, [](std::uint8_t q) { return q == 0xFF; })};
    if (qualUnavailable) {
        buf += '*';
    } else {
        const std::size_t startPos{std::size(buf)};
        buf.resize(startPos + std::size(qual));
        for (std::size_t qi{0}; qi < std::size(qual); ++qi) {
            buf[startPos + qi] = static_cast<char>(qual[qi] + 33);
        }
    }

    for (const auto& [key, value] : record.Tags().Entries()) {
        buf += '\t';
        buf.append(SerializeTagToSam(key, value));
    }

    buf += '\n';
}

void PrintMetricsLine(const BgzfMetrics& m, const BgzfMetrics& prev, double elapsedSec)
{
    const double mbRead{static_cast<double>(m.BytesRead) / (1024.0 * 1024.0)};
    const double mbDecomp{static_cast<double>(m.BytesDecompressed) / (1024.0 * 1024.0)};
    const std::uint64_t deltaRecords{m.RecordsConsumed - prev.RecordsConsumed};
    const double recPerSec{(elapsedSec > 0) ? (static_cast<double>(deltaRecords) / elapsedSec)
                                            : 0.0};
    const double deltaMbDecomp{static_cast<double>(m.BytesDecompressed - prev.BytesDecompressed) /
                               (1024.0 * 1024.0)};
    const double mbPerSec{(elapsedSec > 0) ? (deltaMbDecomp / elapsedSec) : 0.0};
    const std::uint64_t queueDepth{m.RecordsProduced - m.RecordsConsumed};

    std::fprintf(stderr,
                 "[metrics] %.1f MB in, %.1f MB out | "
                 "%.0f rec/s, %.1f MB/s | "
                 "pool: %zu active, queue %zu | "
                 "spsc: depth %llu | "
                 "stalls: io=%llu cons=%llu read=%llu\n",
                 mbRead, mbDecomp, recPerSec, mbPerSec, m.PoolActiveWorkers, m.PoolQueueDepth,
                 static_cast<unsigned long long>(queueDepth),
                 static_cast<unsigned long long>(m.IoStalls),
                 static_cast<unsigned long long>(m.ConsumerStalls),
                 static_cast<unsigned long long>(m.ReaderStalls));
}

void PrintMetricsSummary(const BgzfMetrics& m)
{
    const double mbRead{static_cast<double>(m.BytesRead) / (1024.0 * 1024.0)};
    const double mbDecomp{static_cast<double>(m.BytesDecompressed) / (1024.0 * 1024.0)};
    const double ioMs{static_cast<double>(m.IoReadNs) / 1e6};
    const double decompMs{static_cast<double>(m.DecompressNs) / 1e6};
    const double parseMs{static_cast<double>(m.RecordParseNs) / 1e6};

    std::fprintf(stderr, "\n--- Final Pipeline Metrics ---\n");
    std::fprintf(stderr, "BGZF IO:       %.1f MB read, %llu blocks\n", mbRead,
                 static_cast<unsigned long long>(m.BlocksRead));
    std::fprintf(stderr, "Decompressed:  %.1f MB\n", mbDecomp);
    std::fprintf(stderr, "Records:       %llu produced, %llu consumed\n",
                 static_cast<unsigned long long>(m.RecordsProduced),
                 static_cast<unsigned long long>(m.RecordsConsumed));
    std::fprintf(stderr, "Pool peak:     queue %zu, workers %zu\n", m.PoolPeakQueueDepth,
                 m.PoolPeakActiveWorkers);
    std::fprintf(stderr, "Result peak:   queue %zu\n", m.PoolPeakResultQueueDepth);
    std::fprintf(stderr, "Stalls:        IO=%llu, consumer=%llu, reader=%llu\n",
                 static_cast<unsigned long long>(m.IoStalls),
                 static_cast<unsigned long long>(m.ConsumerStalls),
                 static_cast<unsigned long long>(m.ReaderStalls));
    std::fprintf(stderr, "Timing:        IO %.1fms, decompress %.1fms, parse %.1fms\n", ioMs,
                 decompMs, parseMs);
}

void DumpSam(const std::filesystem::path& path)
{
    static std::array<char, 1 << 20> stdoutBuf;
    std::setvbuf(stdout, std::data(stdoutBuf), _IOFBF, std::size(stdoutBuf));

    SamReader reader{path};
    const SamHeader& header{reader.Header()};

    const std::string headerText{header.ToText()};
    std::fwrite(std::data(headerText), 1, std::size(headerText), stdout);

    std::string buf;
    for (const BamRecord& record : reader.Records()) {
        FormatRecord(record, header, buf);
        std::fwrite(std::data(buf), 1, std::size(buf), stdout);
    }
}

bool MetricsEnabled()
{
    const char* env{std::getenv("PBSAMOA_METRICS")};
    return (env != nullptr) && (std::string_view{env} == "1");
}

void DumpBam(const std::filesystem::path& path, std::size_t numWorkers,
             std::size_t numFormatThreads)
{
    static std::array<char, 1 << 20> stdoutBuf;
    std::setvbuf(stdout, std::data(stdoutBuf), _IOFBF, std::size(stdoutBuf));

    const bool showMetrics{MetricsEnabled()};

    // Pipeline handles BGZF I/O + decompression with its own threads
    BamRawReader reader{path, BamRawReaderConfig{.BgzfWorkers = numWorkers}};

    const auto& header{reader.Header()};

    const std::string headerText{header.ToText()};
    std::fwrite(std::data(headerText), 1, std::size(headerText), stdout);

    // Metrics thread: print every 1s to stderr (only when PBSAMOA_METRICS=1)
    std::jthread metricsThread;
    if (showMetrics) {
        metricsThread = std::jthread{[&reader](std::stop_token stopToken) {
            using namespace std::chrono_literals;
            BgzfMetrics prev{};
            auto lastTime{std::chrono::steady_clock::now()};

            while (!stopToken.stop_requested()) {
                std::this_thread::sleep_for(1s);
                if (stopToken.stop_requested()) {
                    break;
                }

                const auto now{std::chrono::steady_clock::now()};
                const double elapsed{
                    std::chrono::duration_cast<std::chrono::duration<double>>(now - lastTime)
                        .count()};
                const BgzfMetrics current{reader.GetMetrics()};
                PrintMetricsLine(current, prev, elapsed);
                prev = current;
                lastTime = now;
            }
        }};
    }

    PacBio::Parallel::ThreadPool<std::string> formatPool{
        PacBio::Parallel::ThreadPool<std::string>::Config{.NumThreads = numFormatThreads,
                                                          .QueueMultiplier = 5}};
    std::exception_ptr consumerException;
    std::jthread consumerThread{[&formatPool, &consumerException](std::stop_token) {
        try {
            while (formatPool.ConsumeWith([](std::string output) {
                std::fwrite(std::data(output), 1, std::size(output), stdout);
            })) {}
        } catch (...) {
            consumerException = std::current_exception();
        }
    }};

    try {
        // Double-buffer: read next batch while formatting current
        auto currentBatch{reader.ReadBatch()};

        while (currentBatch) {
            // Start reading next batch concurrently
            auto nextBatchFuture{
                std::async(std::launch::async, [&reader]() { return reader.ReadBatch(); })};

            const auto batch{std::make_shared<RawRecordBatch>(std::move(*currentBatch))};
            for (std::size_t recordIdx{0}; recordIdx < batch->RecordCount(); ++recordIdx) {
                formatPool.Submit([batch, &header, recordIdx]() -> std::string {
                    std::string output;
                    const RawRecord view{batch->RecordData(recordIdx)};
                    FormatRecord(view, header, output);
                    return output;
                });
            }

            currentBatch = nextBatchFuture.get();
        }

        formatPool.Finalize();
        consumerThread.join();
        if (showMetrics) {
            metricsThread.request_stop();
            metricsThread.join();
            PrintMetricsSummary(reader.GetMetrics());
        }
    } catch (...) {
        if (metricsThread.joinable()) {
            metricsThread.request_stop();
        }
        try {
            formatPool.Finalize();
        } catch (...) {
        }
        if (consumerThread.joinable()) {
            consumerThread.join();
        }
        if (metricsThread.joinable()) {
            metricsThread.join();
        }
        if (consumerException) {
            std::rethrow_exception(consumerException);
        }
        throw;
    }
    if (consumerException) {
        std::rethrow_exception(consumerException);
    }
}

}  // namespace

int Runner(int argc, char* argv[])
{
    // Parse options
    std::int32_t bgzfOpt{-1};
    std::int32_t formatOpt{0};
    const char* inputFile{nullptr};

    for (int i{0}; i < argc; ++i) {
        const std::string_view arg{argv[i]};
        if ((arg == "--bgzf-threads") && (i + 1 < argc)) {
            bgzfOpt = std::stoi(std::string{argv[++i]});
        } else if ((arg == "--format-threads") && (i + 1 < argc)) {
            formatOpt = std::stoi(std::string{argv[++i]});
        } else if ((arg == "-j") && (i + 1 < argc)) {
            bgzfOpt = std::stoi(std::string{argv[++i]});
        } else if (arg[0] != '-') {
            inputFile = argv[i];
        }
    }

    if (inputFile == nullptr) {
        std::fprintf(stderr, "Usage: pbsamoa dump [--bgzf-threads N] [--format-threads N] INPUT\n");
        return EXIT_FAILURE;
    }

    const std::filesystem::path path{inputFile};
    if (path.extension() == ".sam") {
        DumpSam(path);
        return EXIT_SUCCESS;
    }

    const std::size_t bgzfWorkers{ResolveNumWorkers(bgzfOpt)};

    const std::size_t formatThreads{
        (formatOpt > 0)
            ? static_cast<std::size_t>(formatOpt)
            : static_cast<std::size_t>(std::ranges::max(std::thread::hardware_concurrency(), 4U))};

    DumpBam(path, bgzfWorkers, formatThreads);
    return EXIT_SUCCESS;
}

}  // namespace Dump
}  // namespace Samoa
}  // namespace PacBio
