#include "Dump.hpp"
#include "../../CramInternal.hpp"
#include "../ParseUtils.hpp"

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/Metrics.hpp>
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/Sequence.hpp>
#include <pbsamoa/core/Tags.hpp>
#include <pbsamoa/index/CraiIndex.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/CramReader.hpp>
#include <pbsamoa/io/SamReader.hpp>
#include <print>

#include <parallel/ThreadPool.h>

#include <algorithm>
#include <array>
#include <charconv>
#include <chrono>
#include <exception>
#include <expected>
#include <filesystem>
#include <future>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
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
    return Tools::ResolveNumWorkers(requested, /*explicitCap=*/10);
}

struct CramRegion
{
    bool Unmapped{false};
    std::string RefName;
    std::int32_t Beg{0};  // 0-based inclusive
    std::int32_t End{0};  // 0-based exclusive
};

std::expected<CramRegion, std::string> ParseCramRegion(std::string_view text)
{
    if (text == "*") {
        return CramRegion{.Unmapped = true, .Beg = 0, .End = 0};
    }

    const std::size_t colonPos{text.find(':')};
    if (colonPos == std::string_view::npos) {
        return std::unexpected{
            std::format("invalid region format '{}' (expected ref:start-end or '*')", text)};
    }

    const std::string refName{text.substr(0, colonPos)};
    const std::string_view rest{text.substr(colonPos + 1)};
    const std::size_t dashPos{rest.find('-')};
    if (dashPos == std::string_view::npos) {
        return std::unexpected{
            std::format("invalid region format '{}' (expected ref:start-end or '*')", text)};
    }

    return Tools::ParseInteger<std::int32_t>(rest.substr(0, dashPos), "region start")
        .and_then([&](std::int32_t start) {
            return Tools::ParseInteger<std::int32_t>(rest.substr(dashPos + 1), "region end")
                .and_then([&](std::int32_t end) -> std::expected<CramRegion, std::string> {
                    if ((start < 1) || (end < start)) {
                        return std::unexpected{"region must satisfy start >= 1 and end >= start"};
                    }
                    return CramRegion{
                        .Unmapped = false,
                        .RefName = refName,
                        .Beg = start - 1,
                        .End = end,
                    };
                });
        });
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
    const double mbRead = m.BytesRead / (1024.0 * 1024.0);
    const double mbDecomp = m.BytesDecompressed / (1024.0 * 1024.0);
    const double deltaRecords = m.RecordsConsumed - prev.RecordsConsumed;
    const double recPerSec = (elapsedSec > 0) ? (deltaRecords / elapsedSec) : 0.0;
    const double deltaMbDecomp = (m.BytesDecompressed - prev.BytesDecompressed) / (1024.0 * 1024.0);
    const double mbPerSec{(elapsedSec > 0) ? (deltaMbDecomp / elapsedSec) : 0.0};
    const std::uint64_t queueDepth{m.RecordsProduced - m.RecordsConsumed};

    std::println(stderr,
                 "[metrics] {:.1f} MB in, {:.1f} MB out | {:.0f} rec/s, {:.1f} "
                 "MB/s | pool: {} active, "
                 "queue {} | spsc: depth {} | stalls: io={} cons={} read={}",
                 mbRead, mbDecomp, recPerSec, mbPerSec, m.PoolActiveWorkers, m.PoolQueueDepth,
                 queueDepth, m.IoStalls, m.ConsumerStalls, m.ReaderStalls);
}

void PrintMetricsSummary(const BgzfMetrics& m)
{
    const double mbRead = m.BytesRead / (1024.0 * 1024.0);
    const double mbDecomp = m.BytesDecompressed / (1024.0 * 1024.0);
    const double ioMs = 1.0 * m.IoReadNs / 1e6;
    const double decompMs = 1.0 * m.DecompressNs / 1e6;
    const double parseMs = 1.0 * m.RecordParseNs / 1e6;

    std::println(stderr, "\n--- Final Pipeline Metrics ---");
    std::println(stderr, "BGZF IO:       {:.1f} MB read, {} blocks", mbRead, m.BlocksRead);
    std::println(stderr, "Decompressed:  {:.1f} MB", mbDecomp);
    std::println(stderr, "Records:       {} produced, {} consumed", m.RecordsProduced,
                 m.RecordsConsumed);
    std::println(stderr, "Pool peak:     queue {}, workers {}", m.PoolPeakQueueDepth,
                 m.PoolPeakActiveWorkers);
    std::println(stderr, "Result peak:   queue {}", m.PoolPeakResultQueueDepth);
    std::println(stderr, "Stalls:        IO={}, consumer={}, reader={}", m.IoStalls,
                 m.ConsumerStalls, m.ReaderStalls);
    std::println(stderr, "Timing:        IO {:.1f}ms, decompress {:.1f}ms, parse {:.1f}ms", ioMs,
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

void DumpCram(const std::filesystem::path& path, const std::filesystem::path& referencePath,
              const std::optional<std::string>& regionText,
              const std::optional<std::filesystem::path>& indexPath, std::size_t numWorkers)
{
    static std::array<char, 1 << 20> stdoutBuf;
    std::setvbuf(stdout, std::data(stdoutBuf), _IOFBF, std::size(stdoutBuf));

    CramReaderConfig config;
    config.ReferencePath = referencePath;
    config.DecompressionWorkers = numWorkers;
    CramReader reader{path, config};
    const SamHeader& header{reader.Header()};

    const std::string headerText{header.ToText()};
    std::fwrite(std::data(headerText), 1, std::size(headerText), stdout);

    std::string buf;
    if (!regionText.has_value()) {
        for (const BamRecord& record : reader.Records()) {
            FormatRecord(record, header, buf);
            std::fwrite(std::data(buf), 1, std::size(buf), stdout);
        }
        return;
    }

    const auto parsedRegion = ParseCramRegion(*regionText);
    if (!parsedRegion) {
        throw std::runtime_error{parsedRegion.error()};
    }

    const std::filesystem::path craiPath =
        indexPath.has_value() ? std::filesystem::path{*indexPath} : DefaultCraiPath(path);
    if (!std::filesystem::exists(craiPath)) {
        throw std::runtime_error{std::format("index file not found: {}", craiPath.string())};
    }

    const CraiIndex index = CraiIndex::FromFile(craiPath);
    const std::int32_t refId = parsedRegion->Unmapped ? -1 : [&]() -> std::int32_t {
        const auto resolved = header.ReferenceId(parsedRegion->RefName);
        if (resolved < 0) {
            throw std::runtime_error{
                std::format("reference '{}' not found in CRAM header", parsedRegion->RefName)};
        }
        return resolved;
    }();

    const auto queried = reader.Query(index, refId, parsedRegion->Beg, parsedRegion->End);
    for (const BamRecord& record : queried) {
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
            })) {
            }
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

int Runner(int argc, char** argv)
{
    // Parse options
    std::int32_t bgzfOpt{-1};
    std::int32_t formatOpt{0};
    const char* inputFile{nullptr};
    const char* referenceFile{nullptr};
    const char* regionText{nullptr};
    const char* indexFile{nullptr};

    for (int i{0}; i < argc; ++i) {
        const std::string_view arg{argv[i]};
        if (arg == "--bgzf-threads") {
            if (i + 1 >= argc) {
                throw std::runtime_error{"missing value for --bgzf-threads"};
            }
            const auto parsed{Tools::ParseInteger<std::int32_t>(argv[++i], "bgzf-threads")};
            if (!parsed) {
                throw std::runtime_error{parsed.error()};
            }
            bgzfOpt = *parsed;
        } else if (arg == "--format-threads") {
            if (i + 1 >= argc) {
                throw std::runtime_error{"missing value for --format-threads"};
            }
            const auto parsed{Tools::ParseInteger<std::int32_t>(argv[++i], "format-threads")};
            if (!parsed) {
                throw std::runtime_error{parsed.error()};
            }
            formatOpt = *parsed;
        } else if (arg == "-j") {
            if (i + 1 >= argc) {
                throw std::runtime_error{"missing value for -j"};
            }
            const auto parsed{Tools::ParseInteger<std::int32_t>(argv[++i], "bgzf-threads")};
            if (!parsed) {
                throw std::runtime_error{parsed.error()};
            }
            bgzfOpt = *parsed;
        } else if (arg == "--reference") {
            if (i + 1 >= argc) {
                throw std::runtime_error{"missing value for --reference"};
            }
            referenceFile = argv[++i];
        } else if (arg == "--region") {
            if (i + 1 >= argc) {
                throw std::runtime_error{"missing value for --region"};
            }
            regionText = argv[++i];
        } else if (arg == "--index") {
            if (i + 1 >= argc) {
                throw std::runtime_error{"missing value for --index"};
            }
            indexFile = argv[++i];
        } else if (arg[0] != '-') {
            inputFile = argv[i];
        } else {
            throw std::runtime_error{std::format("unknown option: {}", arg)};
        }
    }

    if (inputFile == nullptr) {
        std::println(stderr,
                     "Usage: pbsamoa dump [--bgzf-threads N] [--format-threads N] "
                     "[--reference ref.fa] [--region ref:start-end|*] [--index file.crai] "
                     "INPUT.(bam|sam|cram)");
        return EXIT_FAILURE;
    }

    const std::filesystem::path path{inputFile};
    const std::size_t bgzfWorkers{ResolveNumWorkers(bgzfOpt)};

    const bool hasRegion = regionText != nullptr;
    const bool hasIndex = indexFile != nullptr;
    if (path.extension() == ".sam") {
        if (hasRegion || hasIndex) {
            throw std::runtime_error{"--region/--index are only supported for CRAM input"};
        }
        DumpSam(path);
        return EXIT_SUCCESS;
    }
    if (path.extension() == ".cram") {
        if (hasIndex && !hasRegion) {
            throw std::runtime_error{"--index requires --region for CRAM input"};
        }
        DumpCram(path,
                 referenceFile != nullptr ? std::filesystem::path{referenceFile}
                                          : std::filesystem::path{},
                 hasRegion ? std::optional<std::string>{regionText} : std::nullopt,
                 hasIndex ? std::optional<std::filesystem::path>{indexFile} : std::nullopt,
                 bgzfWorkers);
        return EXIT_SUCCESS;
    }

    if (hasRegion || hasIndex) {
        throw std::runtime_error{"--region/--index are only supported for CRAM input"};
    }

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
