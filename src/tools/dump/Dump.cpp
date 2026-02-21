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

#include <parallel/ThreadPool.h>

#include <algorithm>
#include <array>
#include <charconv>
#include <chrono>
#include <exception>
#include <expected>
#include <filesystem>
#include <format>
#include <future>
#include <memory>
#include <optional>
#include <print>
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

enum class HeaderMode
{
    Full,
    NoHeader,
    HeaderOnly,
};

struct CramRegion
{
    bool Unmapped{false};
    std::string RefName;
    std::int32_t Beg{0};  // 0-based inclusive
    std::int32_t End{0};  // 0-based exclusive
};

HeaderMode ResolveHeaderMode(bool noHeader, bool headerOnly)
{
    if (noHeader) {
        return HeaderMode::NoHeader;
    }
    if (headerOnly) {
        return HeaderMode::HeaderOnly;
    }
    return HeaderMode::Full;
}

std::int32_t OneBasedPosOrZero(std::int32_t pos)
{
    if (pos < 0) {
        return 0;
    }
    return pos + 1;
}

std::expected<CramRegion, std::string> ParseCramRegion(std::string_view text)
{
    constexpr std::string_view FORMAT_ERROR =
        "invalid region format '{}' (expected ref:start-end or '*')";

    if (text == "*") {
        return CramRegion{.Unmapped = true, .Beg = 0, .End = 0};
    }

    const std::size_t colonPos{text.find(':')};
    if (colonPos == std::string_view::npos) {
        return std::unexpected{std::format(FORMAT_ERROR, text)};
    }

    const std::string refName{text.substr(0, colonPos)};
    const std::string_view rest{text.substr(colonPos + 1)};
    const std::size_t dashPos{rest.find('-')};
    if (dashPos == std::string_view::npos) {
        return std::unexpected{std::format(FORMAT_ERROR, text)};
    }

    const auto startResult{
        Tools::ParseInteger<std::int32_t>(rest.substr(0, dashPos), "region start")};
    if (!startResult) {
        return std::unexpected{startResult.error()};
    }

    const auto endResult{Tools::ParseInteger<std::int32_t>(rest.substr(dashPos + 1), "region end")};
    if (!endResult) {
        return std::unexpected{endResult.error()};
    }

    const std::int32_t start{*startResult};
    const std::int32_t end{*endResult};
    if ((start < 1) || (end < start)) {
        return std::unexpected{"region must satisfy start >= 1 and end >= start"};
    }

    return CramRegion{
        .Unmapped = false,
        .RefName = refName,
        .Beg = start - 1,
        .End = end,
    };
}

void AppendInt(std::string& out, std::int64_t v)
{
    std::array<char, 24> buf{};
    const auto [ptr, ec]{std::to_chars(std::data(buf), std::data(buf) + std::size(buf), v)};
    out.append(std::data(buf), ptr);
}

void AppendReferenceName(std::string& out, const SamHeader& header, std::int32_t refId)
{
    if (refId < 0) {
        out += '*';
        return;
    }
    out.append(header.ReferenceName(refId));
}

void AppendNextReferenceName(std::string& out, const SamHeader& header, std::int32_t refId,
                             std::int32_t nextRefId)
{
    if (nextRefId < 0) {
        out += '*';
        return;
    }
    if (nextRefId == refId) {
        out += '=';
        return;
    }
    out.append(header.ReferenceName(nextRefId));
}

void AppendQualities(std::string& out, std::span<const std::uint8_t> qual)
{
    const bool qualUnavailable{std::empty(qual) ||
                               std::ranges::all_of(qual, [](std::uint8_t q) { return q == 0xFF; })};
    if (qualUnavailable) {
        out += '*';
        return;
    }

    const std::size_t startPos{std::size(out)};
    out.resize(startPos + std::size(qual));
    for (std::size_t qi{0}; qi < std::size(qual); ++qi) {
        out[startPos + qi] = static_cast<char>(qual[qi] + 33);
    }
}

void WriteHeader(const SamHeader& header, HeaderMode headerMode)
{
    if (headerMode == HeaderMode::NoHeader) {
        return;
    }

    const std::string headerText{header.ToText()};
    std::fwrite(std::data(headerText), 1, std::size(headerText), stdout);
}

void FormatRecord(const RawRecord& view, const SamHeader& header, std::string& buf)
{
    buf.clear();

    buf.append(view.Name());
    buf += '\t';

    AppendInt(buf, view.Flag());
    buf += '\t';

    const std::int32_t refId{view.RefId()};
    AppendReferenceName(buf, header, refId);
    buf += '\t';

    const std::int32_t pos{view.Pos()};
    AppendInt(buf, OneBasedPosOrZero(pos));
    buf += '\t';

    AppendInt(buf, view.MapQ());
    buf += '\t';

    WriteCigarTo(view.CigarOps(), buf);
    buf += '\t';

    const std::int32_t nextRefId{view.NextRefId()};
    AppendNextReferenceName(buf, header, refId, nextRefId);
    buf += '\t';

    const std::int32_t nextPos{view.NextPos()};
    AppendInt(buf, OneBasedPosOrZero(nextPos));
    buf += '\t';

    AppendInt(buf, view.Tlen());
    buf += '\t';

    if (view.SeqLength() == 0) {
        buf += '*';
    } else {
        view.Seq().WriteTo(buf);
    }
    buf += '\t';

    AppendQualities(buf, view.Qual());

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
    AppendReferenceName(buf, header, refId);
    buf += '\t';

    const std::int32_t pos{record.Pos()};
    AppendInt(buf, OneBasedPosOrZero(pos));
    buf += '\t';

    AppendInt(buf, record.MapQ());
    buf += '\t';

    WriteCigarTo(record.Cigar(), buf);
    buf += '\t';

    const std::int32_t nextRefId{record.NextRefId()};
    AppendNextReferenceName(buf, header, refId, nextRefId);
    buf += '\t';

    const std::int32_t nextPos{record.NextPos()};
    AppendInt(buf, OneBasedPosOrZero(nextPos));
    buf += '\t';

    AppendInt(buf, record.Tlen());
    buf += '\t';

    if (std::empty(record.Sequence())) {
        buf += '*';
    } else {
        buf.append(record.Sequence());
    }
    buf += '\t';

    AppendQualities(buf, record.Qualities());

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
                 mbRead, mbDecomp, recPerSec, mbPerSec, m.Pool.ActiveWorkers, m.Pool.QueueDepth,
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
    std::println(stderr, "Pool peak:     queue {}, workers {}", m.Pool.PeakQueueDepth,
                 m.Pool.PeakActiveWorkers);
    std::println(stderr, "Result peak:   queue {}", m.Pool.PeakResultQueueDepth);
    std::println(stderr, "Stalls:        IO={}, consumer={}, reader={}", m.IoStalls,
                 m.ConsumerStalls, m.ReaderStalls);
    std::println(stderr, "Timing:        IO {:.1f}ms, decompress {:.1f}ms, parse {:.1f}ms", ioMs,
                 decompMs, parseMs);
}

void DumpSam(const std::filesystem::path& path, HeaderMode headerMode)
{
    SamReader reader{path};
    const SamHeader& header{reader.Header()};

    WriteHeader(header, headerMode);
    if (headerMode == HeaderMode::HeaderOnly) {
        return;
    }

    std::string buf;
    for (const BamRecord& record : reader.Records()) {
        FormatRecord(record, header, buf);
        std::fwrite(std::data(buf), 1, std::size(buf), stdout);
    }
}

void DumpCram(const std::filesystem::path& path, const std::filesystem::path& referencePath,
              const std::optional<std::string>& regionText,
              const std::optional<std::filesystem::path>& indexPath, std::size_t numWorkers,
              HeaderMode headerMode)
{
    CramReaderConfig config;
    config.ReferencePath = referencePath;
    config.DecompressionWorkers = numWorkers;
    CramReader reader{path, config};
    const SamHeader& header{reader.Header()};

    WriteHeader(header, headerMode);
    if (headerMode == HeaderMode::HeaderOnly) {
        return;
    }

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
    std::int32_t refId{-1};
    if (!parsedRegion->Unmapped) {
        refId = header.ReferenceId(parsedRegion->RefName);
        if (refId < 0) {
            throw std::runtime_error{
                std::format("reference '{}' not found in CRAM header", parsedRegion->RefName)};
        }
    }

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
             std::size_t numFormatThreads, HeaderMode headerMode)
{
    const bool showMetrics{MetricsEnabled()};

    // Pipeline handles BGZF I/O + decompression with its own threads
    BamRawReader reader{path, BamRawReaderConfig{.BgzfWorkers = numWorkers}};

    const auto& header{reader.Header()};

    WriteHeader(header, headerMode);
    if (headerMode == HeaderMode::HeaderOnly) {
        return;
    }

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
    static std::array<char, 1 << 20> stdoutBuf;

    // Parse options
    std::int32_t bgzfOpt{-1};
    std::int32_t formatOpt{0};
    const char* inputFile{nullptr};
    const char* referenceFile{nullptr};
    const char* regionText{nullptr};
    const char* indexFile{nullptr};
    bool noHeader{false};
    bool headerOnly{false};
    const auto requireValue = [&](int currentIndex, std::string_view option) -> const char* {
        if (currentIndex + 1 >= argc) {
            throw std::runtime_error{std::format("missing value for {}", option)};
        }
        return argv[++currentIndex];
    };
    const auto parseIntOption = [&](int currentIndex, std::string_view option,
                                    std::string_view parseName) -> std::int32_t {
        return Tools::ParseIntegerOrThrow<std::int32_t>(requireValue(currentIndex, option),
                                                        parseName);
    };

    for (int i{0}; i < argc; ++i) {
        const std::string_view arg{argv[i]};
        if (arg == "--bgzf-threads") {
            bgzfOpt = parseIntOption(i, "--bgzf-threads", "bgzf-threads");
            ++i;
        } else if (arg == "--format-threads") {
            formatOpt = parseIntOption(i, "--format-threads", "format-threads");
            ++i;
        } else if (arg == "-j") {
            bgzfOpt = parseIntOption(i, "-j", "bgzf-threads");
            ++i;
        } else if (arg == "--reference") {
            referenceFile = requireValue(i, "--reference");
            ++i;
        } else if (arg == "--region") {
            regionText = requireValue(i, "--region");
            ++i;
        } else if (arg == "--index") {
            indexFile = requireValue(i, "--index");
            ++i;
        } else if (arg == "--no-header") {
            noHeader = true;
        } else if (arg == "--header-only") {
            headerOnly = true;
        } else if (!std::empty(arg) && arg[0] != '-') {
            inputFile = argv[i];
        } else {
            throw std::runtime_error{std::format("unknown option: {}", arg)};
        }
    }

    if (inputFile == nullptr) {
        std::println(stderr,
                     "Usage: pbsamoa dump [--bgzf-threads N] [--format-threads N] "
                     "[--no-header|--header-only] [--reference ref.fa] "
                     "[--region ref:start-end|*] [--index file.crai] "
                     "INPUT.(bam|sam|cram)");
        return EXIT_FAILURE;
    }

    if (noHeader && headerOnly) {
        throw std::runtime_error{"--no-header and --header-only are mutually exclusive"};
    }

    const HeaderMode headerMode{ResolveHeaderMode(noHeader, headerOnly)};

    const std::filesystem::path path{inputFile};
    const std::size_t bgzfWorkers{ResolveNumWorkers(bgzfOpt)};
    std::setvbuf(stdout, std::data(stdoutBuf), _IOFBF, std::size(stdoutBuf));

    const bool hasRegion = regionText != nullptr;
    const bool hasIndex = indexFile != nullptr;
    if (path.extension() == ".sam") {
        if (hasRegion || hasIndex) {
            throw std::runtime_error{"--region/--index are only supported for CRAM input"};
        }
        DumpSam(path, headerMode);
        return EXIT_SUCCESS;
    }
    if (path.extension() == ".cram") {
        if (hasIndex && !hasRegion) {
            throw std::runtime_error{"--index requires --region for CRAM input"};
        }
        const std::filesystem::path referencePath{referenceFile != nullptr
                                                      ? std::filesystem::path{referenceFile}
                                                      : std::filesystem::path{}};
        const std::optional<std::string> region{hasRegion ? std::optional<std::string>{regionText}
                                                          : std::nullopt};
        const std::optional<std::filesystem::path> indexPath{
            hasIndex ? std::optional<std::filesystem::path>{indexFile} : std::nullopt};
        DumpCram(path, referencePath, region, indexPath, bgzfWorkers, headerMode);
        return EXIT_SUCCESS;
    }

    if (hasRegion || hasIndex) {
        throw std::runtime_error{"--region/--index are only supported for CRAM input"};
    }

    std::size_t formatThreads{};
    if (formatOpt > 0) {
        formatThreads = static_cast<std::size_t>(formatOpt);
    } else {
        formatThreads =
            static_cast<std::size_t>(std::ranges::max(std::thread::hardware_concurrency(), 4U));
    }

    DumpBam(path, bgzfWorkers, formatThreads, headerMode);
    return EXIT_SUCCESS;
}

}  // namespace Dump
}  // namespace Samoa
}  // namespace PacBio
