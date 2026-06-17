#include "Dump.hpp"
#include "../../BinaryUtils.hpp"
#include "../../CramInternal.hpp"
#include "../../SamFieldUtils.hpp"
#include "../../WriterUtils.hpp"
#include "../CliUtils.hpp"
#include "../MetricUtils.hpp"
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

#include <pbcopper/parallel/ThreadPool.h>

#include <algorithm>
#include <array>
#include <charconv>
#include <chrono>
#include <exception>
#include <expected>
#include <filesystem>
#include <format>
#include <functional>
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

std::expected<CramRegion, std::string> ParseCramRegion(std::string_view text)
{
    if (text == "*") {
        return CramRegion{.Unmapped = true};
    }

    const auto parsedRegion{Tools::ParseRegion(text, "ref:start-end or '*'")};
    if (!parsedRegion) {
        return std::unexpected{parsedRegion.error()};
    }

    return CramRegion{
        .RefName = parsedRegion->RefName, .Beg = parsedRegion->Beg, .End = parsedRegion->End};
}

void RequestStop(std::jthread& thread)
{
    if (thread.joinable()) {
        thread.request_stop();
    }
}

void JoinThread(std::jthread& thread)
{
    if (thread.joinable()) {
        thread.join();
    }
}

void AppendInt(std::string& out, std::int64_t v)
{
    std::array<char, 24> buf{};
    const auto result{std::to_chars(std::data(buf), std::data(buf) + std::size(buf), v)};
    out.append(std::data(buf), result.ptr);
}

void AppendReferenceName(std::string& out, const SamHeader& header, std::int32_t refId)
{
    const std::string_view sentinel{RnameSentinel(refId)};
    if (!sentinel.empty()) {
        out.append(sentinel);
        return;
    }
    out.append(header.ReferenceName(refId));
}

void AppendNextReferenceName(std::string& out, const SamHeader& header, std::int32_t refId,
                             std::int32_t nextRefId)
{
    const std::string_view sentinel{RnextSentinel(refId, nextRefId)};
    if (!sentinel.empty()) {
        out.append(sentinel);
        return;
    }
    out.append(header.ReferenceName(nextRefId));
}

void AppendQualities(std::string& out, std::span<const std::uint8_t> qual)
{
    if (IsQualityUnavailable(qual)) {
        out += '*';
        return;
    }

    const std::size_t startPos{std::size(out)};
    out.resize(startPos + std::size(qual));
    for (std::size_t qi{0}; qi < std::size(qual); ++qi) {
        out[startPos + qi] = static_cast<char>(qual[qi] + 33);
    }
}

void AppendSequenceField(std::string& out, const RawRecord& view)
{
    if (view.SeqLength() == 0) {
        out += '*';
        return;
    }
    view.Seq().WriteTo(out);
}

void AppendSequenceField(std::string& out, const BamRecord& record)
{
    if (std::empty(record.Sequence())) {
        out += '*';
        return;
    }
    out.append(record.Sequence());
}

void AppendTagFields(std::string& out, const RawRecord& view)
{
    SerializeRawTagsToSam(view.AuxData(), out);
}

void AppendTagFields(std::string& out, const BamRecord& record)
{
    for (const auto& [key, value] : record.Tags().Entries()) {
        out += '\t';
        out.append(SerializeTagToSam(key, value));
    }
}

void AppendCommonSamFields(std::string& out, const SamHeader& header, std::string_view name,
                           std::uint16_t flag, std::int32_t refId, std::int32_t pos,
                           std::uint8_t mapQ, CigarView cigar, std::int32_t nextRefId,
                           std::int32_t nextPos, std::int32_t tlen)
{
    out.append(name);
    out += '\t';
    AppendInt(out, flag);
    out += '\t';
    AppendReferenceName(out, header, refId);
    out += '\t';
    AppendInt(out, OneBasedPositionOrZero(pos));
    out += '\t';
    AppendInt(out, mapQ);
    out += '\t';
    WriteCigarTo(cigar, out);
    out += '\t';
    AppendNextReferenceName(out, header, refId, nextRefId);
    out += '\t';
    AppendInt(out, OneBasedPositionOrZero(nextPos));
    out += '\t';
    AppendInt(out, tlen);
    out += '\t';
}

template <typename Record>
void FormatRecordImpl(const Record& record, const SamHeader& header, std::string& buf,
                      std::string_view name, std::uint16_t flag, std::int32_t refId,
                      std::int32_t pos, std::uint8_t mapQ, CigarView cigar, std::int32_t nextRefId,
                      std::int32_t nextPos, std::int32_t tlen,
                      std::span<const std::uint8_t> qualities)
{
    buf.clear();
    AppendCommonSamFields(buf, header, name, flag, refId, pos, mapQ, cigar, nextRefId, nextPos,
                          tlen);
    AppendSequenceField(buf, record);
    buf += '\t';
    AppendQualities(buf, qualities);
    AppendTagFields(buf, record);
    buf += '\n';
}

void WriteHeader(const SamHeader& header, HeaderMode headerMode)
{
    if (headerMode == HeaderMode::NoHeader) {
        return;
    }

    const std::string headerText{header.ToText()};
    std::fwrite(std::data(headerText), 1, std::size(headerText), stdout);
}

void WriteStdout(const std::string& output)
{
    std::fwrite(std::data(output), 1, std::size(output), stdout);
}

void FormatRecord(const RawRecord& view, const SamHeader& header, std::string& buf);
void FormatRecord(const BamRecord& record, const SamHeader& header, std::string& buf);
void PrintMetricsLine(const BgzfMetrics& m, const BgzfMetrics& prev, double elapsedSec);

std::string FormatBatchRecord(std::shared_ptr<const RawRecordBatch> batch, const SamHeader& header,
                              std::size_t recordIdx)
{
    std::string output;
    const RawRecord view{batch->RecordData(recordIdx)};
    FormatRecord(view, header, output);
    return output;
}

template <typename Records>
void WriteFormattedRecords(Records&& records, const SamHeader& header, std::string& buffer)
{
    for (const auto& record : records) {
        FormatRecord(record, header, buffer);
        std::fwrite(std::data(buffer), 1, std::size(buffer), stdout);
    }
}

void RunMetricsLoop(std::stop_token stopToken, const BamRawReader& reader)
{
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
            std::chrono::duration_cast<std::chrono::duration<double>>(now - lastTime).count()};
        const BgzfMetrics current{reader.GetMetrics()};
        PrintMetricsLine(current, prev, elapsed);
        prev = current;
        lastTime = now;
    }
}

void ConsumeFormattedOutput(std::stop_token, PacBio::Parallel::ThreadPool<std::string>& formatPool,
                            std::exception_ptr& consumerException)
{
    try {
        while (formatPool.ConsumeWith(WriteStdout)) {
        }
    } catch (...) {
        consumerException = std::current_exception();
    }
}

void FormatRecord(const RawRecord& view, const SamHeader& header, std::string& buf)
{
    const std::int32_t refId{view.RefId()};
    FormatRecordImpl(view, header, buf, view.Name(), view.Flag(), refId, view.Pos(), view.MapQ(),
                     view.CigarOps(), view.NextRefId(), view.NextPos(), view.Tlen(), view.Qual());
}

void FormatRecord(const BamRecord& record, const SamHeader& header, std::string& buf)
{
    const std::int32_t refId{record.RefId()};
    FormatRecordImpl(record, header, buf, record.Name(), record.Flag(), refId, record.Pos(),
                     record.MapQ(), record.Cigar(), record.NextRefId(), record.NextPos(),
                     record.Tlen(), record.Qualities());
}

void PrintMetricsLine(const BgzfMetrics& m, const BgzfMetrics& prev, double elapsedSec)
{
    const double mbRead{Tools::ToMiB(m.BytesRead)};
    const double mbDecomp{Tools::ToMiB(m.BytesDecompressed)};
    const double deltaRecords = m.RecordsConsumed - prev.RecordsConsumed;
    const double recPerSec{Tools::RateOrZero(deltaRecords, elapsedSec)};
    const double deltaMbDecomp{Tools::ToMiB(m.BytesDecompressed - prev.BytesDecompressed)};
    const double mbPerSec{Tools::RateOrZero(deltaMbDecomp, elapsedSec)};
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
    const double mbRead{Tools::ToMiB(m.BytesRead)};
    const double mbDecomp{Tools::ToMiB(m.BytesDecompressed)};
    const double ioMs{Tools::ToMs(m.IoReadNs)};
    const double decompMs{Tools::ToMs(m.DecompressNs)};
    const double parseMs{Tools::ToMs(m.RecordParseNs)};

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
    WriteFormattedRecords(reader.Records(), header, buf);
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
    if (!regionText) {
        WriteFormattedRecords(reader.Records(), header, buf);
        return;
    }

    const auto parsedRegion = ParseCramRegion(*regionText);
    if (!parsedRegion) {
        throw std::runtime_error{parsedRegion.error()};
    }

    std::filesystem::path craiPath{DefaultCraiPath(path)};
    if (indexPath) {
        craiPath = *indexPath;
    }
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
    WriteFormattedRecords(queried, header, buf);
}

void DumpBam(const std::filesystem::path& path, std::size_t numWorkers,
             std::size_t numFormatThreads, HeaderMode headerMode)
{
    const char* const env{std::getenv("PBSAMOA_METRICS")};
    const bool showMetrics{env && (std::string_view{env} == "1")};

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
        metricsThread = std::jthread{RunMetricsLoop, std::cref(reader)};
    }

    PacBio::Parallel::ThreadPool<std::string> formatPool{
        PacBio::Parallel::ThreadPool<std::string>::Config{.NumThreads = numFormatThreads,
                                                          .QueueMultiplier = 5}};
    std::exception_ptr consumerException;
    std::jthread consumerThread{ConsumeFormattedOutput, std::ref(formatPool),
                                std::ref(consumerException)};

    try {
        // Double-buffer: read next batch while formatting current
        auto currentBatch{reader.ReadBatch()};

        while (currentBatch) {
            // Start reading next batch concurrently
            auto nextBatchFuture{
                std::async(std::launch::async, &BamRawReader::ReadBatch, &reader, ByteLimit{})};

            const auto batch{std::make_shared<RawRecordBatch>(std::move(*currentBatch))};
            for (std::size_t recordIdx{0}; recordIdx < batch->RecordCount(); ++recordIdx) {
                formatPool.Submit(FormatBatchRecord, batch, std::cref(header), recordIdx);
            }

            currentBatch = nextBatchFuture.get();
        }

        formatPool.Finalize();
        JoinThread(consumerThread);
        if (showMetrics) {
            RequestStop(metricsThread);
            JoinThread(metricsThread);
            PrintMetricsSummary(reader.GetMetrics());
        }
    } catch (...) {
        RequestStop(metricsThread);
        try {
            formatPool.Finalize();
        } catch (...) {
        }
        JoinThread(consumerThread);
        JoinThread(metricsThread);
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

    for (int i{0}; i < argc; ++i) {
        const std::string_view arg{argv[i]};
        if (arg == "--bgzf-threads") {
            bgzfOpt = Tools::ParseIntOption<std::int32_t>(argc, argv, i, "--bgzf-threads",
                                                          "bgzf-threads");
        } else if (arg == "--format-threads") {
            formatOpt = Tools::ParseIntOption<std::int32_t>(argc, argv, i, "--format-threads",
                                                            "format-threads");
        } else if (arg == "-j") {
            bgzfOpt = Tools::ParseIntOption<std::int32_t>(argc, argv, i, "-j", "bgzf-threads");
        } else if (arg == "--reference") {
            referenceFile = Tools::RequireOptionValue(argc, argv, i, "--reference");
        } else if (arg == "--region") {
            regionText = Tools::RequireOptionValue(argc, argv, i, "--region");
        } else if (arg == "--index") {
            indexFile = Tools::RequireOptionValue(argc, argv, i, "--index");
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

    if (!inputFile) {
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

    HeaderMode headerMode{HeaderMode::Full};
    if (noHeader) {
        headerMode = HeaderMode::NoHeader;
    } else if (headerOnly) {
        headerMode = HeaderMode::HeaderOnly;
    }

    const std::filesystem::path path{inputFile};
    const std::size_t bgzfWorkers{Tools::ResolveNumWorkers(bgzfOpt, /*explicitCap=*/10)};
    std::setvbuf(stdout, std::data(stdoutBuf), _IOFBF, std::size(stdoutBuf));

    if (path.extension() == ".sam") {
        if (regionText || indexFile) {
            throw std::runtime_error{"--region/--index are only supported for CRAM input"};
        }
        DumpSam(path, headerMode);
        return EXIT_SUCCESS;
    }
    if (path.extension() == ".cram") {
        if (indexFile && !regionText) {
            throw std::runtime_error{"--index requires --region for CRAM input"};
        }
        std::filesystem::path referencePath{};
        if (referenceFile) {
            referencePath = referenceFile;
        }

        std::optional<std::string> region{};
        if (regionText) {
            region.emplace(regionText);
        }

        std::optional<std::filesystem::path> indexPath{};
        if (indexFile) {
            indexPath.emplace(indexFile);
        }

        DumpCram(path, referencePath, region, indexPath, bgzfWorkers, headerMode);
        return EXIT_SUCCESS;
    }

    if (regionText || indexFile) {
        throw std::runtime_error{"--region/--index are only supported for CRAM input"};
    }

    std::size_t formatWorkers{
        static_cast<std::size_t>(std::ranges::max(std::thread::hardware_concurrency(), 4U))};
    if (formatOpt > 0) {
        formatWorkers = static_cast<std::size_t>(formatOpt);
    }
    DumpBam(path, bgzfWorkers, formatWorkers, headerMode);
    return EXIT_SUCCESS;
}

}  // namespace Dump
}  // namespace Samoa
}  // namespace PacBio
