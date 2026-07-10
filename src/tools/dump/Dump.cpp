#include "Dump.hpp"
#include "../../BinaryUtils.hpp"
#include "../../CramInternal.hpp"
#include "../../SamFieldUtils.hpp"
#include "../../WriterUtils.hpp"
#include "../MetricUtils.hpp"
#include "../ParseUtils.hpp"

#include <pbsamoa/PbSamoaLibraryInfo.hpp>
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

#include <pbcopper/cli2/Interface.h>
#include <pbcopper/cli2/Option.h>
#include <pbcopper/cli2/PositionalArgument.h>
#include <pbcopper/cli2/Results.h>
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
#include <vector>

#include <cstdint>
#include <cstdio>
#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace Dump {
namespace {

// ---- CLIv2 option / positional constants ----

// --bgzf-threads is long-only: `-j` belongs to the built-in --num-threads (disabled
// below via DisableNumThreadsOption()), whose auto-resolution differs from dump's
// -1 = auto (capped at 10).
const CLI_v2::Option BgzfThreads{
    R"({
    "names" : ["bgzf-threads"],
    "description" : "BGZF decompression threads; -1 = auto (capped at 10).",
    "type" : "integer",
    "default" : -1
})"};

const CLI_v2::Option FormatThreads{
    R"({
    "names" : ["format-threads"],
    "description" : "SAM formatting threads; 0 = auto (max of hardware concurrency and 4). BAM only.",
    "type" : "integer",
    "default" : 0
})"};

const CLI_v2::Option Reference{
    R"({
    "names" : ["reference"],
    "description" : "Reference FASTA path (required for reference-based CRAM).",
    "type" : "file"
})"};

const CLI_v2::Option Region{
    R"({
    "names" : ["region"],
    "description" : "Genomic region to query (ref:start-end or * for unmapped). CRAM only.",
    "type" : "string"
})"};

const CLI_v2::Option Index{
    R"({
    "names" : ["index"],
    "description" : "Index file path (.crai). CRAM only.",
    "type" : "file"
})"};

const CLI_v2::Option NoHeader{
    R"({
    "names" : ["no-header"],
    "description" : "Suppress the SAM header."
})"};

const CLI_v2::Option HeaderOnly{
    R"({
    "names" : ["header-only"],
    "description" : "Output only the SAM header."
})"};

const CLI_v2::PositionalArgument InputFile{
    R"({
    "name" : "input",
    "description" : "Input BAM/SAM/CRAM file.",
    "type" : "file"
})"};

// ---- Internal types and helpers (unchanged from legacy) ----

enum class HeaderMode
{
    Full,
    Suppress,
    OnlyHeader,
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
    if (headerMode == HeaderMode::Suppress) {
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
    if (headerMode == HeaderMode::OnlyHeader) {
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
    if (headerMode == HeaderMode::OnlyHeader) {
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
    if (headerMode == HeaderMode::OnlyHeader) {
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

CLI_v2::Interface CreateInterface()
{
    CLI_v2::Interface interface{"pbsamoa dump",
                                "Convert a BAM/SAM/CRAM file to SAM text on stdout.",
                                LibraryFormattedVersion()};
    // dump defines its own `-j` (bgzf threads); disable the built-in `-j`/`--num-threads`
    // so it does not shadow or conflict with dump's custom `-j`.
    interface.DisableNumThreadsOption();
    interface.AddOptions(
        {BgzfThreads, FormatThreads, Reference, Region, Index, NoHeader, HeaderOnly});
    interface.AddPositionalArguments({InputFile});
    return interface;
}

int Runner(const CLI_v2::Results& results)
{
    static std::array<char, 1 << 20> stdoutBuf;

    // CLIv2 does not enforce required positional count; guard before indexing.
    const std::vector<std::string>& positional{results.PositionalArguments()};
    if (positional.size() != 1) {
        throw std::runtime_error{"dump requires exactly one argument: INPUT"};
    }

    // Numeric reads: brace-init exact-width type (never read into std::size_t directly).
    const std::int32_t bgzfOpt{results[BgzfThreads]};
    const std::int32_t formatOpt{results[FormatThreads]};
    if (bgzfOpt < -1) {
        throw std::runtime_error{"--bgzf-threads must be >= -1"};
    }
    if (formatOpt < 0) {
        throw std::runtime_error{"--format-threads must be >= 0"};
    }

    // String reads: copy-init (brace-init is ambiguous due to Result's char-type conversion).
    const std::string referenceStr = results[Reference];
    const std::string regionStr = results[Region];
    const std::string indexStr = results[Index];

    // Bool flag reads: copy-init.
    const bool noHeader = results[NoHeader];
    const bool headerOnly = results[HeaderOnly];

    if (noHeader && headerOnly) {
        throw std::runtime_error{"--no-header and --header-only are mutually exclusive"};
    }

    HeaderMode headerMode{HeaderMode::Full};
    if (noHeader) {
        headerMode = HeaderMode::Suppress;
    } else if (headerOnly) {
        headerMode = HeaderMode::OnlyHeader;
    }

    const std::filesystem::path path{positional[0]};
    const std::size_t bgzfWorkers{Tools::ResolveNumWorkers(bgzfOpt, /*explicitCap=*/10)};
    std::setvbuf(stdout, std::data(stdoutBuf), _IOFBF, std::size(stdoutBuf));

    const bool hasRegion{!regionStr.empty()};
    const bool hasIndex{!indexStr.empty()};

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
        std::filesystem::path referencePath{};
        if (!referenceStr.empty()) {
            referencePath = referenceStr;
        }

        std::optional<std::string> region{};
        if (hasRegion) {
            region.emplace(regionStr);
        }

        std::optional<std::filesystem::path> indexPath{};
        if (hasIndex) {
            indexPath.emplace(indexStr);
        }

        DumpCram(path, referencePath, region, indexPath, bgzfWorkers, headerMode);
        return EXIT_SUCCESS;
    }

    if (hasRegion || hasIndex) {
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
