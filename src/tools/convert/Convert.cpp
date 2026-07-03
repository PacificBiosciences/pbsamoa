#include "Convert.hpp"

#include "../ParseUtils.hpp"

#include <pbsamoa/PbSamoaLibraryInfo.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamRecordReader.hpp>
#include <pbsamoa/io/CramWriter.hpp>
#include <pbsamoa/io/SamReader.hpp>

#include <pbcopper/cli2/Interface.h>
#include <pbcopper/cli2/Option.h>
#include <pbcopper/cli2/PositionalArgument.h>
#include <pbcopper/cli2/Results.h>

#include <algorithm>
#include <array>
#include <expected>
#include <filesystem>
#include <format>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <cctype>
#include <cstdint>
#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace Convert {
namespace {

// ---- CLIv2 option / positional constants ------------------------------------

const CLI_v2::Option RecordsPerSlice{
    R"({
    "names" : ["records-per-slice"],
    "description" : "Number of records per CRAM slice.",
    "type" : "integer",
    "default" : 0
})"};

// Thread counts use -1 as a sentinel meaning "not user-provided / auto-detect".
const CLI_v2::Option BgzfThreads{
    R"({
    "names" : ["bgzf-threads"],
    "description" : "BGZF decompression threads; 0 = auto.",
    "type" : "integer",
    "default" : -1
})"};

const CLI_v2::Option DecodeThreads{
    R"({
    "names" : ["decode-threads"],
    "description" : "BAM record decode threads; 0 = auto.",
    "type" : "integer",
    "default" : -1
})"};

const CLI_v2::Option ConvertToBamRecord{
    R"({
    "names" : ["convert-to-bam-record"],
    "description" : "Decode BAM records before writing to CRAM."
})"};

const CLI_v2::Option CompressionThreads{
    R"({
    "names" : ["compression-threads"],
    "description" : "CRAM compression threads; 0 = auto.",
    "type" : "integer",
    "default" : -1
})"};

const CLI_v2::Option BlockCompression{
    R"({
    "names" : ["block-compression"],
    "description" : "CRAM block compression method.",
    "type" : "string",
    "choices" : ["raw", "gzip", "bzip2", "lzma", "rans4x8", "rans4x16", "arith", "fqzcomp", "tok"],
    "default" : "rans4x8"
})"};

// Repeatable: --series-compression may appear multiple times.
const CLI_v2::Option SeriesCompression{
    R"({
    "names" : ["series-compression"],
    "description" : "Per-series compression override (SERIES=METHOD); may be repeated.",
    "type" : "string",
    "repeatable" : true
})"};

const CLI_v2::Option WriteCrai{
    R"({
    "names" : ["write-crai"],
    "description" : "Write CRAI index alongside output CRAM."
})"};

const CLI_v2::PositionalArgument Input{
    R"({
    "name" : "input",
    "description" : "Input BAM or SAM file.",
    "type" : "file"
})"};

const CLI_v2::PositionalArgument Output{
    R"({
    "name" : "output",
    "description" : "Output CRAM file.",
    "type" : "file"
})"};

// ---- Domain helpers (unchanged from the original) ---------------------------

std::optional<CramBlockMethod> ParseCompressionMethod(std::string_view arg)
{
    static constexpr std::array METHODS{
        std::pair{"raw", CramBlockMethod::RAW},
        std::pair{"gzip", CramBlockMethod::GZIP},
        std::pair{"bzip2", CramBlockMethod::BZIP2},
        std::pair{"lzma", CramBlockMethod::LZMA},
        std::pair{"rans4x8", CramBlockMethod::RANS4X8},
        std::pair{"rans4x16", CramBlockMethod::RANS4X16},
        std::pair{"arith", CramBlockMethod::ADAPTIVE_ARITH},
        std::pair{"fqzcomp", CramBlockMethod::FQZCOMP},
        std::pair{"tok", CramBlockMethod::NAME_TOKENISER},
    };

    for (const auto& [name, method] : METHODS) {
        if (name == arg) {
            return method;
        }
    }

    const auto parsedId{Tools::ParseInteger<std::int32_t>(arg, "block-compression")};
    if (!parsedId) {
        return std::nullopt;
    }
    if ((*parsedId < 0) || (*parsedId > 8)) {
        return std::nullopt;
    }
    return static_cast<CramBlockMethod>(*parsedId);
}

std::optional<CramDataSeries> ParseDataSeries(std::string_view arg)
{
    if (std::size(arg) != 2) {
        return std::nullopt;
    }

    // Build the enum value directly from the two-character code (matching enum encoding).
    const auto c1{static_cast<char>(std::toupper(static_cast<unsigned char>(arg[0])))};
    const auto c2{static_cast<char>(std::toupper(static_cast<unsigned char>(arg[1])))};
    const auto candidate{static_cast<CramDataSeries>(
        (static_cast<std::uint16_t>(static_cast<std::uint8_t>(c1)) << 8) |
        static_cast<std::uint8_t>(c2))};

    // Validate against known values.
    static constexpr std::array KNOWN = {
        CramDataSeries::BF, CramDataSeries::CF, CramDataSeries::RI, CramDataSeries::RL,
        CramDataSeries::AP, CramDataSeries::RG, CramDataSeries::RN, CramDataSeries::MF,
        CramDataSeries::NS, CramDataSeries::NP, CramDataSeries::TS, CramDataSeries::NF,
        CramDataSeries::TL, CramDataSeries::FN, CramDataSeries::FC, CramDataSeries::FP,
        CramDataSeries::DL, CramDataSeries::BB, CramDataSeries::QQ, CramDataSeries::BS,
        CramDataSeries::IN, CramDataSeries::RS, CramDataSeries::PD, CramDataSeries::HC,
        CramDataSeries::SC, CramDataSeries::MQ, CramDataSeries::BA, CramDataSeries::QS,
    };
    if (std::ranges::find(KNOWN, candidate) != std::ranges::end(KNOWN)) {
        return candidate;
    }
    return std::nullopt;
}

std::expected<std::pair<CramDataSeries, CramBlockMethod>, std::string> ParseSeriesCompression(
    std::string_view arg)
{
    const std::size_t eqPos{arg.find('=')};
    if (eqPos == std::string_view::npos) {
        return std::unexpected{"series-compression must have form SERIES=METHOD"};
    }
    if ((eqPos == 0) || (eqPos + 1 >= std::size(arg))) {
        return std::unexpected{"series-compression must have non-empty SERIES and METHOD"};
    }

    const std::string_view seriesArg{arg.substr(0, eqPos)};
    const std::string_view methodArg{arg.substr(eqPos + 1)};

    const auto parsedSeries{ParseDataSeries(seriesArg)};
    if (!parsedSeries) {
        return std::unexpected{std::format("invalid CRAM data series '{}'", seriesArg)};
    }

    const auto parsedMethod{ParseCompressionMethod(methodArg)};
    if (!parsedMethod) {
        return std::unexpected{std::format("invalid block compression method '{}'", methodArg)};
    }

    return std::pair{*parsedSeries, *parsedMethod};
}

void ConvertBamToCramRaw(const std::filesystem::path& inputPath,
                         const std::filesystem::path& outputPath, const CramWriterConfig& config,
                         std::size_t bgzfWorkers)
{
    BamRawReader reader{inputPath, BamRawReaderConfig{
                                       .BgzfWorkers = bgzfWorkers,
                                   }};
    CramWriter writer{outputPath, reader.Header(), config};
    for (const RawRecord& record : reader.Records()) {
        writer.Write(record);
    }
}

void ConvertBamToCramViaBamRecord(const std::filesystem::path& inputPath,
                                  const std::filesystem::path& outputPath,
                                  const CramWriterConfig& config, std::size_t bgzfWorkers,
                                  std::size_t decodeWorkers)
{
    BamRecordReader reader{inputPath, BamRecordReaderConfig{
                                          .RawReaderConfig =
                                              BamRawReaderConfig{
                                                  .BgzfWorkers = bgzfWorkers,
                                              },
                                          .DecodeWorkers = decodeWorkers,
                                      }};
    CramWriter writer{outputPath, reader.Header(), config};
    for (const BamRecord& record : reader.Records()) {
        writer.Write(record);
    }
}

void ConvertSamToCram(const std::filesystem::path& inputPath,
                      const std::filesystem::path& outputPath, const CramWriterConfig& config)
{
    SamReader reader{inputPath};
    CramWriter writer{outputPath, reader.Header(), config};
    for (const BamRecord& record : reader.Records()) {
        writer.Write(record);
    }
}

}  // namespace

CLI_v2::Interface CreateInterface()
{
    CLI_v2::Interface interface{"pbsamoa convert", "Convert BAM/SAM to CRAM format.",
                                LibraryFormattedVersion()};
    // convert has no --num-threads flag; disable the built-in to avoid shadowing.
    interface.DisableNumThreadsOption();
    interface.AddOptions({RecordsPerSlice, BgzfThreads, DecodeThreads, ConvertToBamRecord,
                          CompressionThreads, BlockCompression, SeriesCompression, WriteCrai});
    interface.AddPositionalArguments({Input, Output});
    return interface;
}

int Runner(const CLI_v2::Results& results)
{
    CramWriterConfig config;

    // --records-per-slice (optional; default 0 = use library default).
    // Use IsUserProvided() to distinguish "not given" (skip) from "user gave 0" (must throw).
    if (results[RecordsPerSlice].IsUserProvided()) {
        const std::int32_t v{results[RecordsPerSlice]};
        if (v < 1) {
            throw std::runtime_error{"records-per-slice must be > 0"};
        }
        config.RecordsPerSlice = v;
    }

    // Thread counts: -1 sentinel = not provided / auto; validate only user-supplied values.
    const std::int32_t bgzfThreadsOpt{results[BgzfThreads]};
    if (results[BgzfThreads].IsUserProvided() && (bgzfThreadsOpt < 0)) {
        throw std::runtime_error{"bgzf-threads must be >= 0"};
    }

    const std::int32_t decodeThreadsOpt{results[DecodeThreads]};
    if (results[DecodeThreads].IsUserProvided() && (decodeThreadsOpt < 0)) {
        throw std::runtime_error{"decode-threads must be >= 0"};
    }

    // Copy-init (not brace-init): Result's bool conversion is ambiguous with brace-init.
    const bool convertToBamRecord = results[ConvertToBamRecord];

    const std::int32_t compressionThreadsOpt{results[CompressionThreads]};
    if (results[CompressionThreads].IsUserProvided() && (compressionThreadsOpt < 0)) {
        throw std::runtime_error{"compression-threads must be >= 0"};
    }

    // --block-compression (default "rans4x8" matches the library's default).
    // choices in the JSON is for --help only; ParseCompressionMethod also accepts 0-8 numerics.
    const std::string blockCompressionStr = results[BlockCompression];
    const auto blockMethod{ParseCompressionMethod(blockCompressionStr)};
    if (!blockMethod) {
        throw std::runtime_error{std::format(
            "invalid block compression method '{}'\n"
            "Valid values: raw,gzip,bzip2,lzma,rans4x8,rans4x16,arith,fqzcomp,tok or numeric 0-8",
            blockCompressionStr)};
    }
    config.BlockCompressionMethod = *blockMethod;

    // --series-compression (repeatable; processed in order of appearance).
    const std::vector<std::string> seriesCompressions =
        results.ToVector<std::string>(SeriesCompression);
    for (const std::string& specStr : seriesCompressions) {
        const auto parsed{ParseSeriesCompression(specStr)};
        if (!parsed) {
            throw std::runtime_error{std::format(
                "{}\n"
                "Valid series include: BF,CF,RI,RL,AP,RG,RN,MF,NS,NP,TS,NF,TL,"
                "FN,FC,FP,MQ,BA,QS,BS,IN,DL,SC,RS,PD,HC,BB,QQ\n"
                "Valid values: raw,gzip,bzip2,lzma,rans4x8,rans4x16,arith,fqzcomp,tok or "
                "numeric 0-8",
                parsed.error())};
        }
        const auto [series, method]{*parsed};
        config.DataSeriesCompressionMethods[series] = method;
    }

    const bool writeCrai = results[WriteCrai];
    config.WriteCrai = writeCrai;

    // CLIv2 does not enforce required positional count — guard before indexing.
    const std::vector<std::string>& positional{results.PositionalArguments()};
    if (positional.size() != 2) {
        throw std::runtime_error{
            "convert requires exactly two arguments: <input.bam|sam> <output.cram>"};
    }

    const std::filesystem::path inputPath{positional[0]};
    const std::filesystem::path outputPath{positional[1]};
    constexpr std::int32_t MAX_IO_WORKERS{16};
    constexpr std::int32_t MAX_COMPRESSION_WORKERS{8};
    const std::size_t bgzfWorkers{Tools::ResolveNumWorkers(bgzfThreadsOpt, MAX_IO_WORKERS)};
    const std::size_t decodeWorkers{Tools::ResolveNumWorkers(decodeThreadsOpt, MAX_IO_WORKERS)};
    config.CompressionWorkers =
        Tools::ResolveNumWorkers(compressionThreadsOpt, MAX_COMPRESSION_WORKERS);

    if (outputPath.extension() != ".cram") {
        throw std::runtime_error{"output must use .cram extension"};
    }

    if (inputPath.extension() == ".bam") {
        if (convertToBamRecord || (decodeThreadsOpt >= 0)) {
            ConvertBamToCramViaBamRecord(inputPath, outputPath, config, bgzfWorkers, decodeWorkers);
        } else {
            ConvertBamToCramRaw(inputPath, outputPath, config, bgzfWorkers);
        }
        return EXIT_SUCCESS;
    }
    if (inputPath.extension() == ".sam") {
        ConvertSamToCram(inputPath, outputPath, config);
        return EXIT_SUCCESS;
    }

    throw std::runtime_error{"input must use .bam or .sam extension"};
}

}  // namespace Convert
}  // namespace Samoa
}  // namespace PacBio
