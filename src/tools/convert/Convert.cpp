#include "Convert.hpp"

#include "../ParseUtils.hpp"

#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamRecordReader.hpp>
#include <pbsamoa/io/CramWriter.hpp>
#include <pbsamoa/io/SamReader.hpp>

#include <algorithm>
#include <array>
#include <filesystem>
#include <format>
#include <optional>
#include <print>
#include <ranges>
#include <string>
#include <string_view>

#include <cctype>
#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace Convert {
namespace {

void PrintUsage()
{
    std::println(stderr,
                 "Usage: pbsamoa convert [--records-per-slice N] "
                 "[--bgzf-threads N] [--convert-to-bam-record] "
                 "[--decode-threads N] "
                 "[--compression-threads N] "
                 "[--block-compression METHOD] "
                 "[--series-compression SERIES=METHOD]... [--write-crai] "
                 "INPUT.(bam|sam) OUTPUT.cram");
}

std::optional<CramBlockMethod> ParseCompressionMethod(const std::string_view arg)
{
    if (arg == "raw") {
        return CramBlockMethod::RAW;
    }
    if (arg == "gzip") {
        return CramBlockMethod::GZIP;
    }
    if (arg == "bzip2") {
        return CramBlockMethod::BZIP2;
    }
    if (arg == "lzma") {
        return CramBlockMethod::LZMA;
    }
    if (arg == "rans4x8") {
        return CramBlockMethod::RANS4X8;
    }
    if (arg == "rans4x16") {
        return CramBlockMethod::RANS4X16;
    }
    if (arg == "arith") {
        return CramBlockMethod::ADAPTIVE_ARITH;
    }
    if (arg == "fqzcomp") {
        return CramBlockMethod::FQZCOMP;
    }
    if (arg == "tok") {
        return CramBlockMethod::NAME_TOKENISER;
    }

    const auto parsedId = Tools::ParseInteger<std::int32_t>(arg, "block-compression");
    if (parsedId && *parsedId >= 0 && *parsedId <= 8) {
        return static_cast<CramBlockMethod>(*parsedId);
    }
    return std::nullopt;
}

std::optional<CramDataSeries> ParseDataSeries(const std::string_view arg)
{
    if (std::size(arg) != 2) {
        return std::nullopt;
    }

    // Build the enum value directly from the two-character code (matching enum encoding).
    const auto c1 = static_cast<char>(std::toupper(static_cast<unsigned char>(arg[0])));
    const auto c2 = static_cast<char>(std::toupper(static_cast<unsigned char>(arg[1])));
    const auto candidate = static_cast<CramDataSeries>(
        (static_cast<std::uint16_t>(static_cast<std::uint8_t>(c1)) << 8) |
        static_cast<std::uint8_t>(c2));

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

std::optional<std::pair<CramDataSeries, CramBlockMethod>> ParseSeriesCompression(
    const std::string_view arg, std::string& error)
{
    const auto eqPos = arg.find('=');
    if (eqPos == std::string_view::npos) {
        error = "series-compression must have form SERIES=METHOD";
        return std::nullopt;
    }
    if (eqPos == 0 || eqPos + 1 >= std::size(arg)) {
        error = "series-compression must have non-empty SERIES and METHOD";
        return std::nullopt;
    }

    const auto seriesArg = arg.substr(0, eqPos);
    const auto methodArg = arg.substr(eqPos + 1);
    const auto series = ParseDataSeries(seriesArg);
    if (!series) {
        error = std::format("invalid CRAM data series '{}'", seriesArg);
        return std::nullopt;
    }

    const auto method = ParseCompressionMethod(methodArg);
    if (!method) {
        error = std::format("invalid block compression method '{}'", methodArg);
        return std::nullopt;
    }

    return std::pair{*series, *method};
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

int Runner(int argc, char** argv)
{
    CramWriterConfig config;
    config.BlockCompressionMethod = CramBlockMethod::RANS4X8;
    std::int32_t bgzfThreadsOpt{-1};
    std::int32_t decodeThreadsOpt{-1};
    std::int32_t compressionThreadsOpt{-1};
    bool convertToBamRecord{false};
    const char* inputFile{nullptr};
    const char* outputFile{nullptr};

    for (int i{0}; i < argc; ++i) {
        const std::string_view arg{argv[i]};

        if ((arg == "--help") || (arg == "-h")) {
            PrintUsage();
            return EXIT_SUCCESS;
        }

        if ((arg == "--records-per-slice") && (i + 1 < argc)) {
            const auto parsed{Tools::ParseInteger<std::int32_t>(argv[++i], "records-per-slice")};
            if (!parsed) {
                std::println(stderr, "Error: {}", parsed.error());
                return EXIT_FAILURE;
            }
            if (*parsed <= 0) {
                std::println(stderr, "Error: records-per-slice must be > 0");
                return EXIT_FAILURE;
            }
            config.RecordsPerSlice = *parsed;
            continue;
        }

        if ((arg == "--bgzf-threads") && (i + 1 < argc)) {
            const auto parsed{Tools::ParseInteger<std::int32_t>(argv[++i], "bgzf-threads")};
            if (!parsed) {
                std::println(stderr, "Error: {}", parsed.error());
                return EXIT_FAILURE;
            }
            if (*parsed < 0) {
                std::println(stderr, "Error: bgzf-threads must be >= 0");
                return EXIT_FAILURE;
            }
            bgzfThreadsOpt = *parsed;
            continue;
        }

        if ((arg == "--decode-threads") && (i + 1 < argc)) {
            const auto parsed{Tools::ParseInteger<std::int32_t>(argv[++i], "decode-threads")};
            if (!parsed) {
                std::println(stderr, "Error: {}", parsed.error());
                return EXIT_FAILURE;
            }
            if (*parsed < 0) {
                std::println(stderr, "Error: decode-threads must be >= 0");
                return EXIT_FAILURE;
            }
            decodeThreadsOpt = *parsed;
            continue;
        }

        if (arg == "--convert-to-bam-record") {
            convertToBamRecord = true;
            continue;
        }

        if ((arg == "--compression-threads") && (i + 1 < argc)) {
            const auto parsed{Tools::ParseInteger<std::int32_t>(argv[++i], "compression-threads")};
            if (!parsed) {
                std::println(stderr, "Error: {}", parsed.error());
                return EXIT_FAILURE;
            }
            if (*parsed < 0) {
                std::println(stderr, "Error: compression-threads must be >= 0");
                return EXIT_FAILURE;
            }
            compressionThreadsOpt = *parsed;
            continue;
        }

        if ((arg == "--block-compression") && (i + 1 < argc)) {
            const std::string_view methodArg{argv[++i]};
            const auto method = ParseCompressionMethod(methodArg);
            if (!method) {
                std::println(stderr, "Error: invalid block compression method '{}'", methodArg);
                std::println(stderr,
                             "Valid values: raw,gzip,bzip2,lzma,rans4x8,rans4x16,arith,"
                             "fqzcomp,tok or numeric 0-8");
                return EXIT_FAILURE;
            }
            config.BlockCompressionMethod = *method;
            continue;
        }

        if ((arg == "--series-compression") && (i + 1 < argc)) {
            const std::string_view specArg{argv[++i]};
            std::string error;
            const auto parsed = ParseSeriesCompression(specArg, error);
            if (!parsed) {
                std::println(stderr, "Error: {}", error);
                std::println(stderr,
                             "Valid series include: BF,CF,RI,RL,AP,RG,RN,MF,NS,NP,TS,NF,TL,"
                             "FN,FC,FP,MQ,BA,QS,BS,IN,DL,SC,RS,PD,HC,BB,QQ");
                std::println(stderr,
                             "Valid methods: raw,gzip,bzip2,lzma,rans4x8,rans4x16,arith,"
                             "fqzcomp,tok or numeric 0-8");
                return EXIT_FAILURE;
            }
            config.DataSeriesCompressionMethods[parsed->first] = parsed->second;
            continue;
        }

        if (arg == "--write-crai") {
            config.WriteCrai = true;
            continue;
        }

        if (!std::empty(arg) && arg[0] == '-') {
            std::println(stderr, "Error: unknown option '{}'", arg);
            PrintUsage();
            return EXIT_FAILURE;
        }

        if (inputFile == nullptr) {
            inputFile = argv[i];
        } else if (outputFile == nullptr) {
            outputFile = argv[i];
        } else {
            std::println(stderr, "Error: too many positional arguments");
            PrintUsage();
            return EXIT_FAILURE;
        }
    }

    if ((inputFile == nullptr) || (outputFile == nullptr)) {
        PrintUsage();
        return EXIT_FAILURE;
    }

    const std::filesystem::path inputPath{inputFile};
    const std::filesystem::path outputPath{outputFile};
    constexpr std::int32_t MAX_IO_WORKERS{16};
    constexpr std::int32_t MAX_COMPRESSION_WORKERS{8};
    const std::size_t bgzfWorkers = Tools::ResolveNumWorkers(bgzfThreadsOpt, MAX_IO_WORKERS);
    const std::size_t decodeWorkers = Tools::ResolveNumWorkers(decodeThreadsOpt, MAX_IO_WORKERS);
    config.CompressionWorkers =
        Tools::ResolveNumWorkers(compressionThreadsOpt, MAX_COMPRESSION_WORKERS);

    if (outputPath.extension() != ".cram") {
        std::println(stderr, "Error: output must use .cram extension");
        return EXIT_FAILURE;
    }

    if (inputPath.extension() == ".bam") {
        if (convertToBamRecord || decodeThreadsOpt >= 0) {
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

    std::println(stderr, "Error: input must use .bam or .sam extension");
    return EXIT_FAILURE;
}

}  // namespace Convert
}  // namespace Samoa
}  // namespace PacBio
