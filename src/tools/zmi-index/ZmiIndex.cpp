#include "ZmiIndex.hpp"

#include "RecordParser.hpp"

#include "../../PathUtils.hpp"
#include "../CliUtils.hpp"
#include "../ParseUtils.hpp"

#include <pbsamoa/core/Bgzf.hpp>
#include <pbsamoa/io/ZmiWriter.hpp>

#include <pbcopper/parallel/ThreadPool.h>

#include <array>
#include <atomic>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <ios>
#include <print>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

#include <cstddef>
#include <cstdlib>
#include <cstring>

namespace PacBio {
namespace Samoa {
namespace ZmiIndex {
namespace {

constexpr std::size_t BGZF_HEADER_PREFIX{18U};
constexpr std::size_t MAX_DECOMPRESSED_SIZE{65536U};

struct CliArgs
{
    std::filesystem::path InputPath;
    std::int32_t Threads{-1};
    bool Quiet{false};
};

struct DecompressedBlock
{
    std::uint64_t FileOffset{0};
    std::vector<std::byte> Data;
};

void PrintUsage()
{
    std::println(stderr,
                 "Usage: pbsamoa zmi-index [--threads N] [--quiet] IN.bam\n"
                 "\n"
                 "Build a .zmi sidecar index for IN.bam without rewriting the BAM.\n"
                 "Output is written to IN.bam.zmi.");
}

CliArgs ParseArgs(int argc, char** argv)
{
    CliArgs result;
    for (int i{0}; i < argc; ++i) {
        const std::string_view arg{argv[i]};
        if ((arg == "--help") || (arg == "-h")) {
            PrintUsage();
            std::exit(EXIT_SUCCESS);
        } else if (arg == "--quiet") {
            result.Quiet = true;
        } else if (arg == "--threads") {
            result.Threads = ::PacBio::Samoa::Tools::ParseIntegerOrThrow<std::int32_t>(
                ::PacBio::Samoa::Tools::RequireOptionValue(argc, argv, i, "--threads"),
                "--threads");
            if (result.Threads < 1) {
                throw std::runtime_error{"zmi-index: --threads must be >= 1"};
            }
        } else if (arg.starts_with("--")) {
            throw std::runtime_error{"zmi-index: unknown flag '" + std::string{arg} + "'"};
        } else if (std::empty(result.InputPath)) {
            result.InputPath = std::filesystem::path{arg};
        } else {
            throw std::runtime_error{"zmi-index: unexpected positional argument '" +
                                     std::string{arg} + "'"};
        }
    }

    if (std::empty(result.InputPath)) {
        throw std::runtime_error{"zmi-index: input BAM path is required"};
    }
    return result;
}

std::size_t ResolveThreadCount(std::int32_t requested)
{
    if (requested >= 1) {
        return static_cast<std::size_t>(requested);
    }
    const unsigned hw{std::thread::hardware_concurrency()};
    return (hw == 0U) ? std::size_t{1U} : static_cast<std::size_t>(hw);
}

DecompressedBlock DecompressOne(std::vector<std::byte> compressed, BgzfBlockInfo info,
                                std::uint64_t fileOffset)
{
    DecompressedBlock result;
    result.FileOffset = fileOffset;
    result.Data.resize(MAX_DECOMPRESSED_SIZE);
    const std::optional<std::size_t> size{DecompressBgzfBlock(compressed, info, result.Data)};
    if (!size) {
        throw std::runtime_error{"zmi-index: BGZF block decompression failed"};
    }
    result.Data.resize(*size);
    return result;
}

}  // namespace

int Runner(int argc, char** argv)
{
    CliArgs args;
    try {
        args = ParseArgs(argc, argv);
    } catch (const std::exception& e) {
        std::println(stderr, "{}", e.what());
        PrintUsage();
        return EXIT_FAILURE;
    }

    if (!std::filesystem::exists(args.InputPath)) {
        std::println(stderr, "zmi-index: input file not found: {}", args.InputPath.string());
        return EXIT_FAILURE;
    }

    const std::filesystem::path zmiPath{SidecarPath(args.InputPath, ".zmi")};
    const std::size_t numWorkers{ResolveThreadCount(args.Threads)};

    using Pool = ::PacBio::Parallel::ThreadPool<DecompressedBlock>;
    Pool pool{Pool::Config{.NumThreads = numWorkers, .QueueMultiplier = 4U}};

    std::atomic<bool> readerFailed{false};
    std::exception_ptr readerError;

    const std::jthread reader{[&]() {
        try {
            std::ifstream input{args.InputPath, std::ios::binary};
            if (!input) {
                throw std::runtime_error{"zmi-index: cannot open input '" +
                                         args.InputPath.string() + "'"};
            }

            std::array<std::byte, BGZF_HEADER_PREFIX> headerBytes{};
            std::uint64_t fileOffset{0};
            while (true) {
                input.read(reinterpret_cast<char*>(std::data(headerBytes)),
                           static_cast<std::streamsize>(BGZF_HEADER_PREFIX));
                const std::streamsize headerRead{input.gcount()};
                if (headerRead == 0) {
                    break;
                }
                if (headerRead != static_cast<std::streamsize>(BGZF_HEADER_PREFIX)) {
                    throw std::runtime_error{"zmi-index: truncated BGZF header"};
                }

                const std::optional<BgzfBlockInfo> info{ParseBgzfBlockHeader(headerBytes)};
                if (!info) {
                    throw std::runtime_error{"zmi-index: invalid BGZF block header"};
                }

                std::vector<std::byte> blockBytes(info->blockSize);
                std::memcpy(std::data(blockBytes), std::data(headerBytes), BGZF_HEADER_PREFIX);
                const std::size_t remaining{info->blockSize - BGZF_HEADER_PREFIX};
                if (remaining > 0U) {
                    input.read(reinterpret_cast<char*>(std::data(blockBytes) + BGZF_HEADER_PREFIX),
                               static_cast<std::streamsize>(remaining));
                    if (input.gcount() != static_cast<std::streamsize>(remaining)) {
                        throw std::runtime_error{"zmi-index: truncated BGZF block"};
                    }
                }

                if (IsBgzfEofMarker(blockBytes)) {
                    fileOffset += info->blockSize;
                    continue;
                }

                pool.Submit([compressed = std::move(blockBytes), blockInfo = *info,
                             blockFileOffset = fileOffset]() mutable {
                    return DecompressOne(std::move(compressed), blockInfo, blockFileOffset);
                });

                fileOffset += info->blockSize;
            }
        } catch (...) {
            readerError = std::current_exception();
            readerFailed.store(true, std::memory_order_release);
        }

        try {
            pool.Finalize();
        } catch (...) {
            if (!readerFailed.load(std::memory_order_acquire)) {
                readerError = std::current_exception();
                readerFailed.store(true, std::memory_order_release);
            }
        }
    }};

    try {
        ZmiWriter writer{zmiPath};
        RecordParser parser;

        while (pool.ConsumeWith([&](DecompressedBlock&& block) {
            parser.Feed(block.FileOffset, block.Data);
            while (auto entry{parser.NextRecord()}) {
                writer.AddRecord(entry->rgId, entry->zmw,
                                 static_cast<std::int64_t>(entry->virtualOffset.Value()));
            }
        })) {
        }

        if (readerFailed.load(std::memory_order_acquire)) {
            std::rethrow_exception(readerError);
        }

        parser.Finish();
        writer.Close();

        if (!args.Quiet) {
            std::println(stderr, "ZMI index written to {} ({} records)", zmiPath.string(),
                         parser.RecordsEmitted());
        }
    } catch (const std::exception& e) {
        std::println(stderr, "zmi-index: {}", e.what());
        std::filesystem::remove(zmiPath);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

}  // namespace ZmiIndex
}  // namespace Samoa
}  // namespace PacBio
