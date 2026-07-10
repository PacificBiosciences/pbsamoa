#include "ZmiIndex.hpp"

#include "RecordParser.hpp"

#include "../../PathUtils.hpp"

#include <pbsamoa/PbSamoaLibraryInfo.hpp>
#include <pbsamoa/core/Bgzf.hpp>
#include <pbsamoa/io/ZmiWriter.hpp>

#include <pbcopper/cli2/Interface.h>
#include <pbcopper/cli2/Option.h>
#include <pbcopper/cli2/PositionalArgument.h>
#include <pbcopper/cli2/Results.h>
#include <pbcopper/parallel/ThreadPool.h>

#include <array>
#include <atomic>
#include <exception>
#include <filesystem>
#include <fstream>
#include <ios>
#include <optional>
#include <print>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>

namespace PacBio {
namespace Samoa {
namespace ZmiIndex {
namespace {

constexpr std::size_t BGZF_HEADER_PREFIX{18U};
constexpr std::size_t MAX_DECOMPRESSED_SIZE{65536U};

struct DecompressedBlock
{
    std::uint64_t FileOffset{0};
    std::vector<std::byte> Data;
};

const CLI_v2::Option Threads{
    R"({
    "names" : ["threads"],
    "description" : "Worker threads; 0 = auto (hardware concurrency).",
    "type" : "integer",
    "default" : 0
})"};

const CLI_v2::Option Quiet{
    R"({
    "names" : ["quiet"],
    "description" : "Suppress the summary line written to stderr."
})"};

const CLI_v2::PositionalArgument Input{
    R"({
    "name" : "input",
    "description" : "Input BAM file.",
    "type" : "file"
})"};

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

CLI_v2::Interface CreateInterface()
{
    CLI_v2::Interface interface{"pbsamoa zmi-index",
                                "Build a .zmi sidecar index for IN.bam without rewriting the BAM. "
                                "Output is written to IN.bam.zmi.",
                                LibraryFormattedVersion()};
    // zmi-index keeps its own --threads (0 = auto via hardware_concurrency); the built-in
    // --num-threads resolves 0 to the raw hardware count with different semantics.
    interface.DisableNumThreadsOption();
    interface.AddOptions({Threads, Quiet});
    interface.AddPositionalArguments({Input});
    return interface;
}

int Runner(const CLI_v2::Results& results)
{
    // Exact-width read: never read into std::size_t (ambiguous conversion).
    const std::int32_t threads{results[Threads]};
    if (threads < 0) {
        throw std::runtime_error{"--threads must be >= 0"};
    }
    // Bool flags use copy-init.
    const bool quiet = results[Quiet];

    // CLIv2 does not enforce the required positional-argument count, so guard before indexing.
    const std::vector<std::string>& pos{results.PositionalArguments()};
    if (pos.size() != 1) {
        throw std::runtime_error{"zmi-index requires exactly one argument: <input>"};
    }

    const std::filesystem::path inputPath{pos[0]};

    if (!std::filesystem::exists(inputPath)) {
        throw std::runtime_error{std::string{"zmi-index: input file not found: "} +
                                 inputPath.string()};
    }

    const std::filesystem::path zmiPath{SidecarPath(inputPath, ".zmi")};
    const std::size_t numWorkers{ResolveThreadCount(threads)};

    using Pool = ::PacBio::Parallel::ThreadPool<DecompressedBlock>;
    Pool pool{Pool::Config{.NumThreads = numWorkers, .QueueMultiplier = 4U}};

    std::atomic<bool> readerFailed{false};
    std::exception_ptr readerError;

    const std::jthread reader{[&]() {
        try {
            std::ifstream input{inputPath, std::ios::binary};
            if (!input) {
                throw std::runtime_error{"zmi-index: cannot open input '" + inputPath.string() +
                                         "'"};
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

    // Remove the partial sidecar on any failure so we don't leave a corrupt file behind.
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

        if (!quiet) {
            std::println(stderr, "ZMI index written to {} ({} records)", zmiPath.string(),
                         parser.RecordsEmitted());
        }
    } catch (...) {
        std::filesystem::remove(zmiPath);
        throw;
    }

    return EXIT_SUCCESS;
}

}  // namespace ZmiIndex
}  // namespace Samoa
}  // namespace PacBio
