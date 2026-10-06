#include "TestTempDir.hpp"

#include <pbsamoa/core/Bgzf.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <format>
#include <mutex>
#include <optional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <thread>
#include <vector>

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

namespace PacBio {
namespace Samoa {
namespace {

struct CallbackCapture
{
    std::mutex* mutex;
    std::vector<std::int64_t>* offsets;
    std::vector<std::vector<std::byte>>* payloads;

    void operator()(std::int64_t offset, std::span<const std::byte> rawData) const
    {
        const std::lock_guard lock{*mutex};
        offsets->push_back(offset);
        payloads->emplace_back(std::ranges::begin(rawData), std::ranges::end(rawData));
    }
};

void CheckSharedOutput(const std::filesystem::path& path, const std::size_t i,
                       const std::size_t records, const std::span<const std::int64_t> offsets,
                       const std::span<const std::byte> markers)
{
    BgzfReader reader{path};
    ASSERT_TRUE(reader.HasEofMarker());
    ASSERT_EQ(std::size(offsets), records);
    std::vector<std::byte> output;
    std::vector<std::byte> block(65536);
    std::vector<std::size_t> recordPositions;
    while (true) {
        const auto size{reader.ReadBlock(block)};
        ASSERT_TRUE(size);
        if (*size == 0) {
            break;
        }
        const VirtualOffset blockOffset{reader.Tell()};
        for (const std::int64_t offset : offsets) {
            const VirtualOffset recordOffset{static_cast<std::uint64_t>(offset)};
            if (recordOffset.BlockOffset() == blockOffset.BlockOffset()) {
                recordPositions.push_back(std::size(output) + recordOffset.WithinBlockOffset());
            }
        }
        output.insert(std::end(output), std::begin(block), std::begin(block) + *size);
    }
    ASSERT_EQ(std::size(output), records * 70000);
    ASSERT_EQ(std::size(recordPositions), records);
    for (std::size_t round{0}; round < records; ++round) {
        EXPECT_EQ(recordPositions[round], round * 70000);
        const std::byte value{static_cast<std::byte>(i + round)};
        EXPECT_EQ(markers[round], value);
        EXPECT_TRUE(std::ranges::all_of(std::span{output}.subspan(round * 70000, 70000),
                                        [value](const std::byte byte) { return byte == value; }));
    }
}

}  // namespace

class BgzfWriterTest : public ::testing::Test
{
protected:
    tests::TempDirGuard tempDir_;
    std::filesystem::path tmpPath_{};

    void SetUp() override
    {
        tempDir_.Reset("bgzf_writer");
        tmpPath_ = tempDir_.File("output.bgzf");
    }
};

TEST(BgzfWriterConfigTest, DefaultValues)
{
    const BgzfWriterConfig config{};
    EXPECT_EQ(config.CompressionLevel, 6);
    EXPECT_EQ(config.BgzfWorkers, 4U);
    EXPECT_EQ(config.InputQueueCapacity, 256U);
    EXPECT_EQ(config.BlocksPerBatch, 32);
    EXPECT_FALSE(config.UseTempFile);
}

TEST(BgzfWriteMetricsTest, DefaultValues)
{
    const BgzfWriteMetrics metrics{};
    EXPECT_EQ(metrics.CallerStalls, 0U);
    EXPECT_EQ(metrics.PackerStalls, 0U);
    EXPECT_EQ(metrics.CompressNs, 0U);
    EXPECT_EQ(metrics.BytesCompressed, 0U);
    EXPECT_EQ(metrics.BlocksWritten, 0U);
    EXPECT_EQ(metrics.IoWriteNs, 0U);
    EXPECT_EQ(metrics.WriterStalls, 0U);
    EXPECT_EQ(metrics.CallbackNs, 0U);
}

TEST_F(BgzfWriterTest, WriteAndReadBack)
{
    const std::string original{"Hello, BGZF round-trip!"};

    {
        BgzfWriter writer{tmpPath_};
        writer.Write(std::as_bytes(std::span{std::data(original), std::size(original)}));
    }

    BgzfReader reader{tmpPath_};
    EXPECT_TRUE(reader.HasEofMarker());

    std::vector<std::byte> buffer{};
    buffer.resize(65536U);
    const std::optional<std::size_t> bytesRead{reader.ReadBlock(buffer)};
    ASSERT_TRUE(bytesRead);
    ASSERT_EQ(*bytesRead, std::size(original));

    const std::string result{reinterpret_cast<const char*>(std::data(buffer)), *bytesRead};
    EXPECT_EQ(result, original);
}

TEST_F(BgzfWriterTest, LargeDataSpansMultipleBlocks)
{
    std::vector<std::byte> input{};
    input.resize(200000U);
    for (std::size_t i{0}; i < std::size(input); ++i) {
        input[i] = static_cast<std::byte>(i & 0xFFU);
    }

    {
        BgzfWriter writer{tmpPath_};
        writer.Write(input);
    }

    BgzfReader reader{tmpPath_};
    std::vector<std::byte> output{};
    std::vector<std::byte> buffer{};
    buffer.resize(65536U);
    while (true) {
        const std::optional<std::size_t> bytesRead{reader.ReadBlock(buffer)};
        if (!bytesRead || (*bytesRead == 0U)) {
            break;
        }
        output.insert(std::ranges::end(output), std::ranges::begin(buffer),
                      std::ranges::begin(buffer) + *bytesRead);
    }

    ASSERT_EQ(std::size(output), std::size(input));
    EXPECT_TRUE(std::ranges::equal(output, input));
}

TEST_F(BgzfWriterTest, EmptyFileHasEofMarker)
{
    {
        const BgzfWriter writer{tmpPath_};
    }

    BgzfReader reader{tmpPath_};
    EXPECT_TRUE(reader.HasEofMarker());

    std::vector<std::byte> buffer{};
    buffer.resize(65536U);
    const std::optional<std::size_t> bytesRead{reader.ReadBlock(buffer)};
    EXPECT_TRUE(!bytesRead || (*bytesRead == 0U));
}

TEST_F(BgzfWriterTest, MultipleSmallWrites)
{
    const std::string part1{"Hello, "};
    const std::string part2{"world!"};

    {
        BgzfWriter writer{tmpPath_};
        writer.Write(std::as_bytes(std::span{std::data(part1), std::size(part1)}));
        writer.Write(std::as_bytes(std::span{std::data(part2), std::size(part2)}));
    }

    BgzfReader reader{tmpPath_};
    std::vector<std::byte> buffer{};
    buffer.resize(65536U);
    const std::optional<std::size_t> bytesRead{reader.ReadBlock(buffer)};
    ASSERT_TRUE(bytesRead);

    const std::string result{reinterpret_cast<const char*>(std::data(buffer)), *bytesRead};
    EXPECT_EQ(result, "Hello, world!");
}

TEST_F(BgzfWriterTest, ConstructWithConfig)
{
    const BgzfWriterConfig config{
        .CompressionLevel = 4,
        .BgzfWorkers = 2,
    };
    BgzfWriter writer{tmpPath_, config};
    writer.Close();
    EXPECT_TRUE(std::filesystem::exists(tmpPath_));
}

TEST_F(BgzfWriterTest, UseTempFileAtomicWrite)
{
    const std::string original{"Hello, temp-file BGZF!"};
    const BgzfWriterConfig config{
        .UseTempFile = true,
    };

    {
        BgzfWriter writer{tmpPath_, config};
        writer.Write(std::as_bytes(std::span{std::data(original), std::size(original)}));
        EXPECT_FALSE(std::filesystem::exists(tmpPath_));
        writer.Close();
    }

    EXPECT_TRUE(std::filesystem::exists(tmpPath_));
}

TEST_F(BgzfWriterTest, WriteWithCallback)
{
    std::vector<std::int64_t> callbackOffsets;
    std::vector<std::vector<std::byte>> callbackData;
    std::mutex mu;

    BgzfWriter writer{tmpPath_, BgzfWriterConfig{.BgzfWorkers = 2}};
    writer.SetCallback(CallbackCapture{&mu, &callbackOffsets, &callbackData});

    for (int i{0}; i < 3; ++i) {
        std::vector<std::byte> record(100);
        std::ranges::fill(record, static_cast<std::byte>(i));
        const PendingCallback cb{
            .rawData = record,
            .active = true,
        };
        writer.Write(record, cb);
    }
    writer.Close();

    const std::lock_guard lock{mu};
    ASSERT_EQ(std::size(callbackOffsets), 3U);
    ASSERT_EQ(std::size(callbackData), 3U);
    for (std::size_t i{1}; i < std::size(callbackOffsets); ++i) {
        EXPECT_GE(callbackOffsets[i], callbackOffsets[i - 1]);
    }
}

TEST_F(BgzfWriterTest, MetricsPopulatedAfterWrite)
{
    BgzfWriter writer{tmpPath_, BgzfWriterConfig{.BgzfWorkers = 2}};
    std::vector<std::byte> chunk(32000U);
    for (int i{0}; i < 20; ++i) {
        writer.Write(chunk);
    }
    writer.Close();

    const BgzfWriteMetrics metrics{writer.GetMetrics()};
    EXPECT_GT(metrics.BlocksWritten, 0U);
    EXPECT_GT(metrics.BytesCompressed, 0U);
    EXPECT_GT(metrics.CompressNs, 0U);
    EXPECT_GT(metrics.IoWriteNs, 0U);
}

TEST_F(BgzfWriterTest, OversizedRecordIsSplitAcrossBlocksAndCallbackFiresOnce)
{
    std::vector<std::byte> input(200000U);
    for (std::size_t i{0}; i < std::size(input); ++i) {
        input[i] = static_cast<std::byte>(i & 0xFFU);
    }

    std::vector<std::vector<std::byte>> callbackData;
    std::vector<std::int64_t> callbackOffsets;
    std::mutex mu;

    BgzfWriter writer{tmpPath_, BgzfWriterConfig{
                                    .BgzfWorkers = 2,
                                    .BlocksPerBatch = 1,
                                }};
    writer.SetCallback(CallbackCapture{&mu, &callbackOffsets, &callbackData});

    const PendingCallback cb{
        .rawData = input,
        .active = true,
    };
    writer.Write(input, cb);
    writer.Close();

    {
        const std::lock_guard lock{mu};
        ASSERT_EQ(std::size(callbackOffsets), 1U);
        ASSERT_EQ(std::size(callbackData), 1U);
        EXPECT_TRUE(std::ranges::equal(callbackData.front(), input));
    }

    BgzfReader reader{tmpPath_};
    std::vector<std::byte> output{};
    std::vector<std::byte> buffer(65536U);
    while (true) {
        const std::optional<std::size_t> bytesRead{reader.ReadBlock(buffer)};
        ASSERT_TRUE(bytesRead);
        if (*bytesRead == 0U) {
            break;
        }
        output.insert(std::ranges::end(output), std::ranges::begin(buffer),
                      std::ranges::begin(buffer) + static_cast<std::ptrdiff_t>(*bytesRead));
    }
    EXPECT_TRUE(std::ranges::equal(output, input));

    const BgzfWriteMetrics metrics{writer.GetMetrics()};
    EXPECT_GE(metrics.BlocksWritten, 2U);
}

TEST_F(BgzfWriterTest, SharedPoolPreservesInterleavedDataAndOffsetsWithBoundedQueues)
{
    for (const std::size_t workers : {1UZ, 3UZ}) {
        const auto pool{std::make_shared<BgzfWriterPool>(workers, 1)};
        std::vector<std::vector<std::int64_t>> offsets(64);
        std::vector<std::vector<std::byte>> markers(64);
        std::vector<std::unique_ptr<BgzfWriter>> writers;
        for (std::size_t i{0}; i < std::size(offsets); ++i) {
            writers.push_back(std::make_unique<BgzfWriter>(
                tempDir_.File(std::format("shared-{}-{}.bgzf", workers, i)),
                BgzfWriterConfig{.CompressionLevel = (i % 2 == 0) ? 1 : 9,
                                 .BgzfWorkers = 0,
                                 .InputQueueCapacity = 0,
                                 .BlocksPerBatch = 1,
                                 .Pool = pool}));
            writers.back()->SetCallback(
                [&, i](const std::int64_t offset, const std::span<const std::byte> marker) {
                    offsets[i].push_back(offset);
                    markers[i].push_back(marker.front());
                });
        }
        // One destination closes early, one stays empty, and records cross block boundaries.
        for (std::size_t round{0}; round < 3; ++round) {
            for (std::size_t i{round == 0 ? 0UZ : 1UZ}; i + 1 < std::size(writers); ++i) {
                const std::byte value{static_cast<std::byte>(i + round)};
                writers[i]->Write(std::vector<std::byte>(70000, value),
                                  PendingCallback{.rawData = {value}, .active = true});
            }
            writers.front()->Close();
        }
        for (std::size_t i{0}; i < std::size(writers); ++i) {
            writers[i]->Close();
            const std::size_t records{(i + 1 == std::size(writers)) ? 0UZ : (i == 0 ? 1UZ : 3UZ)};
            CheckSharedOutput(tempDir_.File(std::format("shared-{}-{}.bgzf", workers, i)), i,
                              records, offsets[i], markers[i]);
            EXPECT_LE(writers[i]->GetMetrics().Pool.PeakActiveWorkers, workers);
        }
    }
}

TEST_F(BgzfWriterTest, SharedPoolSurvivesWriterFailureAndOwnerRelease)
{
    auto pool{std::make_shared<BgzfWriterPool>(1, 1)};
    const auto goodPath{tempDir_.File("good.bgzf")};
    BgzfWriter good{goodPath, {.BlocksPerBatch = 1, .Pool = pool}};
    BgzfWriter bad{tmpPath_, {.BlocksPerBatch = 1, .Pool = pool}};
    pool.reset();
    bad.SetCallback([](std::int64_t, std::span<const std::byte>) {
        throw std::runtime_error{"callback failure"};
    });
    bad.Write(std::vector<std::byte>(70000), {.active = true});
    EXPECT_THROW(bad.Close(), std::runtime_error);
    EXPECT_NO_THROW(bad.Close());
    good.Write(std::vector<std::byte>(70000, std::byte{42}));
    EXPECT_NO_THROW(good.Close());
    EXPECT_TRUE(BgzfReader{goodPath}.HasEofMarker());
    EXPECT_THROW(good.Write(std::vector<std::byte>(1)), std::runtime_error);
}

TEST_F(BgzfWriterTest, SharedPoolMoveAssignmentDrainsPreviousWriter)
{
    const auto pool{std::make_shared<BgzfWriterPool>(1, 1)};
    const auto secondPath{tempDir_.File("second.bgzf")};
    BgzfWriter first{tmpPath_, {.BlocksPerBatch = 1, .Pool = pool}};
    BgzfWriter second{secondPath, {.BlocksPerBatch = 1, .Pool = pool}};
    first.Write(std::vector<std::byte>(70000, std::byte{1}));
    second.Write(std::vector<std::byte>(70000, std::byte{2}));
    first = std::move(second);
    EXPECT_TRUE(BgzfReader{tmpPath_}.HasEofMarker());
    first.Close();
    EXPECT_TRUE(BgzfReader{secondPath}.HasEofMarker());
}

class BgzfBlockWriterTest : public ::testing::Test
{
protected:
    tests::TempDirGuard tempDir_;
    std::filesystem::path tmpPath_{};

    void SetUp() override
    {
        tempDir_.Reset("bgzf_block_writer");
        tmpPath_ = tempDir_.File("output.bgzf");
    }

    std::vector<std::byte> ReadAll() const
    {
        BgzfReader reader{tmpPath_};
        std::vector<std::byte> output{};
        std::vector<std::byte> buffer(65536U);
        while (true) {
            const std::optional<std::size_t> bytesRead{reader.ReadBlock(buffer)};
            if (!bytesRead || (*bytesRead == 0U)) {
                break;
            }
            output.insert(std::ranges::end(output), std::ranges::begin(buffer),
                          std::ranges::begin(buffer) + *bytesRead);
        }
        return output;
    }
};

TEST(BgzfBlockEncoderTest, RejectsInvalidCompressionLevel)
{
    EXPECT_THROW(BgzfBlockEncoder{0}, std::invalid_argument);
    EXPECT_THROW(BgzfBlockEncoder{13}, std::invalid_argument);
}

TEST(BgzfBlockEncoderTest, FinishWithoutDataEmitsNothing)
{
    BgzfBlockEncoder encoder{};
    encoder.Finish();
    EXPECT_TRUE(std::empty(encoder.TakeBlocks()));
}

TEST_F(BgzfBlockWriterTest, BlocksFromSeveralEncodersReadBackInWriteOrder)
{
    std::vector<std::byte> expected{};
    {
        BgzfBlockWriter writer{tmpPath_};
        for (std::size_t part{0}; part < 3; ++part) {
            // 150000 bytes fill two full blocks and a part of a third one.
            std::vector<std::byte> data(150000U);
            for (std::size_t i{0}; i < std::size(data); ++i) {
                data[i] = static_cast<std::byte>((i * 7U + part) & 0xFFU);
            }
            BgzfBlockEncoder encoder{2};
            encoder.Write(data);
            encoder.Finish();
            writer.Write(encoder.TakeBlocks());
            expected.insert(std::ranges::end(expected), std::ranges::begin(data),
                            std::ranges::end(data));
        }
        writer.Close();
    }

    EXPECT_TRUE(BgzfReader{tmpPath_}.HasEofMarker());
    EXPECT_TRUE(std::ranges::equal(ReadAll(), expected));
}

TEST_F(BgzfBlockWriterTest, IncompressibleDataFitsInBlocks)
{
    // A xorshift sequence does not compress; each full block must still fit.
    std::vector<std::byte> data(200000U);
    std::uint64_t state{0x9E3779B97F4A7C15ULL};
    for (std::byte& b : data) {
        state ^= state << 13U;
        state ^= state >> 7U;
        state ^= state << 17U;
        b = static_cast<std::byte>(state & 0xFFU);
    }
    {
        BgzfBlockWriter writer{tmpPath_};
        BgzfBlockEncoder encoder{12};
        encoder.Write(data);
        encoder.Finish();
        writer.Write(encoder.TakeBlocks());
    }

    EXPECT_TRUE(std::ranges::equal(ReadAll(), data));
}

TEST_F(BgzfBlockWriterTest, EmptyFileHasEofMarker)
{
    {
        const BgzfBlockWriter writer{tmpPath_};
    }

    EXPECT_TRUE(BgzfReader{tmpPath_}.HasEofMarker());
    EXPECT_TRUE(std::empty(ReadAll()));
}

TEST_F(BgzfBlockWriterTest, UseTempFileAtomicWrite)
{
    {
        BgzfBlockWriter writer{tmpPath_, BgzfBlockWriterConfig{.UseTempFile = true}};
        EXPECT_FALSE(std::filesystem::exists(tmpPath_));
        writer.Close();
    }

    EXPECT_TRUE(std::filesystem::exists(tmpPath_));
}

TEST_F(BgzfBlockWriterTest, WriteAfterCloseThrows)
{
    BgzfBlockWriter writer{tmpPath_};
    writer.Close();
    EXPECT_THROW(writer.Write(std::vector<std::byte>(4U)), std::runtime_error);
}

TEST_F(BgzfBlockWriterTest, WriteBlocksWhenOutputStalls)
{
    // A FIFO with no active reader stalls the IO thread. The FIFO buffer (64 KiB at most on
    // Linux and macOS) is smaller than one chunk, so at most CAP / CHUNK writes can return.
    constexpr std::size_t CHUNK{std::size_t{256} * 1024};
    constexpr std::size_t CAP{4 * CHUNK};
    constexpr std::int32_t CHUNKS{16};
    const std::filesystem::path fifoPath{tempDir_.File("stall.fifo")};
    ASSERT_EQ(::mkfifo(fifoPath.c_str(), 0600), 0);
    // Open the read end first without blocking, so that the writer can open the FIFO.
    const int readFd{::open(fifoPath.c_str(), O_RDONLY | O_NONBLOCK)};
    ASSERT_GE(readFd, 0);

    std::vector<std::byte> expected{};
    {
        BgzfBlockWriter writer{fifoPath, BgzfBlockWriterConfig{.MaxPendingBytes = CAP}};
        ASSERT_EQ(::fcntl(readFd, F_SETFL, 0), 0);

        std::atomic<std::int32_t> written{0};
        std::thread producer{[&writer, &written] {
            for (std::int32_t i{0}; i < CHUNKS; ++i) {
                writer.Write(std::vector<std::byte>(CHUNK, static_cast<std::byte>(i)));
                written.fetch_add(1, std::memory_order_relaxed);
            }
        }};
        // A longer wait only lets more writes return, so the bound below cannot fail by timing.
        std::this_thread::sleep_for(std::chrono::milliseconds{200});
        EXPECT_LE(written.load(std::memory_order_relaxed), CAP / CHUNK);

        std::vector<std::byte> drained{};
        std::thread reader{[readFd, &drained] {
            std::vector<std::byte> buffer(CHUNK);
            while (true) {
                const ssize_t n{::read(readFd, std::data(buffer), std::size(buffer))};
                if (n <= 0) {
                    break;
                }
                drained.insert(std::ranges::end(drained), std::ranges::begin(buffer),
                               std::ranges::begin(buffer) + n);
            }
        }};
        producer.join();
        EXPECT_EQ(written.load(std::memory_order_relaxed), CHUNKS);
        writer.Close();
        reader.join();

        for (std::int32_t i{0}; i < CHUNKS; ++i) {
            expected.insert(std::ranges::end(expected), CHUNK, static_cast<std::byte>(i));
        }
        expected.insert(std::ranges::end(expected), std::ranges::begin(BGZF_EOF_MARKER),
                        std::ranges::end(BGZF_EOF_MARKER));
        EXPECT_TRUE(std::ranges::equal(drained, expected));
    }
    ::close(readFd);
}

TEST_F(BgzfBlockWriterTest, WriteAboveCapPassesWhenQueueIsEmpty)
{
    {
        BgzfBlockWriter writer{tmpPath_, BgzfBlockWriterConfig{.MaxPendingBytes = 0}};
        for (std::size_t part{0}; part < 3; ++part) {
            BgzfBlockEncoder encoder{};
            encoder.Write(std::vector<std::byte>(100000U, static_cast<std::byte>(part)));
            encoder.Finish();
            writer.Write(encoder.TakeBlocks());
        }
    }

    EXPECT_EQ(std::size(ReadAll()), 300000U);
}

}  // namespace Samoa
}  // namespace PacBio
