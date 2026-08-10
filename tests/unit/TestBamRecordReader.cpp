#include "../../src/ParallelUtils.hpp"
#include "TestData.hpp"

#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamRecordReader.hpp>

#include <gtest/gtest.h>

#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

namespace {

template <typename Source>
void ExpectHeaderForwardedFromSingleInput(Source&& source, const BamRecordReaderConfig& config)
{
    const BamRecordReader reader{std::forward<Source>(source), config};
    EXPECT_EQ(reader.Header().Version(), "1.6");
}

}  // namespace

TEST(BamRecordReader, ReadRecordsMatchViewReaderToOwned)
{
    BamRawReader viewReader{tests::DataDir / "spec_example.bam"};
    BamRecordReader recordReader{tests::DataDir / "spec_example.bam",
                                 BamRecordReaderConfig{.DecodeWorkers = 0}};

    std::size_t count{0};
    while (const auto viewRec = viewReader.ReadRecord()) {
        const BamRecord owned{viewRec->ToOwned()};
        const auto decoded{recordReader.ReadRecord()};
        ASSERT_TRUE(decoded) << "at record " << count;
        EXPECT_EQ(decoded->Name(), owned.Name()) << "at record " << count;
        EXPECT_EQ(decoded->Flag(), owned.Flag()) << "at record " << count;
        EXPECT_EQ(decoded->Pos(), owned.Pos()) << "at record " << count;
        EXPECT_EQ(decoded->RefId(), owned.RefId()) << "at record " << count;
        EXPECT_EQ(decoded->MapQ(), owned.MapQ()) << "at record " << count;
        EXPECT_EQ(decoded->Sequence(), owned.Sequence()) << "at record " << count;
        ++count;
    }
    EXPECT_FALSE(recordReader.ReadRecord());
    EXPECT_GT(count, 0U);
}

TEST(BamRecordReader, ParallelChunkSizeClampsAndRoundsUp)
{
    EXPECT_EQ(detail::ParallelChunkSize(0, 0), 1U);
    EXPECT_EQ(detail::ParallelChunkSize(4, 1), 1U);
    EXPECT_EQ(detail::ParallelChunkSize(5, 1), 2U);
    EXPECT_EQ(detail::ParallelChunkSize(1025, 1), 256U);
}

TEST(BamRecordReader, HeaderForwarded)
{
    const BamRecordReader reader{tests::DataDir / "spec_example.bam",
                                 BamRecordReaderConfig{.DecodeWorkers = 0}};
    EXPECT_EQ(reader.Header().Version(), "1.6");
}

TEST(BamRecordReader, SingleInputConstructorsPreserveSingleFileRawReader)
{
    const std::filesystem::path path{tests::DataDir / "spec_example.bam"};
    const BamRecordReaderConfig config{
        .RawReaderConfig =
            {
                .Whitelist = ZmwWhitelist{std::vector<std::int32_t>{42}},
            },
        .DecodeWorkers = 0,
    };

    EXPECT_NO_THROW(ExpectHeaderForwardedFromSingleInput(BamCollection{path}, config));
    EXPECT_NO_THROW(
        ExpectHeaderForwardedFromSingleInput(std::vector<std::filesystem::path>{path}, config));
    EXPECT_NO_THROW(
        ExpectHeaderForwardedFromSingleInput(std::vector<BamFile>{BamFile{path}}, config));
}

TEST(BamRecordReader, HeaderOnlyBamProducesNoRecords)
{
    BamRecordReader reader{tests::DataDir / "header_only.bam",
                           BamRecordReaderConfig{.DecodeWorkers = 0}};
    EXPECT_FALSE(reader.ReadRecord());
}

TEST(BamRecordReader, RangeInterface)
{
    BamRecordReader reader{tests::DataDir / "spec_example.bam",
                           BamRecordReaderConfig{.DecodeWorkers = 0}};
    std::vector<std::string> names;
    for (const auto& rec : reader.Records()) {
        names.emplace_back(rec.Name());
    }
    ASSERT_EQ(std::size(names), 6u);
    EXPECT_EQ(names[0], "r001");
}

TEST(BamRecordReader, ParallelDecodeMatchesSerial)
{
    const auto path{tests::DataDir / "spec_example.bam"};
    BamRecordReader serial{path, BamRecordReaderConfig{.DecodeWorkers = 0}};
    BamRecordReader parallel{path, BamRecordReaderConfig{
                                       .RawReaderConfig = {.BgzfWorkers = 2},
                                       .DecodeWorkers = 2,
                                   }};

    std::size_t count{0};
    while (true) {
        const auto s{serial.ReadRecord()};
        const auto p{parallel.ReadRecord()};
        ASSERT_EQ(s.has_value(), p.has_value()) << "at record " << count;
        if (!s) {
            break;
        }
        EXPECT_EQ(s->Name(), p->Name()) << "at record " << count;
        EXPECT_EQ(s->Flag(), p->Flag()) << "at record " << count;
        EXPECT_EQ(s->Pos(), p->Pos()) << "at record " << count;
        ++count;
    }
    EXPECT_GT(count, 0U);
}

TEST(BamRecordReader, ManyRecordsParallel)
{
    const auto path{tests::DataDir / "many_records.bam"};
    if (!std::filesystem::exists(path)) {
        GTEST_SKIP() << "many_records.bam not generated";
    }

    BamRecordReader reader{path, BamRecordReaderConfig{
                                     .RawReaderConfig = {.BgzfWorkers = 2},
                                     .DecodeWorkers = 2,
                                 }};

    std::size_t count{0};
    for (const auto& rec : reader.Records()) {
        EXPECT_FALSE(std::empty(rec.Name()));
        ++count;
    }
    EXPECT_EQ(count, 500u);
}

TEST(BamRecordReader, ChunkedDecodeMatchesSerialAcrossWorkerCounts)
{
    // The fixture and small batch budget exercise non-divisible tail chunks across worker counts.
    const auto path{tests::DataDir / "benchmark.bam"};

    std::vector<BamRecord> expected;
    BamRecordReader serial{path, BamRecordReaderConfig{.DecodeWorkers = 0}};
    for (const auto& rec : serial.Records()) {
        expected.push_back(rec);
    }
    ASSERT_EQ(std::size(expected), 10000U);

    for (const std::size_t workers : {1U, 2U, 3U, 5U, 8U}) {
        for (const ByteLimit budget : {ByteLimit{4 * 1024 * 1024}, ByteLimit{64 * 1024}}) {
            BamRecordReader parallel{path, BamRecordReaderConfig{
                                               .DecodeWorkers = workers,
                                               .BatchBudget = budget,
                                           }};

            std::size_t i{0};
            for (const auto& rec : parallel.Records()) {
                ASSERT_LT(i, std::size(expected))
                    << "extra record at " << i << ", workers=" << workers;
                EXPECT_EQ(rec.Name(), expected[i].Name()) << "at " << i << ", workers=" << workers;
                EXPECT_EQ(rec.Flag(), expected[i].Flag()) << "at " << i << ", workers=" << workers;
                EXPECT_EQ(rec.Pos(), expected[i].Pos()) << "at " << i << ", workers=" << workers;
                EXPECT_EQ(rec.RefId(), expected[i].RefId())
                    << "at " << i << ", workers=" << workers;
                EXPECT_EQ(rec.Sequence(), expected[i].Sequence())
                    << "at " << i << ", workers=" << workers;
                ++i;
            }
            EXPECT_EQ(i, std::size(expected)) << "missing records, workers=" << workers;
        }
    }
}

TEST(BamRecordReader, TagFilterDropTags)
{
    BamRecordReader noFilter{tests::DataDir / "spec_example.bam",
                             BamRecordReaderConfig{.DecodeWorkers = 0}};
    BamRecordReader withFilter{tests::DataDir / "spec_example.bam",
                               BamRecordReaderConfig{
                                   .DecodeWorkers = 0,
                                   .TagFilter = DropTags{TagKey{'N', 'M'}},
                               }};

    while (true) {
        const auto full{noFilter.ReadRecord()};
        const auto filtered{withFilter.ReadRecord()};
        ASSERT_EQ(full.has_value(), filtered.has_value());
        if (!full) {
            break;
        }
        // Core fields match
        EXPECT_EQ(full->Name(), filtered->Name());
        EXPECT_EQ(full->Pos(), filtered->Pos());
        // Filtered record should not have the NM tag
        EXPECT_FALSE(filtered->Tags().Contains(TagKey{'N', 'M'}));
    }
}

TEST(BamRecordReader, TagFilterKeepTags)
{
    BamRecordReader noFilter{tests::DataDir / "spec_example.bam",
                             BamRecordReaderConfig{.DecodeWorkers = 0}};
    BamRecordReader withFilter{tests::DataDir / "spec_example.bam",
                               BamRecordReaderConfig{
                                   .DecodeWorkers = 0,
                                   .TagFilter = KeepTags{TagKey{'N', 'M'}},
                               }};

    while (true) {
        const auto full{noFilter.ReadRecord()};
        const auto filtered{withFilter.ReadRecord()};
        ASSERT_EQ(full.has_value(), filtered.has_value());
        if (!full) {
            break;
        }
        // Core fields match
        EXPECT_EQ(full->Name(), filtered->Name());
        EXPECT_EQ(full->Pos(), filtered->Pos());
        // Filtered record should only have the NM tag (if it was present originally)
        if (full->Tags().Contains(TagKey{'N', 'M'})) {
            EXPECT_TRUE(filtered->Tags().Contains(TagKey{'N', 'M'}));
        }
        // Filtered record should NOT have other tags like RG
        if (full->Tags().Contains(TagKey{'R', 'G'})) {
            EXPECT_FALSE(filtered->Tags().Contains(TagKey{'R', 'G'}));
        }
    }
}

TEST(BamRecordReader, EarlyDestruction)
{
    const auto path{tests::DataDir / "many_records.bam"};
    if (!std::filesystem::exists(path)) {
        GTEST_SKIP() << "many_records.bam not generated";
    }

    {
        BamRecordReader reader{path, BamRecordReaderConfig{
                                         .RawReaderConfig = {.BgzfWorkers = 2},
                                         .DecodeWorkers = 2,
                                     }};
        // Read just one record, then destroy
        const auto rec{reader.ReadRecord()};
        EXPECT_TRUE(rec);
    }
    SUCCEED();
}

// Regression: ProducerLoop's try_push retry moved from the record twice,
// producing phantom default-constructed BamRecords when the SPSC queue was full.
TEST(BamRecordReader, NoPhantomRecordsOnQueueBackpressure)
{
    const auto path{tests::DataDir / "benchmark.bam"};
    if (!std::filesystem::exists(path)) {
        GTEST_SKIP() << "benchmark.bam not generated";
    }

    // Ground truth via raw reader
    std::size_t rawCount{0};
    {
        BamRawReader raw{path};
        while (raw.ReadRecord()) {
            ++rawCount;
        }
    }
    ASSERT_GT(rawCount, 0U);

    // Tiny queue (capacity=2) forces heavy try_push backpressure.
    // Before the fix, every failed try_push would move-from the record,
    // and the retry would push an empty (phantom) record.
    BamRecordReader reader{path, BamRecordReaderConfig{
                                     .DecodeWorkers = 2,
                                     .OutputCapacity = 2,
                                 }};

    std::size_t count{0};
    while (const auto rec = reader.ReadRecord()) {
        ASSERT_FALSE(std::empty(rec->Name())) << "phantom record at index " << count
                                              << " (empty name, tags=" << rec->Tags().Size() << ")";
        ++count;
    }
    EXPECT_EQ(count, rawCount);
}

TEST(BamRecordReader, NumZmwsWithoutZmiReturnsNegativeOne)
{
    // spec_example.bam has no ZMI, so NumZmws() should return -1.
    const BamRecordReader reader{tests::DataDir / "spec_example.bam"};
    EXPECT_EQ(reader.NumZmws(), -1);
}

}  // namespace Samoa
}  // namespace PacBio
