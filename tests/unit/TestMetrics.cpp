#include "TestData.hpp"

#include <pbsamoa/core/Metrics.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamRecordReader.hpp>

#include <gtest/gtest.h>

#include <filesystem>

#include <cstddef>

namespace PacBio {
namespace Samoa {
namespace {

template <typename Reader>
std::size_t ReadAllRecords(Reader& reader)
{
    std::size_t count{0};
    while (reader.ReadRecord().has_value()) {
        ++count;
    }
    return count;
}

}  // namespace

TEST(Metrics, BgzfMetricsDefaultZero)
{
    const BgzfMetrics m{};
    EXPECT_EQ(m.BytesRead, 0U);
    EXPECT_EQ(m.BlocksRead, 0U);
    EXPECT_EQ(m.RecordsProduced, 0U);
    EXPECT_EQ(m.RecordsConsumed, 0U);
    EXPECT_EQ(m.IoStalls, 0U);
    EXPECT_EQ(m.ConsumerStalls, 0U);
    EXPECT_EQ(m.ReaderStalls, 0U);
    EXPECT_EQ(m.IoReadNs, 0U);
    EXPECT_EQ(m.DecompressNs, 0U);
    EXPECT_EQ(m.RecordParseNs, 0U);
}

TEST(Metrics, SyncModeGetMetrics)
{
    BamRawReader reader{tests::DataDir / "spec_example.bam"};

    // Before reading, counters should be zero
    const BgzfMetrics before{reader.GetMetrics()};
    EXPECT_EQ(before.RecordsConsumed, 0U);

    // Read all records
    const std::size_t count{ReadAllRecords(reader)};
    EXPECT_GT(count, 0U);

    // After reading, RecordsConsumed should match
    const BgzfMetrics after{reader.GetMetrics()};
    EXPECT_EQ(after.RecordsConsumed, count);
}

TEST(Metrics, PipelineModeGetMetrics)
{
    BamRawReader reader{tests::DataDir / "spec_example.bam", BamRawReaderConfig{.BgzfWorkers = 2}};

    // Read all records
    const std::size_t count{ReadAllRecords(reader)};
    EXPECT_GT(count, 0U);

    const BgzfMetrics m{reader.GetMetrics()};

    // Throughput counters must be positive
    EXPECT_GT(m.BytesRead, 0U);
    EXPECT_GT(m.BlocksRead, 0U);
    EXPECT_GT(m.BytesDecompressed, 0U);

    // All records should have been produced and consumed
    EXPECT_EQ(m.RecordsProduced, count);
    EXPECT_EQ(m.RecordsConsumed, count);

    // Timing must be positive (at least some IO and decompression happened)
    EXPECT_GT(m.IoReadNs, 0U);
    EXPECT_GT(m.DecompressNs, 0U);
    EXPECT_GT(m.RecordParseNs, 0U);
}

TEST(Metrics, ReaderMetricsComposite)
{
    BamRecordReader reader{tests::DataDir / "spec_example.bam",
                           BamRecordReaderConfig{
                               .RawReaderConfig = {.BgzfWorkers = 2},
                               .DecodeWorkers = 2,
                           }};

    // Read all records
    const std::size_t count{ReadAllRecords(reader)};
    EXPECT_GT(count, 0U);

    const ReaderMetrics m{reader.GetMetrics()};

    // Flags should reflect configuration
    EXPECT_TRUE(m.ParallelBgzf);
    EXPECT_TRUE(m.ParallelDecode);
    EXPECT_EQ(m.TotalRecordsRead, count);

    // BGZF layer should have positive counters
    EXPECT_GT(m.Bgzf.BytesRead, 0U);
    EXPECT_GT(m.Bgzf.RecordsProduced, 0U);

    // Decode layer should have positive counters
    EXPECT_GT(m.Decode.BatchesDecoded, 0U);
    EXPECT_EQ(m.Decode.RecordsDecoded, count);
    EXPECT_EQ(m.Decode.RecordsConsumed, count);
}

TEST(Metrics, ReaderMetricsSerialMode)
{
    BamRecordReader reader{tests::DataDir / "spec_example.bam", BamRecordReaderConfig{
                                                                    .DecodeWorkers = 0,
                                                                }};

    const std::size_t count{ReadAllRecords(reader)};
    EXPECT_GT(count, 0U);

    const ReaderMetrics m{reader.GetMetrics()};

    // Flags should reflect serial configuration
    EXPECT_FALSE(m.ParallelBgzf);
    EXPECT_FALSE(m.ParallelDecode);
    EXPECT_EQ(m.TotalRecordsRead, count);
}

TEST(Metrics, LivePollingDuringRead)
{
    const auto path{tests::DataDir / "many_records.bam"};
    if (!std::filesystem::exists(path)) {
        GTEST_SKIP() << "many_records.bam not generated";
    }

    BamRecordReader reader{path, BamRecordReaderConfig{
                                     .RawReaderConfig = {.BgzfWorkers = 2},
                                     .DecodeWorkers = 2,
                                 }};

    // Read a few records and check metrics mid-stream
    for (int i{0}; i < 10; ++i) {
        const auto rec{reader.ReadRecord()};
        ASSERT_TRUE(rec.has_value());
    }

    const ReaderMetrics mid{reader.GetMetrics()};
    EXPECT_GE(mid.TotalRecordsRead, 10U);
    EXPECT_GE(mid.Decode.RecordsDecoded, 10U);

    // Finish reading
    while (reader.ReadRecord().has_value()) {
    }

    const ReaderMetrics finalMetrics{reader.GetMetrics()};
    EXPECT_GT(finalMetrics.TotalRecordsRead, mid.TotalRecordsRead);
}

}  // namespace Samoa
}  // namespace PacBio
