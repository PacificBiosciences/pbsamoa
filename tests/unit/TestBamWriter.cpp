#include "TestData.hpp"

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/Bgzf.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamWriter.hpp>

#include <gtest/gtest.h>

#include <filesystem>
#include <format>
#include <string>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {

class BamWriterTest : public ::testing::Test
{
protected:
    std::filesystem::path tmpPath;

    void SetUp() override
    {
        tmpPath = std::filesystem::temp_directory_path() / "pbsamoa_writer_test.bam";
    }

    void TearDown() override { std::filesystem::remove(tmpPath); }

    static SamHeader MakeMinimalHeader()
    {
        SamHeader header;
        header.SetVersion("1.6");
        header.SetSortOrder("unknown");
        header.AddReferenceSequence(ReferenceSequence{"ref", 1000});
        return header;
    }

    static BamRecord MakeTestRecord()
    {
        BamRecord record;
        record.Name("read1")
            .Flag(0)
            .RefId(0)
            .Pos(100)
            .MapQ(30)
            .Cigar({CigarOp{CigarOpType::M, 10}})
            .NextRefId(-1)
            .NextPos(-1)
            .Tlen(0)
            .Sequence("ACGTACGTAC")
            .Qualities({30, 30, 30, 30, 30, 30, 30, 30, 30, 30});
        return record;
    }
};

TEST(BamWriterConfigTest, CompositeDefaults)
{
    const BamWriterConfig config{};
    EXPECT_EQ(config.BgzfConfig.CompressionLevel, 6);
    EXPECT_EQ(config.BgzfConfig.BgzfWorkers, 4U);
}

TEST_F(BamWriterTest, WriteHeaderOnly)
{
    const SamHeader header = MakeMinimalHeader();
    {
        BamWriter writer{tmpPath, header};
        writer.Close();
    }

    BamRawReader reader{tmpPath};
    EXPECT_EQ(reader.Header().Version(), "1.6");
    EXPECT_EQ(std::size(reader.Header().ReferenceSequences()), 1u);
    EXPECT_FALSE(reader.ReadRecord().has_value());
}

TEST_F(BamWriterTest, WriteAndReadBackRecord)
{
    const SamHeader header = MakeMinimalHeader();
    {
        BamWriter writer{tmpPath, header};
        writer.Write(MakeTestRecord());
    }

    BamRawReader reader{tmpPath};
    const auto view = reader.ReadRecord();
    ASSERT_TRUE(view.has_value());
    EXPECT_EQ(view->Name(), "read1");
    EXPECT_EQ(view->Flag(), 0u);
    EXPECT_EQ(view->RefId(), 0);
    EXPECT_EQ(view->Pos(), 100);
    EXPECT_EQ(view->MapQ(), 30u);
    EXPECT_EQ(view->Seq().ToString(), "ACGTACGTAC");
    EXPECT_FALSE(reader.ReadRecord().has_value());
}

TEST_F(BamWriterTest, WriteMultipleRecords)
{
    const SamHeader header = MakeMinimalHeader();
    {
        BamWriter writer{tmpPath, header};
        for (int i{0}; i < 10; ++i) {
            BamRecord rec = MakeTestRecord();
            rec.Name(std::format("read{}", i));
            writer.Write(rec);
        }
    }

    BamRawReader reader{tmpPath};
    std::size_t count{0};
    for ([[maybe_unused]] const auto& view : reader.Records()) {
        ++count;
    }
    EXPECT_EQ(count, 10u);
}

TEST_F(BamWriterTest, WrittenBamHasEofMarker)
{
    const SamHeader header = MakeMinimalHeader();
    {
        BamWriter writer{tmpPath, header};
        writer.Write(MakeTestRecord());
    }

    const BgzfReader bgzf{tmpPath};
    EXPECT_TRUE(bgzf.HasEofMarker());
}

TEST_F(BamWriterTest, WriteViewZeroCopy)
{
    // Read views from spec_example.bam, write them to new file, read back
    std::vector<std::string> originalNames;
    SamHeader header;

    {
        BamRawReader reader{tests::DataDir / "spec_example.bam"};
        header = reader.Header();
        BamWriter writer{tmpPath, header};

        for (const auto& view : reader.Records()) {
            originalNames.emplace_back(view.Name());
            writer.Write(view);
        }
    }

    BamRawReader reader2{tmpPath};
    std::vector<std::string> writtenNames;
    for (const auto& view : reader2.Records()) {
        writtenNames.emplace_back(view.Name());
    }

    ASSERT_EQ(std::size(writtenNames), std::size(originalNames));
    for (std::size_t i{0}; i < std::size(originalNames); ++i) {
        EXPECT_EQ(writtenNames[i], originalNames[i]);
    }
}

TEST_F(BamWriterTest, MixOwnedAndViewWrites)
{
    const SamHeader header = MakeMinimalHeader();
    {
        BamWriter writer{tmpPath, header};

        // Write one owned record
        writer.Write(MakeTestRecord());

        // Write one view from spec_example.bam
        BamRawReader reader{tests::DataDir / "spec_example.bam"};
        const auto view = reader.ReadRecord();
        ASSERT_TRUE(view.has_value());
        writer.Write(*view);
    }

    BamRawReader reader2{tmpPath};
    std::size_t count{0};
    for ([[maybe_unused]] const auto& view : reader2.Records()) {
        ++count;
    }
    EXPECT_EQ(count, 2u);
}

TEST_F(BamWriterTest, WriteBatch)
{
    SamHeader header;
    std::size_t originalCount{0};

    {
        BamRawReader reader{tests::DataDir / "spec_example.bam"};
        header = reader.Header();
        BamWriter writer{tmpPath, header};

        while (const auto batch = reader.ReadBatch(ByteLimit{1024U * 1024U})) {
            writer.WriteBatch(*batch);
            originalCount += batch->RecordCount();
        }
    }

    BamRawReader reader2{tmpPath};
    std::size_t writtenCount{0};
    for ([[maybe_unused]] const auto& view : reader2.Records()) {
        ++writtenCount;
    }
    EXPECT_EQ(writtenCount, originalCount);
}

TEST_F(BamWriterTest, FullRoundTripFromSpecExample)
{
    // Read all records as owned, write as BamRecord, read back, compare fields
    std::vector<BamRecord> originals;
    SamHeader header;

    {
        BamRawReader reader{tests::DataDir / "spec_example.bam"};
        header = reader.Header();
        for (const auto& view : reader.Records()) {
            originals.push_back(view.ToOwned());
        }
    }

    {
        BamWriter writer{tmpPath, header};
        for (const auto& rec : originals) {
            writer.Write(rec);
        }
    }

    BamRawReader reader2{tmpPath};
    std::size_t idx{0};
    for (const auto& view : reader2.Records()) {
        ASSERT_LT(idx, std::size(originals));
        const auto& orig = originals[idx];
        EXPECT_EQ(view.Name(), orig.Name());
        EXPECT_EQ(view.Flag(), orig.Flag());
        EXPECT_EQ(view.RefId(), orig.RefId());
        EXPECT_EQ(view.Pos(), orig.Pos());
        EXPECT_EQ(view.MapQ(), orig.MapQ());
        EXPECT_EQ(view.NextRefId(), orig.NextRefId());
        EXPECT_EQ(view.NextPos(), orig.NextPos());
        EXPECT_EQ(view.Tlen(), orig.Tlen());
        EXPECT_EQ(view.Seq().ToString(), orig.Sequence());
        ++idx;
    }
    EXPECT_EQ(idx, std::size(originals));
}

TEST_F(BamWriterTest, IndexCallbackInvokedPerRecord)
{
    const SamHeader header = MakeMinimalHeader();
    std::vector<std::int64_t> capturedOffsets;

    {
        auto callback = [&capturedOffsets](std::int64_t offset,
                                           std::span<const std::byte> /*rawData*/) {
            capturedOffsets.push_back(offset);
        };

        const BamWriterConfig config{
            .BgzfConfig = {.BgzfWorkers = 4},
        };
        BamWriter writer{tmpPath, header, config, std::move(callback)};
        for (int i{0}; i < 5; ++i) {
            BamRecord rec = MakeTestRecord();
            rec.Name(std::format("read{}", i));
            writer.Write(rec);
        }
    }

    EXPECT_EQ(std::size(capturedOffsets), 5u);
    // Offsets should be monotonically increasing
    for (std::size_t i{1}; i < std::size(capturedOffsets); ++i) {
        EXPECT_GT(capturedOffsets[i], capturedOffsets[i - 1]);
    }
}

TEST_F(BamWriterTest, NoCallbackStillWorks)
{
    const SamHeader header = MakeMinimalHeader();
    {
        BamWriter writer{tmpPath, header};
        writer.Write(MakeTestRecord());
    }

    BamRawReader reader{tmpPath};
    const auto view = reader.ReadRecord();
    ASSERT_TRUE(view.has_value());
    EXPECT_EQ(view->Name(), "read1");
}

TEST_F(BamWriterTest, ConstructWithConfig)
{
    const BamWriterConfig config{
        .BgzfConfig = {.CompressionLevel = 4, .BgzfWorkers = 2},
    };
    BamWriter writer{tmpPath, MakeMinimalHeader(), config};
    writer.Write(MakeTestRecord());
    writer.Close();

    BamRawReader reader{tmpPath};
    EXPECT_EQ(reader.Header().Version(), "1.6");
}

TEST_F(BamWriterTest, MetricsPopulatedAfterWrite)
{
    const BamWriterConfig config{
        .BgzfConfig = {.BgzfWorkers = 2},
    };
    BamWriter writer{tmpPath, MakeMinimalHeader(), config};
    for (int i{0}; i < 20; ++i) {
        BamRecord rec = MakeTestRecord();
        rec.Name(std::format("read{}", i));
        writer.Write(rec);
    }
    writer.Close();

    const WriterMetrics metrics{writer.GetMetrics()};
    EXPECT_EQ(metrics.TotalRecordsWritten, 20U);
    EXPECT_GT(metrics.Bgzf.BlocksWritten, 0U);
}

}  // namespace Samoa
}  // namespace PacBio
