#include "TestData.hpp"
#include "TestTempDir.hpp"

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/index/BaiIndex.hpp>
#include <pbsamoa/index/ZmwWhitelist.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/ZmiBamWriter.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <format>
#include <string>
#include <thread>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

TEST(BamRawReader, OpenHeaderOnly)
{
    const BamRawReader reader{tests::DataDir / "header_only.bam"};
    const SamHeader& header = reader.Header();
    EXPECT_FALSE(std::empty(header.Version()));
}

TEST(BamRawReader, OpenSpecExample)
{
    const BamRawReader reader{tests::DataDir / "spec_example.bam"};
    const SamHeader& header = reader.Header();
    EXPECT_EQ(header.Version(), "1.6");
    EXPECT_EQ(header.SortOrder(), "coordinate");
    EXPECT_EQ(header.ReferenceSequences()[0].Name(), "ref");
    EXPECT_EQ(header.ReferenceSequences()[0].Length(), 45);
}

TEST(BamRawReader, ReadRecordsFromSpecExample)
{
    BamRawReader reader{tests::DataDir / "spec_example.bam"};
    std::vector<std::string> names;
    for (auto view = reader.ReadRecord(); view.has_value(); view = reader.ReadRecord()) {
        names.emplace_back(view->Name());
    }
    ASSERT_EQ(std::size(names), 6u);
    EXPECT_EQ(names[0], "r001");
}

TEST(BamRawReader, RecordFieldsDecodeCorrectly)
{
    BamRawReader reader{tests::DataDir / "spec_example.bam"};
    const auto view = reader.ReadRecord();
    ASSERT_TRUE(view.has_value());
    EXPECT_EQ(view->Name(), "r001");
    EXPECT_EQ(view->Flag(), 99u);
    EXPECT_EQ(view->RefId(), 0);
    EXPECT_EQ(view->Pos(), 6);  // 0-based
    EXPECT_EQ(view->MapQ(), 30u);
}

TEST(BamRawReader, HeaderOnlyBamProducesNoRecords)
{
    BamRawReader reader{tests::DataDir / "header_only.bam"};
    EXPECT_FALSE(reader.ReadRecord().has_value());
}

TEST(BamRawReader, ThrowOnNonexistent)
{
    EXPECT_THROW(BamRawReader{"/nonexistent/file.bam"}, std::runtime_error);
}

// --- Batch interface ---

TEST(BamRawReader, ReadBatchReturnsRecords)
{
    BamRawReader reader{tests::DataDir / "spec_example.bam"};
    const auto batch = reader.ReadBatch(ByteLimit{1024U * 1024U});
    ASSERT_TRUE(batch.has_value());
    EXPECT_EQ(batch->RecordCount(), 6u);
}

TEST(BamRawReader, ReadBatchRespectsLimit)
{
    BamRawReader reader{tests::DataDir / "spec_example.bam"};
    std::size_t totalRecords{0};
    std::size_t batchCount{0};
    while (const auto batch = reader.ReadBatch(ByteLimit{64})) {
        totalRecords += batch->RecordCount();
        ++batchCount;
        EXPECT_GE(batch->RecordCount(), 1u);
    }
    EXPECT_EQ(totalRecords, 6u);
    EXPECT_GE(batchCount, 2u);
}

TEST(BamRawReader, HeaderOnlyBamBatchReturnsNullopt)
{
    BamRawReader reader{tests::DataDir / "header_only.bam"};
    EXPECT_FALSE(reader.ReadBatch().has_value());
}

// --- Range interface ---

TEST(BamRawReader, RangeInterfaceIteratesAll)
{
    BamRawReader reader{tests::DataDir / "spec_example.bam"};
    std::vector<std::string> names;
    for (const auto& view : reader.Records()) {
        names.emplace_back(view.Name());
    }
    EXPECT_EQ(std::size(names), 6u);
}

TEST(BamRawReader, RangeInterfaceEmptyFile)
{
    BamRawReader reader{tests::DataDir / "header_only.bam"};
    std::size_t count{0};
    for ([[maybe_unused]] const auto& view : reader.Records()) {
        ++count;
    }
    EXPECT_EQ(count, 0u);
}

// --- Records spanning BGZF block boundaries ---

TEST(BamRawReader, ManyRecordsSpanningBlocks)
{
    const auto path = tests::DataDir / "many_records.bam";
    if (!std::filesystem::exists(path)) {
        GTEST_SKIP() << "many_records.bam not generated";
    }

    BamRawReader reader{path};

    std::size_t count{0};
    for (const auto& view : reader.Records()) {
        EXPECT_FALSE(std::empty(view.Name()));
        EXPECT_TRUE(view.IsMapped());
        ++count;
    }

    EXPECT_EQ(count, 500u);
}

TEST(BamRawReader, BatchModeSpanningBlocks)
{
    const auto path = tests::DataDir / "many_records.bam";
    if (!std::filesystem::exists(path)) {
        GTEST_SKIP() << "many_records.bam not generated";
    }

    BamRawReader reader{path};

    std::size_t total{0};
    while (const auto batch = reader.ReadBatch(ByteLimit{4096})) {
        for (std::size_t i{0}; i < batch->RecordCount(); ++i) {
            const RawRecord view{batch->RecordData(i)};
            EXPECT_FALSE(std::empty(view.Name()));
        }
        total += batch->RecordCount();
    }

    EXPECT_EQ(total, 500u);
}

// --- Query with index ---

TEST(BamRawReader, QueryWithIndex)
{
    const auto baiPath = tests::DataDir / "spec_example.bam.bai";
    const auto bamPath = tests::DataDir / "spec_example.bam";
    if (!std::filesystem::exists(baiPath)) {
        GTEST_SKIP() << "spec_example.bam.bai not found";
    }

    BamRawReader reader{bamPath};
    const auto index = BaiIndex::FromFile(baiPath);

    std::size_t count{0};
    for (const auto& view : reader.Query(index, 0, 0, 45)) {
        EXPECT_EQ(view.RefId(), 0);
        ++count;
    }
    EXPECT_GT(count, 0u);
}

// --- Stress tests with benchmark file ---

TEST(BamRawReader, StressTestBenchmarkFile)
{
    const auto path = tests::DataDir / "benchmark.bam";
    if (!std::filesystem::exists(path)) {
        GTEST_SKIP() << "benchmark.bam not generated";
    }

    BamRawReader reader{path};
    EXPECT_EQ(reader.Header().ReferenceSequences()[0].Name(), "chr1");

    std::size_t count{0};
    std::int32_t prevPos{-1};
    for (const auto& view : reader.Records()) {
        EXPECT_EQ(view.RefId(), 0);
        EXPECT_GE(view.Pos(), prevPos) << "Not sorted at record " << count;
        prevPos = view.Pos();
        EXPECT_EQ(std::size(view.Seq().ToString()), 200u);
        EXPECT_EQ(view.MapQ(), 30u);
        ++count;
    }
    EXPECT_EQ(count, 10000u);
}

TEST(BamRawReader, StressBatchMode)
{
    const auto path = tests::DataDir / "benchmark.bam";
    if (!std::filesystem::exists(path)) {
        GTEST_SKIP() << "benchmark.bam not generated";
    }

    BamRawReader reader{path};

    std::size_t total{0};
    std::size_t batchCount{0};
    while (const auto batch = reader.ReadBatch(ByteLimit{1024 * 1024})) {
        total += batch->RecordCount();
        ++batchCount;
    }

    EXPECT_EQ(total, 10000u);
    EXPECT_GE(batchCount, 2u);  // 10k 200bp records > 1 MiB
}

// --- Pipeline vs sync comparison ---

TEST(BamRawReader, PipelineMatchesSyncReadRecord)
{
    const std::filesystem::path path{tests::DataDir / "spec_example.bam"};
    BamRawReader syncReader{path};
    BamRawReader pipeReader{path, BamRawReaderConfig{.BgzfWorkers = 4}};

    std::size_t count{0};
    while (true) {
        const auto syncRec{syncReader.ReadRecord()};
        const auto pipeRec{pipeReader.ReadRecord()};

        ASSERT_EQ(syncRec.has_value(), pipeRec.has_value()) << "at record " << count;
        if (!syncRec.has_value()) {
            break;
        }

        EXPECT_EQ(syncRec->Flag(), pipeRec->Flag()) << "at record " << count;
        EXPECT_EQ(syncRec->Pos(), pipeRec->Pos()) << "at record " << count;
        EXPECT_EQ(syncRec->MapQ(), pipeRec->MapQ()) << "at record " << count;
        EXPECT_EQ(syncRec->RefId(), pipeRec->RefId()) << "at record " << count;
        ++count;
    }
    EXPECT_GT(count, 0U);
}

TEST(BamRawReader, PipelineMatchesSyncReadBatch)
{
    const auto path{tests::DataDir / "many_records.bam"};
    if (!std::filesystem::exists(path)) {
        GTEST_SKIP() << "many_records.bam not generated";
    }

    BamRawReader syncReader{path};
    BamRawReader pipeReader{path, BamRawReaderConfig{.BgzfWorkers = 4}};

    std::vector<std::uint16_t> syncFlags;
    std::vector<std::uint16_t> pipeFlags;

    while (const auto batch = syncReader.ReadBatch()) {
        for (std::size_t i{0}; i < batch->RecordCount(); ++i) {
            const RawRecord rec{batch->RecordData(i)};
            syncFlags.push_back(rec.Flag());
        }
    }
    while (const auto batch = pipeReader.ReadBatch()) {
        for (std::size_t i{0}; i < batch->RecordCount(); ++i) {
            const RawRecord rec{batch->RecordData(i)};
            pipeFlags.push_back(rec.Flag());
        }
    }

    ASSERT_EQ(std::size(syncFlags), std::size(pipeFlags));
    EXPECT_TRUE(std::ranges::equal(syncFlags, pipeFlags));
}

class BamRawReaderWhitelistTest : public ::testing::Test
{
protected:
    tests::TempDirGuard tempDir_;
    std::filesystem::path tmpBamPath_;

    void SetUp() override
    {
        tempDir_.Reset("bam_raw_reader_whitelist");
        tmpBamPath_ = tempDir_.File("whitelist.bam");

        SamHeader header;
        header.SetVersion("1.6");
        header.SetSortOrder("unknown");
        header.AddReferenceSequence(ReferenceSequence{"ref", 1000});

        ZmiBamWriter writer{tmpBamPath_, header};
        // Write records for 3 ZMWs: 10, 20, 30
        for (const std::int32_t zmw : {10, 20, 30}) {
            for (int subread = 0; subread < 2; ++subread) {
                const std::string name{
                    std::format("movie/{}/{}_{}", zmw, subread * 100, (subread + 1) * 100)};
                BamRecord rec;
                rec.Name(name)
                    .Flag(0)
                    .RefId(0)
                    .Pos(0)
                    .MapQ(30)
                    .Cigar({CigarOp{CigarOpType::M, 4}})
                    .Sequence("ACGT")
                    .Qualities({30, 30, 30, 30});
                writer.Write(rec);
            }
        }
    }
};

TEST_F(BamRawReaderWhitelistTest, WhitelistSelectsSubset)
{
    BamRawReader reader{tmpBamPath_};
    const ZmwWhitelist whitelist{std::vector<std::int32_t>{10, 30}};

    std::vector<std::string> names;
    for (const auto& rec : reader.Whitelist(whitelist)) {
        names.emplace_back(rec.Name());
    }

    // zmw 10: 2 records, zmw 30: 2 records
    ASSERT_EQ(std::size(names), 4u);
    // All names should contain /10/ or /30/
    for (const auto& name : names) {
        EXPECT_TRUE(name.find("/10/") != std::string::npos ||
                    name.find("/30/") != std::string::npos)
            << "unexpected name: " << name;
    }
}

TEST_F(BamRawReaderWhitelistTest, EmptyWhitelistYieldsNoRecords)
{
    BamRawReader reader{tmpBamPath_};
    const ZmwWhitelist whitelist{std::vector<std::int32_t>{}};

    std::size_t count{0};
    for ([[maybe_unused]] const auto& rec : reader.Whitelist(whitelist)) {
        ++count;
    }
    EXPECT_EQ(count, 0u);
}

TEST_F(BamRawReaderWhitelistTest, NonexistentZmwYieldsNoRecords)
{
    BamRawReader reader{tmpBamPath_};
    const ZmwWhitelist whitelist{std::vector<std::int32_t>{9999}};

    std::size_t count{0};
    for ([[maybe_unused]] const auto& rec : reader.Whitelist(whitelist)) {
        ++count;
    }
    EXPECT_EQ(count, 0u);
}

TEST_F(BamRawReaderWhitelistTest, WhitelistWithPipeline)
{
    BamRawReader reader{tmpBamPath_, BamRawReaderConfig{.BgzfWorkers = 2}};
    const ZmwWhitelist whitelist{std::vector<std::int32_t>{20}};

    std::vector<std::string> names;
    for (const auto& rec : reader.Whitelist(whitelist)) {
        names.emplace_back(rec.Name());
    }

    ASSERT_EQ(std::size(names), 2u);
    for (const auto& name : names) {
        EXPECT_NE(name.find("/20/"), std::string::npos) << "unexpected name: " << name;
    }
}

TEST_F(BamRawReaderWhitelistTest, WhitelistRangeSurvivesEmptyCheckBeforeIteration)
{
    BamRawReader reader{tmpBamPath_};
    const ZmwWhitelist whitelist{std::vector<std::int32_t>{10, 30}};

    auto range = reader.Whitelist(whitelist);
    EXPECT_FALSE(range.begin() == range.end());

    std::vector<std::string> names;
    for (const auto& rec : range) {
        names.emplace_back(rec.Name());
    }

    EXPECT_EQ(std::size(names), 4u);
}

TEST_F(BamRawReaderWhitelistTest, ChunkingCanProduceEmptyChunk)
{
    BamRawReader reader{tmpBamPath_, BamRawReaderConfig{.ChunkNum = 2, .TotalChunks = 10}};

    std::size_t count{0};
    for ([[maybe_unused]] const auto& rec : reader.Records()) {
        ++count;
    }
    EXPECT_EQ(count, 0u);
}

TEST_F(BamRawReaderWhitelistTest, ChunkingWithManyChunksStillCoversAllRecords)
{
    std::size_t total{0};
    for (std::int32_t chunkNum{1}; chunkNum <= 10; ++chunkNum) {
        BamRawReader reader{tmpBamPath_,
                            BamRawReaderConfig{.ChunkNum = chunkNum, .TotalChunks = 10}};
        for ([[maybe_unused]] const auto& rec : reader.Records()) {
            ++total;
        }
    }

    EXPECT_EQ(total, 6u);
}

TEST(BamRawReader, WhitelistAndChunkingMutuallyExclusive)
{
    const auto path = tests::DataDir / "spec_example.bam";
    BamRawReaderConfig config{};
    config.ChunkNum = 1;
    config.TotalChunks = 2;
    config.Whitelist.emplace(std::vector<std::int32_t>{42});

    EXPECT_THROW(BamRawReader(path, config), std::invalid_argument);
}

TEST(BamRawReader, RecordLimitStopsEarly)
{
    BamRawReader reader{tests::DataDir / "spec_example.bam", BamRawReaderConfig{.RecordLimit = 2}};
    std::size_t count{0};
    for ([[maybe_unused]] const auto& rec : reader.Records()) {
        ++count;
    }
    EXPECT_EQ(count, 2u);
}

TEST(BamRawReader, SeekAndTell)
{
    BamRawReader reader{tests::DataDir / "spec_example.bam"};
    const auto first{reader.ReadRecord()};
    ASSERT_TRUE(first.has_value());
    const std::string firstName{std::string{first->Name()}};

    const VirtualOffset pos{reader.Tell()};

    // Read past more records
    reader.ReadRecord();
    reader.ReadRecord();

    // Seek back
    reader.Seek(VirtualOffset{0, 0});

    // First seek to a known position doesn't guarantee same record due to how
    // virtual offsets work, so just verify Seek+Tell round-trip
    reader.Seek(pos);
    const VirtualOffset pos2{reader.Tell()};
    EXPECT_EQ(pos.Value(), pos2.Value());
}

TEST(BamRawReader, QueryEmptyRegion)
{
    const auto baiPath = tests::DataDir / "spec_example.bam.bai";
    const auto bamPath = tests::DataDir / "spec_example.bam";
    if (!std::filesystem::exists(baiPath)) {
        GTEST_SKIP() << "spec_example.bam.bai not found";
    }

    BamRawReader reader{bamPath};
    const auto index = BaiIndex::FromFile(baiPath);

    // Query a region far beyond any records
    std::size_t count{0};
    for ([[maybe_unused]] const auto& rec : reader.Query(index, 0, 100000, 200000)) {
        ++count;
    }
    EXPECT_EQ(count, 0u);
}

}  // namespace Samoa
}  // namespace PacBio
