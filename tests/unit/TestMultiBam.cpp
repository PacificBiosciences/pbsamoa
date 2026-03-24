#include "TestTempDir.hpp"

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/GenomicInterval.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/io/BamCollection.hpp>
#include <pbsamoa/io/BamFile.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamRecordReader.hpp>
#include <pbsamoa/io/BamWriter.hpp>
#include <pbsamoa/io/BamZmwReader.hpp>

#include <gtest/gtest.h>

#include <filesystem>
#include <format>
#include <string>
#include <utility>
#include <vector>

namespace PacBio {
namespace Samoa {
namespace {

SamHeader MakeHeader(std::string rgId = {}, std::string rgSample = {}, std::string pgId = {},
                     std::vector<std::string> comments = {}, std::int32_t refLength = 1000,
                     std::string sortOrder = "coordinate")
{
    SamHeader header;
    header.SetVersion("1.6");
    header.SetSortOrder(std::move(sortOrder));
    header.AddReferenceSequence(ReferenceSequence{"ref", refLength});

    if (!rgId.empty()) {
        ReadGroup rg{std::move(rgId)};
        if (!rgSample.empty()) {
            rg.SetTag("SM", std::move(rgSample));
        }
        header.AddReadGroup(std::move(rg));
    }

    if (!pgId.empty()) {
        ProgramRecord pg{std::move(pgId)};
        header.AddProgramRecord(std::move(pg));
    }

    for (auto& comment : comments) {
        header.AddComment(std::move(comment));
    }

    return header;
}

BamRecord MakeRecord(std::string name, std::int32_t pos)
{
    BamRecord record;
    record.Name(std::move(name))
        .Flag(0)
        .RefId(0)
        .Pos(pos)
        .MapQ(60)
        .Cigar({CigarOp{CigarOpType::M, 10}})
        .NextRefId(-1)
        .NextPos(-1)
        .Tlen(0)
        .Sequence("ACGTACGTAC")
        .Qualities({30, 30, 30, 30, 30, 30, 30, 30, 30, 30});
    return record;
}

void WriteBam(const std::filesystem::path& path, const SamHeader& header,
              std::vector<BamRecord> records)
{
    BamWriter writer{path, header};
    for (auto& record : records) {
        writer.Write(record);
    }
}

}  // namespace

class MultiBamTest : public ::testing::Test
{
protected:
    tests::TempDirGuard tempDir_;
    std::filesystem::path bam1_;
    std::filesystem::path bam2_;
    std::filesystem::path bai1_;
    std::filesystem::path bai2_;

    void SetUp() override
    {
        tempDir_.Reset("multi_bam");
        bam1_ = tempDir_.File("multi_1.bam");
        bam2_ = tempDir_.File("multi_2.bam");
        bai1_ = std::filesystem::path{bam1_.string() + ".bai"};
        bai2_ = std::filesystem::path{bam2_.string() + ".bai"};
    }
};

TEST_F(MultiBamTest, BamCollectionMergesCompatibleHeaders)
{
    WriteBam(bam1_, MakeHeader("rg1", "sample1", "pg1", {"first"}), {MakeRecord("f1_r1", 10)});
    WriteBam(bam2_, MakeHeader("rg2", "sample2", "pg2", {"first", "second"}),
             {MakeRecord("f2_r1", 20)});

    const BamCollection collection{std::vector<std::filesystem::path>{bam1_, bam2_}};
    const SamHeader& header{collection.Header()};

    EXPECT_EQ(collection.Size(), 2u);
    EXPECT_EQ(header.Version(), "1.6");
    EXPECT_EQ(header.SortOrder(), "unknown");
    ASSERT_EQ(std::size(header.ReferenceSequences()), 1u);
    EXPECT_EQ(header.ReferenceSequences()[0].Name(), "ref");
    ASSERT_EQ(std::size(header.ReadGroups()), 2u);
    EXPECT_EQ(header.ReadGroups()[0].Id(), "rg1");
    EXPECT_EQ(header.ReadGroups()[1].Id(), "rg2");
    ASSERT_EQ(std::size(header.ProgramRecords()), 2u);
    EXPECT_EQ(header.ProgramRecords()[0].Id(), "pg1");
    EXPECT_EQ(header.ProgramRecords()[1].Id(), "pg2");
    ASSERT_EQ(std::size(header.Comments()), 2u);
    EXPECT_EQ(header.Comments()[0], "first");
    EXPECT_EQ(header.Comments()[1], "second");
}

TEST_F(MultiBamTest, BamCollectionRejectsConflictingReferenceDictionary)
{
    WriteBam(bam1_, MakeHeader("rg1", "sample1", "pg1", {}, 1000), {MakeRecord("f1_r1", 10)});
    WriteBam(bam2_, MakeHeader("rg2", "sample2", "pg2", {}, 2000), {MakeRecord("f2_r1", 20)});

    EXPECT_THROW((BamCollection{std::vector<std::filesystem::path>{bam1_, bam2_}}),
                 std::runtime_error);
}

TEST_F(MultiBamTest, BamCollectionRejectsConflictingReadGroups)
{
    WriteBam(bam1_, MakeHeader("rg1", "sample1"), {MakeRecord("f1_r1", 10)});
    WriteBam(bam2_, MakeHeader("rg1", "sample2"), {MakeRecord("f2_r1", 20)});

    EXPECT_THROW((BamCollection{std::vector<std::filesystem::path>{bam1_, bam2_}}),
                 std::runtime_error);
}

TEST_F(MultiBamTest, BamRawReaderSupportsSequentialMultiBamInput)
{
    WriteBam(bam1_, MakeHeader("rg1", "sample1"),
             {MakeRecord("f1_r1", 10), MakeRecord("f1_r2", 30)});
    WriteBam(bam2_, MakeHeader("rg2", "sample2"), {MakeRecord("f2_r1", 20)});

    BamRawReader reader{std::vector<std::filesystem::path>{bam1_, bam2_}};
    std::vector<std::string> names;
    for (const auto& record : reader.Records()) {
        names.emplace_back(record.Name());
    }

    EXPECT_EQ((std::vector<std::string>{"f1_r1", "f1_r2", "f2_r1"}), names);
    EXPECT_EQ(std::size(reader.Header().ReadGroups()), 2u);
    EXPECT_EQ(reader.GetMetrics().RecordsConsumed, 3u);
    EXPECT_THROW(reader.Tell(), std::logic_error);
}

TEST_F(MultiBamTest, BamRecordReaderSupportsSequentialAndQueriedMultiBamInput)
{
    WriteBam(bam1_, MakeHeader("rg1", "sample1"),
             {MakeRecord("f1_r1", 10), MakeRecord("f1_r2", 30)});
    WriteBam(bam2_, MakeHeader("rg2", "sample2"),
             {MakeRecord("f2_r1", 20), MakeRecord("f2_r2", 40)});
    BamFile{bam1_}.CreateStandardIndex();
    BamFile{bam2_}.CreateStandardIndex();

    BamRecordReader sequential{std::vector<std::filesystem::path>{bam1_, bam2_},
                               BamRecordReaderConfig{.DecodeWorkers = 0}};
    std::vector<std::string> names;
    for (const auto& record : sequential.Records()) {
        names.emplace_back(record.Name());
    }
    EXPECT_EQ((std::vector<std::string>{"f1_r1", "f1_r2", "f2_r1", "f2_r2"}), names);
    EXPECT_EQ(sequential.GetMetrics().TotalRecordsRead, 4u);

    BamRecordReader queried{std::vector<std::filesystem::path>{bam1_, bam2_},
                            BamRecordReaderConfig{.DecodeWorkers = 0}};
    std::vector<std::int32_t> positions;
    for (const auto& record : queried.Query(GenomicInterval{"ref", 0, 100})) {
        positions.push_back(record.Pos());
    }
    EXPECT_EQ((std::vector<std::int32_t>{10, 20, 30, 40}), positions);
}

TEST_F(MultiBamTest, BamRecordReaderMultiBamQueryThrowsOnMissingIndex)
{
    WriteBam(bam1_, MakeHeader("rg1", "sample1"), {MakeRecord("f1_r1", 10)});
    WriteBam(bam2_, MakeHeader("rg2", "sample2"), {MakeRecord("f2_r1", 20)});
    BamFile{bam1_}.CreateStandardIndex();

    BamRecordReader reader{std::vector<std::filesystem::path>{bam1_, bam2_},
                           BamRecordReaderConfig{.DecodeWorkers = 0}};
    EXPECT_THROW(reader.Query("ref", 0, 100), std::runtime_error);
}

TEST_F(MultiBamTest, BamZmwReaderSupportsMultiBamInput)
{
    WriteBam(bam1_, MakeHeader("rg1", "sample1"), {MakeRecord("movie/42/0_100", 10)});
    WriteBam(bam2_, MakeHeader("rg2", "sample2"),
             {MakeRecord("movie/42/100_200", 20), MakeRecord("movie/99/0_100", 30)});

    BamZmwReader reader{std::vector<std::filesystem::path>{bam1_, bam2_},
                        BamZmwReaderConfig{
                            .Reader = BamRecordReaderConfig{.DecodeWorkers = 0},
                            .PrefetchCapacityZmws = 2,
                        }};
    std::vector<BamRecord> group;

    ASSERT_TRUE(reader.GetNext(group));
    EXPECT_EQ(reader.CurrentZmw().zmw, 42);
    EXPECT_EQ(std::size(group), 2u);

    ASSERT_TRUE(reader.GetNext(group));
    EXPECT_EQ(reader.CurrentZmw().zmw, 99);
    EXPECT_EQ(std::size(group), 1u);

    EXPECT_FALSE(reader.GetNext(group));
}

}  // namespace Samoa
}  // namespace PacBio
