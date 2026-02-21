#include "TestData.hpp"
#include "TestTempDir.hpp"

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/SamWriter.hpp>

#include <gtest/gtest.h>

#include <filesystem>
#include <format>
#include <fstream>
#include <string>
#include <vector>

namespace PacBio {
namespace Samoa {

class SamWriterTest : public ::testing::Test
{
protected:
    tests::TempDirGuard tempDir_;
    std::filesystem::path tmpPath;

    void SetUp() override
    {
        tempDir_.Reset("sam_writer");
        tmpPath = tempDir_.File("output.sam");
    }

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

    static std::vector<std::string> ReadLines(const std::filesystem::path& path)
    {
        std::vector<std::string> lines;
        std::ifstream in{path};
        std::string line;
        while (std::getline(in, line)) {
            lines.push_back(line);
        }
        return lines;
    }

    /// \brief Return non-header lines from a SAM file.
    static std::vector<std::string> ReadRecordLines(const std::filesystem::path& path)
    {
        std::vector<std::string> result;
        for (const auto& line : ReadLines(path)) {
            if (!line.starts_with('@')) {
                result.push_back(line);
            }
        }
        return result;
    }
};

TEST(SamWriterConfigTest, Defaults)
{
    const SamWriterConfig config{};
    EXPECT_FALSE(config.UseTempFile);
}

TEST_F(SamWriterTest, WriteHeaderOnly)
{
    const SamHeader header = MakeMinimalHeader();
    {
        SamWriter writer{tmpPath, header};
        writer.Close();
    }
    const auto lines = ReadLines(tmpPath);
    ASSERT_GE(std::size(lines), 1u);
    // Should have @HD and @SQ lines
    bool hasHD{false};
    bool hasSQ{false};
    for (const auto& line : lines) {
        if (line.starts_with("@HD")) {
            hasHD = true;
        }
        if (line.starts_with("@SQ")) {
            hasSQ = true;
        }
    }
    EXPECT_TRUE(hasHD);
    EXPECT_TRUE(hasSQ);
}

TEST_F(SamWriterTest, WriteAndReadBack)
{
    const SamHeader header = MakeMinimalHeader();
    {
        SamWriter writer{tmpPath, header};
        writer.Write(MakeTestRecord());
    }

    const auto records{ReadRecordLines(tmpPath)};
    ASSERT_EQ(std::size(records), 1u);

    // Verify the record line contains expected fields (tab-separated)
    const auto& line{records[0]};
    EXPECT_TRUE(line.starts_with("read1\t"));
    EXPECT_NE(line.find("ACGTACGTAC"), std::string::npos);
}

TEST_F(SamWriterTest, WriteMultipleRecords)
{
    const SamHeader header = MakeMinimalHeader();
    {
        SamWriter writer{tmpPath, header};
        for (int i{0}; i < 5; ++i) {
            BamRecord rec = MakeTestRecord();
            rec.Name(std::format("read{}", i));
            writer.Write(rec);
        }
    }

    const auto records{ReadRecordLines(tmpPath)};
    EXPECT_EQ(std::size(records), 5u);
}

TEST_F(SamWriterTest, WriteViewFromBam)
{
    BamRawReader bamReader{tests::DataDir / "spec_example.bam"};
    const SamHeader& header = bamReader.Header();

    std::vector<std::string> originalNames;
    {
        SamWriter writer{tmpPath, header};
        for (const auto& view : bamReader.Records()) {
            originalNames.emplace_back(view.Name());
            writer.Write(view);
        }
    }

    const auto records{ReadRecordLines(tmpPath)};
    ASSERT_EQ(std::size(records), std::size(originalNames));
    for (std::size_t i{0}; i < std::size(originalNames); ++i) {
        EXPECT_TRUE(records[i].starts_with(originalNames[i] + '\t'));
    }
}

TEST_F(SamWriterTest, WriteBatch)
{
    SamHeader header;
    std::size_t originalCount{0};

    {
        BamRawReader reader{tests::DataDir / "spec_example.bam"};
        header = reader.Header();
        SamWriter writer{tmpPath, header};

        while (const auto batch = reader.ReadBatch(ByteLimit{1024U * 1024U})) {
            writer.WriteBatch(*batch);
            originalCount += batch->RecordCount();
        }
    }

    const auto records{ReadRecordLines(tmpPath)};
    EXPECT_EQ(std::size(records), originalCount);
    EXPECT_GT(originalCount, 0u);
}

TEST_F(SamWriterTest, WriteBatchPreservesNames)
{
    BamRawReader reader{tests::DataDir / "spec_example.bam"};
    const SamHeader& header = reader.Header();

    std::vector<std::string> originalNames;
    for (const auto& view : reader.Records()) {
        originalNames.emplace_back(view.Name());
    }

    // Re-read with batch API and write
    BamRawReader reader2{tests::DataDir / "spec_example.bam"};
    {
        SamWriter writer{tmpPath, header};
        while (const auto batch = reader2.ReadBatch(ByteLimit{1024U * 1024U})) {
            writer.WriteBatch(*batch);
        }
    }

    const auto records{ReadRecordLines(tmpPath)};
    ASSERT_EQ(std::size(records), std::size(originalNames));
    for (std::size_t i{0}; i < std::size(originalNames); ++i) {
        EXPECT_TRUE(records[i].starts_with(originalNames[i] + '\t'));
    }
}

TEST_F(SamWriterTest, DoubleCloseIsSafe)
{
    const SamHeader header = MakeMinimalHeader();
    SamWriter writer{tmpPath, header};
    writer.Write(MakeTestRecord());
    writer.Close();
    writer.Close();

    const auto records{ReadRecordLines(tmpPath)};
    EXPECT_EQ(std::size(records), 1u);
}

TEST_F(SamWriterTest, UseTempFileAtomicWrite)
{
    const SamHeader header = MakeMinimalHeader();
    const SamWriterConfig config{.UseTempFile = true};

    {
        SamWriter writer{tmpPath, header, config};
        writer.Write(MakeTestRecord());
        EXPECT_FALSE(std::filesystem::exists(tmpPath));
        writer.Close();
    }

    EXPECT_TRUE(std::filesystem::exists(tmpPath));
    const auto records{ReadRecordLines(tmpPath)};
    ASSERT_EQ(std::size(records), 1u);
    EXPECT_TRUE(records.front().starts_with("read1\t"));
}

}  // namespace Samoa
}  // namespace PacBio
