#include "TestData.hpp"
#include "TestTempDir.hpp"

#include <pbsamoa/core/GenomicInterval.hpp>
#include <pbsamoa/index/BaiIndex.hpp>
#include <pbsamoa/io/BamFile.hpp>
#include <pbsamoa/io/BamRawReader.hpp>

#include <gtest/gtest.h>

#include <filesystem>
#include <format>

namespace PacBio {
namespace Samoa {

class BamFileTest : public ::testing::Test
{
protected:
    tests::TempDirGuard tempDir_;
    std::filesystem::path tmpBamPath_;
    std::filesystem::path tmpBaiPath_;

    void SetUp() override
    {
        tempDir_.Reset("bam_file");
        tmpBamPath_ = tempDir_.File("input.bam");
        tmpBaiPath_ = std::filesystem::path{tmpBamPath_.string() + ".bai"};
        std::filesystem::copy_file(tests::DataDir / "spec_example.bam", tmpBamPath_,
                                   std::filesystem::copy_options::overwrite_existing);
        std::filesystem::remove(tmpBaiPath_);
    }
};

TEST(BamFile, ThrowsOnMissingFile)
{
    EXPECT_THROW(BamFile{"does_not_exist.bam"}, std::runtime_error);
}

TEST(BamFile, LoadsHeaderAndReferenceMetadata)
{
    const auto bamPath = tests::DataDir / "spec_example.bam";
    BamRawReader reader{bamPath};
    const SamHeader& header = reader.Header();

    const BamFile file{bamPath};
    ASSERT_EQ(file.Header().ToText(), header.ToText());
    ASSERT_EQ(std::size(file.Header().ReferenceSequences()),
              std::size(header.ReferenceSequences()));

    const auto& ref = header.ReferenceSequences().front();
    EXPECT_TRUE(file.HasReference(ref.Name()));
    EXPECT_EQ(file.ReferenceId(ref.Name()), 0);
    EXPECT_EQ(file.ReferenceName(0), ref.Name());
    EXPECT_EQ(file.ReferenceLength(ref.Name()), static_cast<std::uint32_t>(ref.Length()));
    EXPECT_EQ(file.ReferenceLength(0), static_cast<std::uint32_t>(ref.Length()));
    EXPECT_FALSE(file.HasReference("does_not_exist"));
    EXPECT_EQ(file.ReferenceName(99), "");
    EXPECT_EQ(file.ReferenceLength(99), 0u);
}

TEST(BamFile, ReportsEofAndExistingIndexMetadata)
{
    const BamFile file{tests::DataDir / "spec_example.bam"};

    EXPECT_TRUE(file.HasEOF());
    EXPECT_TRUE(file.StandardIndexExists());
    EXPECT_EQ(file.StandardIndexFilename(), tests::DataDir / "spec_example.bam.bai");
}

TEST_F(BamFileTest, CanCreateAndEnsureStandardIndex)
{
    BamFile file{tmpBamPath_};
    EXPECT_FALSE(file.StandardIndexExists());

    file.CreateStandardIndex();
    EXPECT_TRUE(file.StandardIndexExists());
    EXPECT_TRUE(file.StandardIndexIsNewer());
    EXPECT_NO_THROW(BaiIndex::FromFile(tmpBaiPath_));

    std::filesystem::remove(tmpBaiPath_);
    EXPECT_FALSE(file.StandardIndexExists());

    file.EnsureStandardIndexExists();
    EXPECT_TRUE(file.StandardIndexExists());
    EXPECT_NO_THROW(BaiIndex::FromFile(tmpBaiPath_));
}

TEST(GenomicInterval, ValidatesCoordinatesAndMutators)
{
    GenomicInterval interval{"chr1", 10, 20};
    EXPECT_EQ(interval.Name(), "chr1");
    EXPECT_EQ(interval.Start(), 10);
    EXPECT_EQ(interval.Stop(), 20);
    EXPECT_FALSE(interval.Empty());

    EXPECT_THROW(interval.Start(25), std::invalid_argument);
    EXPECT_EQ(interval.Start(), 10);
    EXPECT_EQ(interval.Stop(), 20);

    interval.Name("chr2").Stop(25).Start(25);
    EXPECT_EQ(interval.Name(), "chr2");
    EXPECT_EQ(interval.Start(), 25);
    EXPECT_EQ(interval.Stop(), 25);
    EXPECT_TRUE(interval.Empty());

    EXPECT_THROW(GenomicInterval(std::string{"chr1"}, -1, 5), std::invalid_argument);
    EXPECT_THROW(GenomicInterval(std::string{"chr1"}, 10, 9), std::invalid_argument);
    EXPECT_THROW(interval.Start(-1), std::invalid_argument);
    EXPECT_THROW(interval.Stop(24), std::invalid_argument);
}

}  // namespace Samoa
}  // namespace PacBio
