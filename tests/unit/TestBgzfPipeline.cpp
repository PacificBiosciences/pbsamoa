#include "TestData.hpp"

#include <pbsamoa/core/Bgzf.hpp>

#include <pbsamoa/io/BamRawReader.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <filesystem>
#include <optional>
#include <span>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

TEST(BgzfReaderBam, OpenValidBam)
{
    const std::filesystem::path path{tests::DataDir / "header_only.bam"};
    const BgzfReader pipeline{path, 2};
    EXPECT_TRUE(pipeline.HasEofMarker());
}

TEST(BgzfReaderBam, ThrowOnNonexistentFile)
{
    EXPECT_THROW(BgzfReader("/nonexistent/file.bam", 2), std::runtime_error);
}

TEST(BgzfReaderBam, SyncModeMatchesBgzfReader)
{
    // Compare from file start: BAM-mode reader parses header eagerly.
    const std::filesystem::path path{tests::DataDir / "spec_example.bam"};
    BgzfReader syncReader{path, 0};
    BgzfReader blockReader{path};
    syncReader.Seek(VirtualOffset{0U, 0U});

    std::vector<std::byte> pBuf(65536);
    std::vector<std::byte> rBuf(65536);

    while (true) {
        const auto pResult{syncReader.ReadBlock(pBuf)};
        const auto rResult{blockReader.ReadBlock(rBuf)};

        ASSERT_EQ(pResult.has_value(), rResult.has_value());
        if (!pResult.has_value()) {
            break;
        }
        ASSERT_EQ(*pResult, *rResult);
        if (*pResult == 0) {
            break;
        }
        EXPECT_TRUE(std::ranges::equal(std::span<const std::byte>{pBuf}.first(*pResult),
                                       std::span<const std::byte>{rBuf}.first(*rResult)));
    }
}

TEST(BgzfReaderBam, PipelineModeMatchesSyncMode)
{
    // Pipeline with N workers must produce identical block-by-block output
    const std::filesystem::path path{tests::DataDir / "spec_example.bam"};
    BgzfReader sync{path, 0};
    BgzfReader parallel{path, 4};

    std::vector<std::byte> sBuf(65536);
    std::vector<std::byte> pBuf(65536);
    std::size_t blockIdx{0};

    while (true) {
        const auto sResult{sync.ReadBlock(sBuf)};
        const auto pResult{parallel.ReadBlock(pBuf)};

        ASSERT_EQ(sResult.has_value(), pResult.has_value()) << "Mismatch at block " << blockIdx;
        if (!sResult.has_value()) {
            break;
        }
        ASSERT_EQ(*sResult, *pResult) << "Size mismatch at block " << blockIdx;
        if (*sResult == 0) {
            break;
        }
        EXPECT_TRUE(std::ranges::equal(std::span<const std::byte>{sBuf}.first(*sResult),
                                       std::span<const std::byte>{pBuf}.first(*pResult)))
            << "Content mismatch at block " << blockIdx;
        ++blockIdx;
    }
    EXPECT_GT(blockIdx, 0U);
}

TEST(BgzfReaderBam, ManyRecordsBamCorrectness)
{
    // Larger file with many BGZF blocks
    const std::filesystem::path path{tests::DataDir / "many_records.bam"};
    BgzfReader sync{path, 0};
    BgzfReader parallel{path, 8};

    std::vector<std::byte> sBuf(65536);
    std::vector<std::byte> pBuf(65536);
    std::size_t blockCount{0};

    while (true) {
        const auto sResult{sync.ReadBlock(sBuf)};
        const auto pResult{parallel.ReadBlock(pBuf)};

        ASSERT_EQ(sResult.has_value(), pResult.has_value());
        if (!sResult.has_value() || *sResult == 0) {
            break;
        }
        ASSERT_EQ(*sResult, *pResult);
        EXPECT_TRUE(std::ranges::equal(std::span<const std::byte>{sBuf}.first(*sResult),
                                       std::span<const std::byte>{pBuf}.first(*pResult)));
        ++blockCount;
    }
    EXPECT_GT(blockCount, 1U);
}

TEST(BgzfReaderBam, EofOnSmallFile)
{
    // BAM-mode reader parses header in ctor, so next block is EOF.
    const std::filesystem::path path{tests::DataDir / "header_only.bam"};
    BgzfReader pipeline{path, 2};

    std::vector<std::byte> buf(65536);
    const auto first{pipeline.ReadBlock(buf)};
    ASSERT_TRUE(first.has_value());
    EXPECT_EQ(*first, 0U);
}

TEST(BgzfReaderBam, SeekBackToStart)
{
    const std::filesystem::path path{tests::DataDir / "spec_example.bam"};
    BgzfReader pipeline{path, 4};  // starts after header
    BgzfReader blockReader{path};

    std::vector<std::byte> expectedBuf(65536);
    const auto expectedResult{blockReader.ReadBlock(expectedBuf)};
    ASSERT_TRUE(expectedResult.has_value());
    ASSERT_GT(*expectedResult, 0U);

    std::vector<std::byte> firstBuf(65536);
    const auto firstResult{pipeline.ReadBlock(firstBuf)};
    ASSERT_TRUE(firstResult.has_value());
    ASSERT_GT(*firstResult, 0U);

    // Seek back to start
    pipeline.Seek(VirtualOffset{0U, 0U});

    std::vector<std::byte> secondBuf(65536);
    const auto secondResult{pipeline.ReadBlock(secondBuf)};
    ASSERT_TRUE(secondResult.has_value());
    ASSERT_EQ(*expectedResult, *secondResult);
    EXPECT_TRUE(std::ranges::equal(std::span<const std::byte>{expectedBuf}.first(*expectedResult),
                                   std::span<const std::byte>{secondBuf}.first(*secondResult)));
}

TEST(BgzfReaderBam, TellTracksPosition)
{
    const std::filesystem::path path{tests::DataDir / "spec_example.bam"};
    BgzfReader pipeline{path, 2};

    // Before first read, tell should be 0
    EXPECT_EQ(pipeline.Tell().BlockOffset(), 0U);

    std::vector<std::byte> buf(65536);
    const auto result{pipeline.ReadBlock(buf)};
    ASSERT_TRUE(result.has_value());
    ASSERT_GT(*result, 0U);

    // After reading first block, tell should be at block 0 offset
    // (consumer reports offset of last consumed block)
    const std::uint64_t afterFirst{pipeline.Tell().BlockOffset()};

    // Read another block
    pipeline.ReadBlock(buf);
    // Tell should have advanced (or be at EOF)
    // For this file we just verify it doesn't regress
    EXPECT_GE(pipeline.Tell().BlockOffset(), afterFirst);
}

TEST(BgzfReaderBam, SingleWorkerCorrectness)
{
    const std::filesystem::path path{tests::DataDir / "spec_example.bam"};
    BgzfReader one{path, 1};
    BgzfReader many{path, 8};

    std::vector<std::byte> oneBuf(65536);
    std::vector<std::byte> manyBuf(65536);

    while (true) {
        const auto oneR{one.ReadBlock(oneBuf)};
        const auto manyR{many.ReadBlock(manyBuf)};

        ASSERT_EQ(oneR.has_value(), manyR.has_value());
        if (!oneR.has_value() || *oneR == 0) {
            break;
        }
        ASSERT_EQ(*oneR, *manyR);
        EXPECT_TRUE(std::ranges::equal(std::span<const std::byte>{oneBuf}.first(*oneR),
                                       std::span<const std::byte>{manyBuf}.first(*manyR)));
    }
}

// --- ReadRecord tests (v2 pipeline) ---

TEST(BgzfReaderBam, ReadRecordMatchesSyncBamRawReader)
{
    // Pipeline ReadRecord() must produce identical records to sync BamRawReader
    const std::filesystem::path path{tests::DataDir / "spec_example.bam"};

    // Reference: sync BamRawReader reads all records
    BamRawReader syncReader{path};
    std::vector<std::vector<std::byte>> refRecords;
    while (const auto view = syncReader.ReadRecord()) {
        refRecords.emplace_back(std::ranges::begin(view->RawData()),
                                std::ranges::end(view->RawData()));
    }
    ASSERT_FALSE(refRecords.empty());

    // Pipeline: reads records via new API
    BgzfReader pipeline{path, 4};
    std::size_t idx{0};
    while (const auto rec = pipeline.ReadRecord()) {
        ASSERT_LT(idx, refRecords.size()) << "Pipeline produced more records than sync reader";
        EXPECT_TRUE(std::ranges::equal(rec->RawData(), std::span<const std::byte>{refRecords[idx]}))
            << "Record " << idx << " differs";
        ++idx;
    }
    EXPECT_EQ(idx, refRecords.size()) << "Pipeline produced fewer records";
}

TEST(BgzfReaderBam, ReadRecordManyRecordsBam)
{
    const std::filesystem::path path{tests::DataDir / "many_records.bam"};

    BamRawReader syncReader{path};
    std::vector<std::vector<std::byte>> refRecords;
    while (const auto view = syncReader.ReadRecord()) {
        refRecords.emplace_back(std::ranges::begin(view->RawData()),
                                std::ranges::end(view->RawData()));
    }
    ASSERT_FALSE(refRecords.empty());

    BgzfReader pipeline{path, 4};
    std::size_t idx{0};
    while (const auto rec = pipeline.ReadRecord()) {
        ASSERT_LT(idx, refRecords.size());
        EXPECT_TRUE(
            std::ranges::equal(rec->RawData(), std::span<const std::byte>{refRecords[idx]}));
        ++idx;
    }
    EXPECT_EQ(idx, refRecords.size());
}

TEST(BgzfReaderBam, ReadRecordSingleWorker)
{
    const std::filesystem::path path{tests::DataDir / "spec_example.bam"};

    BamRawReader syncReader{path};
    std::vector<std::vector<std::byte>> refRecords;
    while (const auto view = syncReader.ReadRecord()) {
        refRecords.emplace_back(std::ranges::begin(view->RawData()),
                                std::ranges::end(view->RawData()));
    }
    ASSERT_FALSE(refRecords.empty());

    BgzfReader pipeline{path, 1};
    std::size_t idx{0};
    while (const auto rec = pipeline.ReadRecord()) {
        ASSERT_LT(idx, refRecords.size());
        EXPECT_TRUE(
            std::ranges::equal(rec->RawData(), std::span<const std::byte>{refRecords[idx]}));
        ++idx;
    }
    EXPECT_EQ(idx, refRecords.size());
}

TEST(BgzfReaderBam, ReadRecordEofOnHeaderOnly)
{
    const std::filesystem::path path{tests::DataDir / "header_only.bam"};
    BgzfReader pipeline{path, 2};

    // Should immediately return nullopt — no records in header-only file
    const auto rec{pipeline.ReadRecord()};
    EXPECT_FALSE(rec.has_value());
}

TEST(BgzfReaderBam, ParseHeaderMatchesBamRawReader)
{
    const std::filesystem::path path{tests::DataDir / "spec_example.bam"};

    const BamRawReader syncReader{path};
    const BgzfReader pipeline{path, 4};

    EXPECT_EQ(pipeline.Header().NumReferences(), syncReader.Header().NumReferences());
    EXPECT_EQ(pipeline.Header().Version(), syncReader.Header().Version());
}

}  // namespace Samoa
}  // namespace PacBio
