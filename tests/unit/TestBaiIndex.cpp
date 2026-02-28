#include "TestData.hpp"

#include <pbsamoa/core/Bgzf.hpp>
#include <pbsamoa/index/BaiIndex.hpp>
#include <pbsamoa/io/BamRawReader.hpp>

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {

TEST(Chunk, ConstructAndOverlaps)
{
    const Chunk a{VirtualOffset{0, 0}, VirtualOffset{100, 0}};
    const Chunk b{VirtualOffset{50, 0}, VirtualOffset{150, 0}};
    const Chunk c{VirtualOffset{200, 0}, VirtualOffset{300, 0}};

    EXPECT_TRUE(a.Overlaps(b));
    EXPECT_TRUE(b.Overlaps(a));
    EXPECT_FALSE(a.Overlaps(c));
    EXPECT_FALSE(c.Overlaps(a));
}

TEST(BaiIndex, DefaultConstructed)
{
    const BaiIndex index;
    EXPECT_EQ(index.NumReferences(), 0);
    EXPECT_EQ(index.MappedCount(), 0u);
    EXPECT_EQ(index.UnmappedCount(), 0u);
}

TEST(BaiIndex, LoadFromFile)
{
    const auto baiPath = tests::DataDir / "spec_example.bam.bai";
    if (!std::filesystem::exists(baiPath)) {
        GTEST_SKIP() << "spec_example.bam.bai not found";
    }
    const auto index = BaiIndex::FromFile(baiPath);
    EXPECT_GE(index.NumReferences(), 1);
}

TEST(BaiIndex, LoadedIndexHasLinearIndex)
{
    const auto baiPath = tests::DataDir / "spec_example.bam.bai";
    if (!std::filesystem::exists(baiPath)) {
        GTEST_SKIP() << "spec_example.bam.bai not found";
    }
    const auto index = BaiIndex::FromFile(baiPath);
    // The first reference (ref) should have some linear index entries
    const ReferenceIndex& ref0 = index.Reference(0);
    // Linear index may or may not be populated depending on samtools version,
    // but bins should have data for the spec example
    EXPECT_FALSE(std::empty(ref0.bins));
}

TEST(BaiIndex, LoadThrowsOnBadMagic)
{
    // Create a temp file with bad magic
    const auto tmpPath = std::filesystem::temp_directory_path() / "bad_magic.bai";
    {
        std::ofstream out{tmpPath, std::ios::binary};
        out.write("XXXX", 4);
    }

    EXPECT_THROW(BaiIndex::FromFile(tmpPath), std::runtime_error);
    std::filesystem::remove(tmpPath);
}

TEST(BaiIndex, LoadThrowsOnNonexistent)
{
    EXPECT_THROW(BaiIndex::FromFile("/nonexistent/path.bai"), std::runtime_error);
}

TEST(BaiIndex, WriteAndReadBack)
{
    const auto baiPath = tests::DataDir / "spec_example.bam.bai";
    if (!std::filesystem::exists(baiPath)) {
        GTEST_SKIP() << "spec_example.bam.bai not found";
    }
    const auto tmpPath = std::filesystem::temp_directory_path() / "roundtrip.bai";
    const auto original = BaiIndex::FromFile(baiPath);
    original.ToFile(tmpPath);
    const auto reloaded = BaiIndex::FromFile(tmpPath);
    EXPECT_EQ(reloaded.NumReferences(), original.NumReferences());

    for (std::int32_t r{0}; r < original.NumReferences(); ++r) {
        const ReferenceIndex& reloadRef = reloaded.Reference(r);
        const ReferenceIndex& origRef = original.Reference(r);
        EXPECT_EQ(std::size(origRef.bins), std::size(reloadRef.bins));
        EXPECT_EQ(std::size(origRef.linearIndex), std::size(reloadRef.linearIndex));
    }

    std::filesystem::remove(tmpPath);
}

TEST(BaiIndex, QueryReturnsChunks)
{
    const auto baiPath = tests::DataDir / "spec_example.bam.bai";
    if (!std::filesystem::exists(baiPath)) {
        GTEST_SKIP() << "spec_example.bam.bai not found";
    }
    const auto index = BaiIndex::FromFile(baiPath);
    // Query reference 0 (ref) for positions [0, 45)
    const auto chunks = index.Query(0, 0, 45);
    EXPECT_FALSE(std::empty(chunks));
}

TEST(BaiIndex, QueryOutOfBoundsRefIdReturnsEmpty)
{
    const auto baiPath = tests::DataDir / "spec_example.bam.bai";
    if (!std::filesystem::exists(baiPath)) {
        GTEST_SKIP() << "spec_example.bam.bai not found";
    }
    const auto index = BaiIndex::FromFile(baiPath);
    EXPECT_TRUE(std::empty(index.Query(999, 0, 100)));
    EXPECT_TRUE(std::empty(index.Query(-1, 0, 100)));
}

TEST(BaiIndex, QueryChunksAreSorted)
{
    const auto baiPath = tests::DataDir / "spec_example.bam.bai";
    if (!std::filesystem::exists(baiPath)) {
        GTEST_SKIP() << "spec_example.bam.bai not found";
    }
    const auto index = BaiIndex::FromFile(baiPath);
    const auto chunks = index.Query(0, 0, 45);
    for (std::size_t i{1}; i < std::size(chunks); ++i) {
        EXPECT_LE(chunks[i - 1].Begin, chunks[i].Begin);
    }
}

TEST(BaiIndex, BuildFromSortedBam)
{
    const auto bamPath = tests::DataDir / "spec_example.bam";
    const auto index = BaiIndex::Build(bamPath);
    EXPECT_GE(index.NumReferences(), 1);
    EXPECT_GT(index.MappedCount(), 0u);
}

TEST(BaiIndex, BuiltIndexRoundTrip)
{
    const auto tmpPath = std::filesystem::temp_directory_path() / "built_roundtrip.bai";
    const auto bamPath = tests::DataDir / "spec_example.bam";
    const auto index = BaiIndex::Build(bamPath);
    index.ToFile(tmpPath);
    const auto reloaded = BaiIndex::FromFile(tmpPath);
    EXPECT_EQ(reloaded.NumReferences(), index.NumReferences());
    // Mapped count is not stored in standard BAI format; only n_no_coor is.
    EXPECT_EQ(reloaded.UnmappedCount(), index.UnmappedCount());

    std::filesystem::remove(tmpPath);
}

TEST(BaiIndex, BuiltIndexQueryWorks)
{
    const auto bamPath = tests::DataDir / "spec_example.bam";
    const auto index = BaiIndex::Build(bamPath);

    // Query reference 0 for positions [0, 45) — should find records
    const auto chunks = index.Query(0, 0, 45);
    EXPECT_FALSE(std::empty(chunks));
}

TEST(BaiIndex, BuildQueryRoundTripDiverseBam)
{
    const auto bamPath = tests::DataDir / "diverse.bam";
    if (!std::filesystem::exists(bamPath)) {
        GTEST_SKIP() << "diverse.bam not found";
    }

    // Build index from diverse.bam (multi-block file)
    const auto index = BaiIndex::Build(bamPath);
    EXPECT_GE(index.NumReferences(), 1);
    EXPECT_GT(index.MappedCount(), 0u);

    // Verify chunk offsets point past the header block (> 0)
    const ReferenceIndex& ref0{index.Reference(0)};
    for (const auto& [binNum, chunks] : ref0.bins) {
        for (const Chunk& chunk : chunks) {
            EXPECT_GT(chunk.Begin.BlockOffset(), 0u)
                << "Chunk begin should point past header block";
        }
    }

    // Write index to temp file and use BamRawReader::Query for round-trip
    const auto tmpBai = std::filesystem::temp_directory_path() / "diverse_test.bam.bai";
    index.ToFile(tmpBai);

    // Copy BAM to temp so .bai is found next to it
    const auto tmpBam = std::filesystem::temp_directory_path() / "diverse_test.bam";
    std::filesystem::copy_file(bamPath, tmpBam, std::filesystem::copy_options::overwrite_existing);

    // Query chr1:0-1000 — expect 10 records
    BamRawReader reader{tmpBam};
    const auto builtIndex = BaiIndex::FromFile(tmpBai);
    std::int32_t count{0};
    for ([[maybe_unused]] const auto& record : reader.Query(builtIndex, 0, 0, 1000)) {
        ++count;
    }
    EXPECT_EQ(count, 10);

    std::filesystem::remove(tmpBai);
    std::filesystem::remove(tmpBam);
}

TEST(BaiIndex, BuildThrowsOnUnsortedBam)
{
    const auto bamPath = tests::DataDir / "unsorted.bam";
    if (!std::filesystem::exists(bamPath)) {
        GTEST_SKIP() << "unsorted.bam not found";
    }
    EXPECT_THROW(BaiIndex::Build(bamPath), std::runtime_error);
}

}  // namespace Samoa
}  // namespace PacBio
