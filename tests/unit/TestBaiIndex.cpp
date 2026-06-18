#include "TestData.hpp"
#include "TestTempDir.hpp"

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/Bgzf.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/index/BaiIndex.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamWriter.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <bit>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {

namespace {

void WriteU32LE(std::ofstream& out, std::uint32_t value)
{
    const char bytes[4] = {
        static_cast<char>(value & 0xFFU),
        static_cast<char>((value >> 8U) & 0xFFU),
        static_cast<char>((value >> 16U) & 0xFFU),
        static_cast<char>((value >> 24U) & 0xFFU),
    };
    out.write(bytes, 4);
}

void WriteI32LE(std::ofstream& out, std::int32_t value)
{
    WriteU32LE(out, std::bit_cast<std::uint32_t>(value));
}

SamHeader MakeCoordinateHeader()
{
    SamHeader header;
    header.SetVersion("1.6");
    header.SetSortOrder("coordinate");
    header.AddReferenceSequence(ReferenceSequence{"ref", 1000});
    return header;
}

BamRecord MakeMappedRecord(std::string name, std::int32_t pos)
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

BamRecord MakeUnmappedRecord(std::string name)
{
    BamRecord record;
    record.Name(std::move(name))
        .Flag(0x4)
        .RefId(-1)
        .Pos(-1)
        .MapQ(0)
        .Cigar({})
        .NextRefId(-1)
        .NextPos(-1)
        .Tlen(0)
        .Sequence("ACGTACGTAC")
        .Qualities({30, 30, 30, 30, 30, 30, 30, 30, 30, 30});
    return record;
}

}  // namespace

TEST(BaiIndex, DefaultConstructed)
{
    const BaiIndex index;
    EXPECT_EQ(index.NumReferences(), 0);
    EXPECT_EQ(index.MappedCount(), 0u);
    EXPECT_EQ(index.UnmappedCount(), 0u);
}

TEST(BaiIndex, LoadFromFile)
{
    const auto baiPath = ::tests::DataDir / "spec_example.bam.bai";
    if (!std::filesystem::exists(baiPath)) {
        GTEST_SKIP() << "spec_example.bam.bai not found";
    }
    const auto index = BaiIndex::FromFile(baiPath);
    EXPECT_GE(index.NumReferences(), 1);
}

TEST(BaiIndex, LoadedIndexHasLinearIndex)
{
    const auto baiPath = ::tests::DataDir / "spec_example.bam.bai";
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
    const tests::TempDirGuard tempDir{"bai_bad_magic"};
    const auto tmpPath = tempDir.File("bad_magic.bai");
    {
        std::ofstream out{tmpPath, std::ios::binary};
        out.write("XXXX", 4);
    }

    EXPECT_THROW(BaiIndex::FromFile(tmpPath), std::runtime_error);
}

TEST(BaiIndex, LoadThrowsOnNegativeReferenceCount)
{
    const tests::TempDirGuard tempDir{"bai_negative_reference_count"};
    const auto tmpPath = tempDir.File("negative_n_ref.bai");
    {
        std::ofstream out{tmpPath, std::ios::binary};
        out.write("BAI\1", 4);
        WriteI32LE(out, -1);
    }

    EXPECT_THROW(BaiIndex::FromFile(tmpPath), std::runtime_error);
}

TEST(BaiIndex, LoadThrowsOnNegativeChunkCount)
{
    const tests::TempDirGuard tempDir{"bai_negative_chunk_count"};
    const auto tmpPath = tempDir.File("negative_n_chunks.bai");
    {
        std::ofstream out{tmpPath, std::ios::binary};
        out.write("BAI\1", 4);
        WriteI32LE(out, 1);   // n_ref
        WriteI32LE(out, 1);   // n_bin
        WriteU32LE(out, 0);   // bin
        WriteI32LE(out, -1);  // n_chunks
    }

    EXPECT_THROW(BaiIndex::FromFile(tmpPath), std::runtime_error);
}

TEST(BaiIndex, LoadThrowsOnNegativeLinearIndexCount)
{
    const tests::TempDirGuard tempDir{"bai_negative_linear_index_count"};
    const auto tmpPath = tempDir.File("negative_n_intv.bai");
    {
        std::ofstream out{tmpPath, std::ios::binary};
        out.write("BAI\1", 4);
        WriteI32LE(out, 1);   // n_ref
        WriteI32LE(out, 0);   // n_bin
        WriteI32LE(out, -1);  // n_intv
    }

    EXPECT_THROW(BaiIndex::FromFile(tmpPath), std::runtime_error);
}

TEST(BaiIndex, LoadThrowsOnTruncatedOptionalNoCoorCount)
{
    const tests::TempDirGuard tempDir{"bai_truncated_no_coor"};
    const auto tmpPath = tempDir.File("truncated_n_no_coor.bai");
    {
        std::ofstream out{tmpPath, std::ios::binary};
        out.write("BAI\1", 4);
        WriteI32LE(out, 0);            // n_ref
        WriteU32LE(out, 0x12345678U);  // partial n_no_coor (4/8 bytes)
    }

    EXPECT_THROW(BaiIndex::FromFile(tmpPath), std::runtime_error);
}

TEST(BaiIndex, LoadThrowsOnNonexistent)
{
    EXPECT_THROW(BaiIndex::FromFile("/nonexistent/path.bai"), std::runtime_error);
}

TEST(BaiIndex, WriteAndReadBack)
{
    const auto baiPath = ::tests::DataDir / "spec_example.bam.bai";
    if (!std::filesystem::exists(baiPath)) {
        GTEST_SKIP() << "spec_example.bam.bai not found";
    }
    const tests::TempDirGuard tempDir{"bai_roundtrip"};
    const auto tmpPath = tempDir.File("roundtrip.bai");
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
}

TEST(BaiIndex, QueryReturnsChunks)
{
    const auto baiPath = ::tests::DataDir / "spec_example.bam.bai";
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
    const auto baiPath = ::tests::DataDir / "spec_example.bam.bai";
    if (!std::filesystem::exists(baiPath)) {
        GTEST_SKIP() << "spec_example.bam.bai not found";
    }
    const auto index = BaiIndex::FromFile(baiPath);
    EXPECT_TRUE(std::empty(index.Query(999, 0, 100)));
    EXPECT_TRUE(std::empty(index.Query(-1, 0, 100)));
}

TEST(BaiIndex, QueryInvalidRangesReturnEmpty)
{
    const auto baiPath = ::tests::DataDir / "spec_example.bam.bai";
    if (!std::filesystem::exists(baiPath)) {
        GTEST_SKIP() << "spec_example.bam.bai not found";
    }
    const auto index = BaiIndex::FromFile(baiPath);
    EXPECT_TRUE(std::empty(index.Query(0, 10, 10)));
    EXPECT_TRUE(std::empty(index.Query(0, 10, 5)));
    EXPECT_TRUE(std::empty(index.Query(0, -10, 0)));
    EXPECT_TRUE(std::empty(index.Query(0, -10, -1)));
}

TEST(BaiIndex, QueryChunksAreSorted)
{
    const auto baiPath = ::tests::DataDir / "spec_example.bam.bai";
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
    const auto bamPath = ::tests::DataDir / "spec_example.bam";
    const auto index = BaiIndex::Build(bamPath);
    EXPECT_GE(index.NumReferences(), 1);
    EXPECT_GT(index.MappedCount(), 0u);
}

TEST(BaiIndex, BuiltIndexRoundTrip)
{
    const tests::TempDirGuard tempDir{"bai_built_roundtrip"};
    const auto tmpPath = tempDir.File("built_roundtrip.bai");
    const auto bamPath = ::tests::DataDir / "spec_example.bam";
    const auto index = BaiIndex::Build(bamPath);
    index.ToFile(tmpPath);
    const auto reloaded = BaiIndex::FromFile(tmpPath);
    EXPECT_EQ(reloaded.NumReferences(), index.NumReferences());
    // Mapped count is not stored in standard BAI format; only n_no_coor is.
    EXPECT_EQ(reloaded.UnmappedCount(), index.UnmappedCount());
}

TEST(BaiIndex, BuiltIndexQueryWorks)
{
    const auto bamPath = ::tests::DataDir / "spec_example.bam";
    const auto index = BaiIndex::Build(bamPath);

    // Query reference 0 for positions [0, 45) — should find records
    const auto chunks = index.Query(0, 0, 45);
    EXPECT_FALSE(std::empty(chunks));
}

TEST(BaiIndex, BuildQueryRoundTripDiverseBam)
{
    const auto bamPath = ::tests::DataDir / "diverse.bam";
    if (!std::filesystem::exists(bamPath)) {
        GTEST_SKIP() << "diverse.bam not found";
    }

    // Build index from diverse.bam (multi-block file)
    const auto index = BaiIndex::Build(bamPath);
    EXPECT_GE(index.NumReferences(), 1);
    EXPECT_GT(index.MappedCount(), 0u);

    // Verify chunk offsets point past the header block (> 0). The metadata pseudo-bin
    // (37450) carries statistics, not file offsets, so it is excluded.
    const ReferenceIndex& ref0{index.Reference(0)};
    for (const auto& [binNum, chunks] : ref0.bins) {
        if (binNum == 37450U) {
            continue;
        }
        for (const Chunk& chunk : chunks) {
            EXPECT_GT(chunk.Begin.BlockOffset(), 0u)
                << "Chunk begin should point past header block";
        }
    }

    // Write index to temp file and use BamRawReader::Query for round-trip
    const tests::TempDirGuard tempDir{"bai_diverse_roundtrip"};
    const auto tmpBai = tempDir.File("diverse_test.bam.bai");
    index.ToFile(tmpBai);

    // Copy BAM to temp so .bai is found next to it
    const auto tmpBam = tempDir.File("diverse_test.bam");
    std::filesystem::copy_file(bamPath, tmpBam, std::filesystem::copy_options::overwrite_existing);

    // Query chr1:0-1000 — expect 10 records
    BamRawReader reader{tmpBam};
    const auto builtIndex = BaiIndex::FromFile(tmpBai);
    // The metadata pseudo-bin must round-trip the mapped/unmapped counts (samtools idxstats).
    EXPECT_EQ(builtIndex.MappedCount(), index.MappedCount());
    EXPECT_EQ(builtIndex.UnmappedCount(), index.UnmappedCount());
    std::int32_t count{0};
    for ([[maybe_unused]] const auto& record : reader.Query(builtIndex, 0, 0, 1000)) {
        ++count;
    }
    EXPECT_EQ(count, 10);
}

TEST(BaiIndex, BuildThrowsOnUnsortedBam)
{
    const auto bamPath = ::tests::DataDir / "unsorted.bam";
    if (!std::filesystem::exists(bamPath)) {
        GTEST_SKIP() << "unsorted.bam not found";
    }
    EXPECT_THROW(BaiIndex::Build(bamPath), std::runtime_error);
}

TEST(BaiIndex, BuildThrowsOnPayloadThatDisagreesWithCoordinateHeader)
{
    const tests::TempDirGuard tempDir{"bai_unsorted_payload_coordinate_header"};
    const auto bamPath = tempDir.File("unsorted_payload_coordinate_header.bam");
    {
        BamWriter writer{bamPath, MakeCoordinateHeader()};
        writer.Write(MakeMappedRecord("later", 100));
        writer.Write(MakeMappedRecord("earlier", 50));
    }

    EXPECT_THROW(BaiIndex::Build(bamPath), std::runtime_error);
}

TEST(BaiIndex, BuildThrowsWhenMappedRecordsResumeAfterUnmappedTail)
{
    const tests::TempDirGuard tempDir{"bai_mapped_after_unmapped"};
    const auto bamPath = tempDir.File("mapped_after_unmapped.bam");
    {
        BamWriter writer{bamPath, MakeCoordinateHeader()};
        writer.Write(MakeMappedRecord("mapped-1", 10));
        writer.Write(MakeUnmappedRecord("unmapped"));
        writer.Write(MakeMappedRecord("mapped-2", 20));
    }

    EXPECT_THROW(BaiIndex::Build(bamPath), std::runtime_error);
}

TEST(BaiIndex, BuildParallelMatchesSerialBytes)
{
    // Synthesize a coordinate-sorted BAM large enough to span many BGZF blocks so
    // the parallel path crosses multiple inflate windows (and records that span a
    // window boundary). pos is non-decreasing and clamped within the reference.
    const tests::TempDirGuard tempDir{"bai_parallel_equiv"};
    const auto bamPath = tempDir.File("parallel_equiv.bam");
    {
        BamWriter writer{bamPath, MakeCoordinateHeader()};
        for (std::int32_t i{0}; i < 50000; ++i) {
            const std::int32_t pos{std::min<std::int32_t>(i / 100, 980)};
            writer.Write(MakeMappedRecord("read-" + std::to_string(i), pos));
        }
    }

    const auto readBytes = [](const std::filesystem::path& path) {
        std::ifstream in{path, std::ios::binary | std::ios::ate};
        const std::streamsize size{in.tellg()};
        in.seekg(0);
        std::vector<char> bytes(static_cast<std::size_t>(size));
        in.read(std::data(bytes), size);
        return bytes;
    };

    const BaiIndex serial{BaiIndex::Build(bamPath, 1)};
    const auto serialBai = tempDir.File("serial.bai");
    serial.ToFile(serialBai);
    const std::vector<char> serialBytes{readBytes(serialBai)};

    // The index must be byte-identical regardless of worker count; determinism is
    // the core correctness guarantee of the parallel build.
    for (const std::size_t workers : {std::size_t{2}, std::size_t{3}, std::size_t{8}}) {
        const BaiIndex parallel{BaiIndex::Build(bamPath, workers)};
        const auto parallelBai = tempDir.File("parallel_" + std::to_string(workers) + ".bai");
        parallel.ToFile(parallelBai);

        EXPECT_EQ(readBytes(parallelBai), serialBytes) << "workers=" << workers;
        EXPECT_EQ(parallel.NumReferences(), serial.NumReferences()) << "workers=" << workers;
        EXPECT_EQ(parallel.MappedCount(), serial.MappedCount()) << "workers=" << workers;
        EXPECT_EQ(parallel.UnmappedCount(), serial.UnmappedCount()) << "workers=" << workers;
    }
}

}  // namespace Samoa
}  // namespace PacBio
