#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/Bgzf.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/index/ZmwIndex.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/ZmiBamWriter.hpp>
#include <pbsamoa/io/ZmiWriter.hpp>

#include <gtest/gtest.h>

#include <array>
#include <filesystem>
#include <format>
#include <ranges>
#include <span>
#include <thread>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {
namespace {

template <typename T>
void WriteAt(std::span<std::byte> dst, std::size_t offset, T value)
{
    std::span<std::byte> field{dst.subspan(offset, sizeof(value))};
    const std::byte* src{reinterpret_cast<const std::byte*>(&value)};
    std::ranges::copy_n(src, sizeof(value), std::begin(field));
}

}  // namespace

class ZmiWriterTest : public ::testing::Test
{
protected:
    std::filesystem::path tmpPath_;

    void SetUp() override
    {
        const std::size_t threadHash{std::hash<std::thread::id>{}(std::this_thread::get_id())};
        tmpPath_ =
            std::filesystem::temp_directory_path() / std::format("pbsamoa_test_{}.zmi", threadHash);
    }

    void TearDown() override { std::filesystem::remove(tmpPath_); }
};

TEST_F(ZmiWriterTest, WriteEmptyIndex)
{
    {
        const ZmiWriter writer{tmpPath_};
    }
    // File should exist and have at least 64 bytes (header)
    EXPECT_TRUE(std::filesystem::exists(tmpPath_));
    EXPECT_GE(std::filesystem::file_size(tmpPath_), 64u);
}

TEST_F(ZmiWriterTest, WriteEntries)
{
    {
        ZmiWriter writer{tmpPath_};
        writer.AddRecord(100, 42, 1000);
        writer.AddRecord(100, 42, 2000);
        writer.AddRecord(100, 99, 3000);
    }
    EXPECT_TRUE(std::filesystem::exists(tmpPath_));
}

TEST_F(ZmiWriterTest, RoundTrip)
{
    // Write
    {
        ZmiWriter writer{tmpPath_};
        writer.AddRecord(100, 42, 1000);
        writer.AddRecord(100, 42, 2000);
        writer.AddRecord(100, 99, 3000);
        writer.AddRecord(200, 42, 4000);  // same zmw, different rgId
    }

    // Read
    const ZmwIndex index{ZmwIndex::FromZmi(tmpPath_)};
    EXPECT_EQ(index.NumRecords(), 4u);

    // Find by zmw only (returns all rgIds)
    const auto zmw42{index.Find(42)};
    ASSERT_EQ(std::size(zmw42), 3u);  // records at 1000, 2000, 4000
    EXPECT_EQ(zmw42[0], 1000);
    EXPECT_EQ(zmw42[1], 2000);
    EXPECT_EQ(zmw42[2], 4000);

    // Find by (rgId, zmw)
    const auto exact{index.Find(ZmwIdentity{100, 42})};
    ASSERT_EQ(std::size(exact), 2u);
    EXPECT_EQ(exact[0], 1000);
    EXPECT_EQ(exact[1], 2000);

    // Find nonexistent
    const auto missing{index.Find(9999)};
    EXPECT_TRUE(std::empty(missing));
}

TEST_F(ZmiWriterTest, UniqueZmws)
{
    {
        ZmiWriter writer{tmpPath_};
        writer.AddRecord(100, 42, 1000);
        writer.AddRecord(100, 42, 2000);
        writer.AddRecord(100, 99, 3000);
        writer.AddRecord(200, 42, 4000);
    }

    const ZmwIndex index{ZmwIndex::FromZmi(tmpPath_)};
    const auto unique{index.UniqueZmws()};
    ASSERT_EQ(std::size(unique), 3u);

    // File order: (100,42), (100,99), (200,42)
    EXPECT_EQ(unique[0].rgId, 100);
    EXPECT_EQ(unique[0].zmw, 42);
    EXPECT_EQ(unique[1].rgId, 100);
    EXPECT_EQ(unique[1].zmw, 99);
    EXPECT_EQ(unique[2].rgId, 200);
    EXPECT_EQ(unique[2].zmw, 42);

    EXPECT_EQ(index.NumZmws(), 3u);
}

class ZmiBamWriterTest : public ::testing::Test
{
protected:
    std::filesystem::path tmpBamPath_;

    void SetUp() override
    {
        const std::size_t threadHash{std::hash<std::thread::id>{}(std::this_thread::get_id())};
        tmpBamPath_ = std::filesystem::temp_directory_path() /
                      std::format("pbsamoa_zmiwriter_{}.bam", threadHash);
    }

    void TearDown() override
    {
        std::filesystem::remove(tmpBamPath_);
        std::filesystem::remove(std::filesystem::path{tmpBamPath_.string() + ".zmi"});
    }

    static SamHeader MakeMinimalHeader()
    {
        SamHeader header;
        header.SetVersion("1.6");
        header.SetSortOrder("unknown");
        header.AddReferenceSequence(ReferenceSequence{"ref", 1000});
        return header;
    }
};

TEST_F(ZmiBamWriterTest, WritesBamAndZmi)
{
    const SamHeader header = MakeMinimalHeader();
    {
        ZmiBamWriter writer{tmpBamPath_, header};
        BamRecord rec;
        rec.Name("movie/42/0_100")
            .Flag(0)
            .RefId(0)
            .Pos(0)
            .MapQ(30)
            .Cigar({CigarOp{CigarOpType::M, 10}})
            .Sequence("ACGTACGTAC")
            .Qualities({30, 30, 30, 30, 30, 30, 30, 30, 30, 30});
        writer.Write(rec);
    }

    // BAM should be readable
    BamRawReader reader{tmpBamPath_};
    EXPECT_TRUE(reader.ReadRecord().has_value());

    // ZMI should exist alongside BAM
    const std::filesystem::path zmiPath{tmpBamPath_.string() + ".zmi"};
    EXPECT_TRUE(std::filesystem::exists(zmiPath));

    // ZMI should have one record
    const ZmwIndex index{ZmwIndex::FromZmi(zmiPath)};
    EXPECT_EQ(index.NumRecords(), 1u);

    // The zmw should be 42 (parsed from "movie/42/0_100")
    const auto offsets{index.Find(42)};
    ASSERT_EQ(std::size(offsets), 1u);
    EXPECT_GT(offsets[0], 0);  // offset should be positive (after header)
}

TEST_F(ZmiWriterTest, FromPbiRoundTrip)
{
    // Construct a minimal synthetic PBI file with 3 records.
    // PBI layout: BGZF([32-byte header][rgId*3][qStart*3][qEnd*3][holeNumber*3][readQual*3][ctxtFlag*3][fileOffset*3])

    const std::filesystem::path pbiPath{
        std::filesystem::temp_directory_path() /
        std::format("pbsamoa_test_{}.pbi",
                    std::hash<std::thread::id>{}(std::this_thread::get_id()))};

    {
        // Build raw PBI data in a buffer
        std::vector<std::byte> raw;

        // Header (32 bytes)
        raw.resize(32, std::byte{0});
        // Magic "PBI\1"
        raw[0] = std::byte{'P'};
        raw[1] = std::byte{'B'};
        raw[2] = std::byte{'I'};
        raw[3] = std::byte{'\1'};
        // Version at offset 4 (uint32 LE) — 0x030000
        const std::uint32_t version{0x030000};
        WriteAt(raw, 4, version);
        // pbiFlags at offset 8 (uint16 LE) — 0x0001 (BasicData present)
        const std::uint16_t flags{1};
        WriteAt(raw, 8, flags);
        // numReads at offset 10 (uint32 LE) — 3
        const std::uint32_t numReads{3};
        WriteAt(raw, 10, numReads);

        // BasicData columns (3 records each):
        // rgId: [100, 100, 200]
        const std::array<std::int32_t, 3> rgIds{100, 100, 200};
        const std::size_t intColSize{3 * sizeof(std::int32_t)};
        raw.insert(std::end(raw), reinterpret_cast<const std::byte*>(std::data(rgIds)),
                   reinterpret_cast<const std::byte*>(std::data(rgIds)) + intColSize);

        // qStart: [0, 0, 0] (skip)
        const std::array<std::int32_t, 3> qStart{0, 0, 0};
        raw.insert(std::end(raw), reinterpret_cast<const std::byte*>(std::data(qStart)),
                   reinterpret_cast<const std::byte*>(std::data(qStart)) + intColSize);

        // qEnd: [100, 200, 300] (skip)
        const std::array<std::int32_t, 3> qEnd{100, 200, 300};
        raw.insert(std::end(raw), reinterpret_cast<const std::byte*>(std::data(qEnd)),
                   reinterpret_cast<const std::byte*>(std::data(qEnd)) + intColSize);

        // holeNumber: [42, 99, 42]
        const std::array<std::int32_t, 3> holeNumbers{42, 99, 42};
        raw.insert(std::end(raw), reinterpret_cast<const std::byte*>(std::data(holeNumbers)),
                   reinterpret_cast<const std::byte*>(std::data(holeNumbers)) + intColSize);

        // readQual: [0.9f, 0.8f, 0.7f] (skip)
        const std::array<float, 3> readQual{0.9f, 0.8f, 0.7f};
        raw.insert(std::end(raw), reinterpret_cast<const std::byte*>(std::data(readQual)),
                   reinterpret_cast<const std::byte*>(std::data(readQual)) + 3 * sizeof(float));

        // ctxtFlag: [0, 0, 0] (skip)
        const std::array<std::uint8_t, 3> ctxtFlag{0, 0, 0};
        raw.insert(std::end(raw), reinterpret_cast<const std::byte*>(std::data(ctxtFlag)),
                   reinterpret_cast<const std::byte*>(std::data(ctxtFlag)) + 3);

        // fileOffset: [1000, 2000, 3000]
        const std::array<std::int64_t, 3> offsets{1000, 2000, 3000};
        raw.insert(
            std::end(raw), reinterpret_cast<const std::byte*>(std::data(offsets)),
            reinterpret_cast<const std::byte*>(std::data(offsets)) + 3 * sizeof(std::int64_t));

        // Write as BGZF
        BgzfWriter bgzf{pbiPath};
        bgzf.Write(raw);
    }

    // Load via FromPbi
    const ZmwIndex index{ZmwIndex::FromPbi(pbiPath)};
    EXPECT_EQ(index.NumRecords(), 3u);

    // Check Find by zmw
    const auto zmw42{index.Find(42)};
    ASSERT_EQ(std::size(zmw42), 2u);
    EXPECT_EQ(zmw42[0], 1000);
    EXPECT_EQ(zmw42[1], 3000);

    // Check Find by identity
    const auto exact{index.Find(ZmwIdentity{100, 99})};
    ASSERT_EQ(std::size(exact), 1u);
    EXPECT_EQ(exact[0], 2000);

    std::filesystem::remove(pbiPath);
}

TEST_F(ZmiWriterTest, OpenAutoDetectsZmi)
{
    // Write a .zmi file
    {
        ZmiWriter writer{tmpPath_};
        writer.AddRecord(100, 42, 1000);
    }

    // Open() with a "BAM path" that matches (strip .zmi, add .zmi back)
    const std::filesystem::path fakeBamPath{
        tmpPath_.string().substr(0, tmpPath_.string().size() - 4)};  // remove ".zmi"
    const ZmwIndex index{ZmwIndex::Open(fakeBamPath)};
    EXPECT_EQ(index.NumRecords(), 1u);
}

TEST_F(ZmiWriterTest, OpenThrowsWhenNeitherExists)
{
    const std::filesystem::path noSuchBam{"/tmp/nonexistent_pbsamoa_test.bam"};
    EXPECT_THROW(ZmwIndex::Open(noSuchBam), std::runtime_error);
}

TEST_F(ZmiWriterTest, FirstOffsetByZmw)
{
    {
        ZmiWriter writer{tmpPath_};
        writer.AddRecord(100, 42, 1000);
        writer.AddRecord(100, 42, 2000);
        writer.AddRecord(100, 99, 3000);
    }
    const ZmwIndex index{ZmwIndex::FromZmi(tmpPath_)};
    EXPECT_EQ(index.FirstOffset(42), 1000);
    EXPECT_EQ(index.FirstOffset(99), 3000);
}

TEST_F(ZmiWriterTest, FirstOffsetByIdentity)
{
    {
        ZmiWriter writer{tmpPath_};
        writer.AddRecord(100, 42, 1000);
        writer.AddRecord(200, 42, 4000);
    }
    const ZmwIndex index{ZmwIndex::FromZmi(tmpPath_)};
    EXPECT_EQ(index.FirstOffset(ZmwIdentity{100, 42}), 1000);
    EXPECT_EQ(index.FirstOffset(ZmwIdentity{200, 42}), 4000);
}

TEST_F(ZmiWriterTest, FirstOffsetThrowsForMissing)
{
    {
        ZmiWriter writer{tmpPath_};
        writer.AddRecord(100, 42, 1000);
    }
    const ZmwIndex index{ZmwIndex::FromZmi(tmpPath_)};
    EXPECT_THROW(index.FirstOffset(9999), std::runtime_error);
    EXPECT_THROW(index.FirstOffset(ZmwIdentity{999, 9999}), std::runtime_error);
}

TEST_F(ZmiBamWriterTest, WriteRawRecord)
{
    const SamHeader header = MakeMinimalHeader();

    // First write a BamRecord, read back as RawRecord, write to a second file
    const std::filesystem::path tmpBam2{tmpBamPath_.string() + ".copy.bam"};
    {
        ZmiBamWriter writer1{tmpBamPath_, header};
        BamRecord rec;
        rec.Name("movie/42/0_100")
            .Flag(0)
            .RefId(0)
            .Pos(0)
            .MapQ(30)
            .Cigar({CigarOp{CigarOpType::M, 4}})
            .Sequence("ACGT")
            .Qualities({30, 30, 30, 30});
        writer1.Write(rec);
    }

    // Read back as RawRecord and write to second file
    {
        BamRawReader reader{tmpBamPath_};
        ZmiBamWriter writer2{tmpBam2, reader.Header()};
        for (const auto& raw : reader.Records()) {
            writer2.Write(raw);
        }
    }

    // Verify second file has the record
    BamRawReader reader2{tmpBam2};
    const auto rec2{reader2.ReadRecord()};
    ASSERT_TRUE(rec2.has_value());
    EXPECT_EQ(rec2->Name(), "movie/42/0_100");
    EXPECT_FALSE(reader2.ReadRecord().has_value());

    std::filesystem::remove(tmpBam2);
    std::filesystem::remove(std::filesystem::path{tmpBam2.string() + ".zmi"});
}

}  // namespace Samoa
}  // namespace PacBio
