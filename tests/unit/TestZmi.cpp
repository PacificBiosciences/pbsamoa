#include "TestTempDir.hpp"

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
#include <bit>
#include <filesystem>
#include <format>
#include <ranges>
#include <span>
#include <thread>
#include <type_traits>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {
namespace {

template <typename T>
void WriteAt(std::span<std::byte> dst, std::size_t offset, T value)
{
    static_assert(std::is_integral_v<T>, "WriteAt requires an integral type");
    using UnsignedT = std::make_unsigned_t<T>;
    std::span<std::byte> field{dst.subspan(offset, sizeof(value))};
    const UnsignedT bits{std::bit_cast<UnsignedT>(value)};
    for (std::size_t i{0}; i < sizeof(value); ++i) {
        field[i] = static_cast<std::byte>((bits >> (8U * i)) & static_cast<UnsignedT>(0xFFU));
    }
}

}  // namespace

class ZmiWriterTest : public ::testing::Test
{
protected:
    tests::TempDirGuard tempDir_;
    std::filesystem::path tmpPath_;

    void SetUp() override
    {
        tempDir_.Reset("zmi_writer");
        tmpPath_ = tempDir_.File("index.zmi");
    }
};

TEST(ZmiWriterConfigTest, Defaults)
{
    const ZmiWriterConfig config{};
    EXPECT_FALSE(config.UseTempFile);
}

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

TEST_F(ZmiWriterTest, UseTempFileAtomicWrite)
{
    const ZmiWriterConfig config{.UseTempFile = true};
    {
        ZmiWriter writer{tmpPath_, config};
        writer.AddRecord(100, 42, 1000);
        EXPECT_FALSE(std::filesystem::exists(tmpPath_));
        writer.Close();
    }
    EXPECT_TRUE(std::filesystem::exists(tmpPath_));

    const ZmwIndex index{ZmwIndex::FromZmi(tmpPath_)};
    EXPECT_EQ(index.NumRecords(), 1u);
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

TEST_F(ZmiWriterTest, UniqueZmwsRemovesNonConsecutiveDuplicates)
{
    {
        ZmiWriter writer{tmpPath_};
        writer.AddRecord(100, 42, 1000);
        writer.AddRecord(100, 99, 2000);
        writer.AddRecord(100, 42, 3000);  // duplicate identity, non-consecutive
        writer.AddRecord(200, 7, 4000);
        writer.AddRecord(100, 99, 5000);  // duplicate identity, non-consecutive
    }

    const ZmwIndex index{ZmwIndex::FromZmi(tmpPath_)};
    const auto unique{index.UniqueZmws()};

    ASSERT_EQ(std::size(unique), 3u);
    // Preserve first-seen file order while removing duplicates.
    EXPECT_EQ(unique[0], (ZmwIdentity{100, 42}));
    EXPECT_EQ(unique[1], (ZmwIdentity{100, 99}));
    EXPECT_EQ(unique[2], (ZmwIdentity{200, 7}));
    EXPECT_EQ(index.NumZmws(), std::size(unique));
}

class ZmiBamWriterTest : public ::testing::Test
{
protected:
    tests::TempDirGuard tempDir_;
    std::filesystem::path tmpBamPath_;

    void SetUp() override
    {
        tempDir_.Reset("zmi_bam_writer");
        tmpBamPath_ = tempDir_.File("reads.bam");
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

TEST_F(ZmiBamWriterTest, WritesReadGroupIdFromRgTag)
{
    SamHeader header = MakeMinimalHeader();
    header.AddReadGroup(ReadGroup{"0000002a"});
    header.AddReadGroup(ReadGroup{"0000002b"});

    {
        ZmiBamWriter writer{tmpBamPath_, header};

        BamRecord rec1;
        TagMap tags1;
        tags1.Set(TagKey{'R', 'G'}, std::string{"0000002a"});
        rec1.Name("movie/42/0_100")
            .Flag(0)
            .RefId(0)
            .Pos(0)
            .MapQ(30)
            .Cigar({CigarOp{CigarOpType::M, 10}})
            .Sequence("ACGTACGTAC")
            .Qualities({30, 30, 30, 30, 30, 30, 30, 30, 30, 30})
            .Tags(std::move(tags1));
        writer.Write(rec1);

        BamRecord rec2;
        TagMap tags2;
        tags2.Set(TagKey{'R', 'G'}, std::string{"0000002b/0--0"});
        rec2.Name("movie/42/100_200")
            .Flag(0)
            .RefId(0)
            .Pos(1)
            .MapQ(30)
            .Cigar({CigarOp{CigarOpType::M, 10}})
            .Sequence("ACGTACGTAC")
            .Qualities({30, 30, 30, 30, 30, 30, 30, 30, 30, 30})
            .Tags(std::move(tags2));
        writer.Write(rec2);
    }

    const ZmwIndex index{ZmwIndex::FromZmi(std::filesystem::path{tmpBamPath_.string() + ".zmi"})};
    const auto rg42{index.Find(ZmwIdentity{42, 42})};
    const auto rg43{index.Find(ZmwIdentity{43, 42})};

    ASSERT_EQ(std::size(rg42), 1u);
    ASSERT_EQ(std::size(rg43), 1u);
    EXPECT_NE(rg42[0], rg43[0]);
}

TEST_F(ZmiBamWriterTest, UseTempFileAtomicWrite)
{
    const SamHeader header = MakeMinimalHeader();
    BamWriterConfig config{};
    config.UseTempFile = true;

    {
        ZmiBamWriter writer{tmpBamPath_, header, config};
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
        EXPECT_FALSE(std::filesystem::exists(tmpBamPath_));
        EXPECT_FALSE(std::filesystem::exists(std::filesystem::path{tmpBamPath_.string() + ".zmi"}));
        writer.Close();
    }

    EXPECT_TRUE(std::filesystem::exists(tmpBamPath_));
    EXPECT_TRUE(std::filesystem::exists(std::filesystem::path{tmpBamPath_.string() + ".zmi"}));
}

TEST_F(ZmiWriterTest, FromPbiRoundTrip)
{
    // Construct a minimal synthetic PBI file with 3 records.
    // PBI layout: BGZF([32-byte header][rgId*3][qStart*3][qEnd*3][holeNumber*3][readQual*3][ctxtFlag*3][fileOffset*3])

    const std::filesystem::path pbiPath{tempDir_.File("roundtrip.pbi")};

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
}

TEST_F(ZmiWriterTest, FromPbiThrowsOnIncompleteBasicData)
{
    const std::filesystem::path pbiPath{tempDir_.File("incomplete_basic_data.pbi")};

    {
        std::vector<std::byte> raw(32, std::byte{0});
        raw[0] = std::byte{'P'};
        raw[1] = std::byte{'B'};
        raw[2] = std::byte{'I'};
        raw[3] = std::byte{'\1'};

        const std::uint32_t version{0x030000};
        const std::uint16_t flags{1};
        const std::uint32_t numReads{1};
        WriteAt(raw, 4, version);
        WriteAt(raw, 8, flags);
        WriteAt(raw, 10, numReads);

        BgzfWriter bgzf{pbiPath};
        bgzf.Write(raw);
    }

    EXPECT_THROW(ZmwIndex::FromPbi(pbiPath), std::runtime_error);
}

TEST_F(ZmiWriterTest, FromPbiThrowsWhenBasicDataFlagMissing)
{
    const std::filesystem::path pbiPath{tempDir_.File("missing_basic_data_flags.pbi")};

    {
        std::vector<std::byte> raw(32, std::byte{0});
        raw[0] = std::byte{'P'};
        raw[1] = std::byte{'B'};
        raw[2] = std::byte{'I'};
        raw[3] = std::byte{'\1'};

        const std::uint32_t version{0x030000};
        const std::uint16_t flags{0};  // BasicData missing
        const std::uint32_t numReads{0};
        WriteAt(raw, 4, version);
        WriteAt(raw, 8, flags);
        WriteAt(raw, 10, numReads);

        BgzfWriter bgzf{pbiPath};
        bgzf.Write(raw);
    }

    EXPECT_THROW(ZmwIndex::FromPbi(pbiPath), std::runtime_error);
}

TEST_F(ZmiWriterTest, FromZmiThrowsOnTruncatedEntryData)
{
    {
        std::vector<std::byte> raw(64, std::byte{0});
        raw[0] = std::byte{'Z'};
        raw[1] = std::byte{'M'};
        raw[2] = std::byte{'I'};
        raw[3] = std::byte{'\1'};

        const std::uint32_t version{0x010000};
        const std::uint16_t entrySize{16};
        WriteAt(raw, 4, version);
        WriteAt(raw, 8, entrySize);

        // One complete entry (16 bytes)
        const std::int32_t rgId{100};
        const std::int32_t zmw{42};
        const std::int64_t voffset{1000};
        const std::size_t base{raw.size()};
        raw.resize(base + 16, std::byte{0});
        WriteAt(raw, base, rgId);
        WriteAt(raw, base + 4, zmw);
        WriteAt(raw, base + 8, voffset);

        // Plus a truncated partial entry (4 bytes)
        raw.resize(raw.size() + 4, std::byte{0});
        WriteAt(raw, raw.size() - 4, std::int32_t{200});

        BgzfWriter bgzf{tmpPath_};
        bgzf.Write(raw);
    }

    EXPECT_THROW(ZmwIndex::FromZmi(tmpPath_), std::runtime_error);
}

TEST_F(ZmiWriterTest, FromZmiThrowsOnNumRecordsHeaderMismatch)
{
    {
        std::vector<std::byte> raw(64, std::byte{0});
        raw[0] = std::byte{'Z'};
        raw[1] = std::byte{'M'};
        raw[2] = std::byte{'I'};
        raw[3] = std::byte{'\1'};

        const std::uint32_t version{0x010000};
        const std::uint16_t entrySize{16};
        const std::uint64_t numRecords{2};  // mismatch: we'll write one entry
        WriteAt(raw, 4, version);
        WriteAt(raw, 8, entrySize);
        WriteAt(raw, 12, numRecords);

        const std::int32_t rgId{100};
        const std::int32_t zmw{42};
        const std::int64_t voffset{1000};
        const std::size_t base{raw.size()};
        raw.resize(base + 16, std::byte{0});
        WriteAt(raw, base, rgId);
        WriteAt(raw, base + 4, zmw);
        WriteAt(raw, base + 8, voffset);

        BgzfWriter bgzf{tmpPath_};
        bgzf.Write(raw);
    }

    EXPECT_THROW(ZmwIndex::FromZmi(tmpPath_), std::runtime_error);
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
