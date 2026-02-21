#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/io/BamRecordReader.hpp>
#include <pbsamoa/io/BamZmwReader.hpp>
#include <pbsamoa/io/ZmiBamWriter.hpp>

#include <gtest/gtest.h>

#include <filesystem>
#include <format>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

class BamZmwReaderTest : public ::testing::Test
{
protected:
    std::filesystem::path tmpBamPath_;

    void SetUp() override
    {
        const std::size_t threadHash{std::hash<std::thread::id>{}(std::this_thread::get_id())};
        tmpBamPath_ = std::filesystem::temp_directory_path() /
                      std::format("pbsamoa_zmwreader_{}.bam", threadHash);
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
        header.AddReferenceSequence(ReferenceSequence{"ref", 10000});
        return header;
    }

    void WriteBamWithZmws(const std::vector<std::pair<std::int32_t, std::int32_t>>& zmwCounts)
    {
        const SamHeader header{MakeMinimalHeader()};
        ZmiBamWriter writer{tmpBamPath_, header};
        for (const auto& [zmw, count] : zmwCounts) {
            for (std::int32_t r{0}; r < count; ++r) {
                BamRecord rec;
                const std::string name{std::format("movie/{}/{}_{}", zmw, r * 100, (r + 1) * 100)};
                rec.Name(name)
                    .Flag(0)
                    .RefId(0)
                    .Pos(zmw * 100 + r)
                    .MapQ(30)
                    .Cigar({CigarOp{CigarOpType::M, 10}})
                    .Sequence("ACGTACGTAC")
                    .Qualities({30, 30, 30, 30, 30, 30, 30, 30, 30, 30});
                writer.Write(rec);
            }
        }
    }

    BamRecordReader MakeReader()
    {
        return BamRecordReader{tmpBamPath_, BamRecordReaderConfig{.DecodeWorkers = 0}};
    }
};

TEST_F(BamZmwReaderTest, GroupsByZmw)
{
    WriteBamWithZmws({{42, 2}, {99, 3}, {7, 1}});
    BamZmwReader zmwReader{MakeReader()};

    std::vector<BamRecord> group;

    ASSERT_TRUE(zmwReader.GetNext(group));
    EXPECT_EQ(zmwReader.CurrentZmw().zmw, 42);
    EXPECT_EQ(std::size(group), 2u);

    ASSERT_TRUE(zmwReader.GetNext(group));
    EXPECT_EQ(zmwReader.CurrentZmw().zmw, 99);
    EXPECT_EQ(std::size(group), 3u);

    ASSERT_TRUE(zmwReader.GetNext(group));
    EXPECT_EQ(zmwReader.CurrentZmw().zmw, 7);
    EXPECT_EQ(std::size(group), 1u);

    EXPECT_FALSE(zmwReader.GetNext(group));
}

TEST_F(BamZmwReaderTest, SingleZmw)
{
    WriteBamWithZmws({{42, 5}});
    BamZmwReader zmwReader{MakeReader()};

    std::vector<BamRecord> group;
    ASSERT_TRUE(zmwReader.GetNext(group));
    EXPECT_EQ(zmwReader.CurrentZmw().zmw, 42);
    EXPECT_EQ(std::size(group), 5u);

    EXPECT_FALSE(zmwReader.GetNext(group));
}

TEST_F(BamZmwReaderTest, EmptyBam)
{
    WriteBamWithZmws({});
    BamZmwReader zmwReader{MakeReader()};

    std::vector<BamRecord> group;
    EXPECT_FALSE(zmwReader.GetNext(group));
}

TEST_F(BamZmwReaderTest, RecordNamesPreserved)
{
    WriteBamWithZmws({{42, 3}});
    BamZmwReader zmwReader{MakeReader()};

    std::vector<BamRecord> group;
    ASSERT_TRUE(zmwReader.GetNext(group));
    ASSERT_EQ(std::size(group), 3u);

    for (const auto& rec : group) {
        const std::string_view name{rec.Name()};
        EXPECT_NE(name.find("movie/42/"), std::string_view::npos);
    }
}

TEST_F(BamZmwReaderTest, RgIdIsZero)
{
    WriteBamWithZmws({{42, 1}, {99, 1}});
    BamZmwReader zmwReader{MakeReader()};

    std::vector<BamRecord> group;

    ASSERT_TRUE(zmwReader.GetNext(group));
    EXPECT_EQ(zmwReader.CurrentZmw().rgId, 0);
    EXPECT_EQ(zmwReader.CurrentZmw().zmw, 42);

    ASSERT_TRUE(zmwReader.GetNext(group));
    EXPECT_EQ(zmwReader.CurrentZmw().rgId, 0);
    EXPECT_EQ(zmwReader.CurrentZmw().zmw, 99);
}

TEST_F(BamZmwReaderTest, HeaderForwarded)
{
    WriteBamWithZmws({{42, 1}});
    const BamZmwReader zmwReader{MakeReader()};
    EXPECT_EQ(zmwReader.Header().Version(), "1.6");
}

TEST_F(BamZmwReaderTest, SplitZmwAcrossCalls)
{
    WriteBamWithZmws({{42, 2}, {99, 1}, {42, 1}});
    BamZmwReader zmwReader{MakeReader()};

    std::vector<BamRecord> group;

    ASSERT_TRUE(zmwReader.GetNext(group));
    EXPECT_EQ(zmwReader.CurrentZmw().zmw, 42);

    ASSERT_TRUE(zmwReader.GetNext(group));
    EXPECT_EQ(zmwReader.CurrentZmw().zmw, 99);

    ASSERT_TRUE(zmwReader.GetNext(group));
    EXPECT_EQ(zmwReader.CurrentZmw().zmw, 42);
    EXPECT_EQ(std::size(group), 1u);
}

TEST_F(BamZmwReaderTest, DifferentMoviesSameZmwAreSeparateGroups)
{
    // Write records from two different movies with the same ZMW hole number.
    // They should be separate groups because the movie name differs.
    const SamHeader header{MakeMinimalHeader()};
    ZmiBamWriter writer{tmpBamPath_, header};

    // Two records from "movieA" with zmw 42
    for (int r{0}; r < 2; ++r) {
        BamRecord rec;
        rec.Name(std::format("movieA/{}/{}_{}", 42, r * 100, (r + 1) * 100))
            .Flag(0)
            .RefId(0)
            .Pos(r)
            .MapQ(30)
            .Cigar({CigarOp{CigarOpType::M, 10}})
            .Sequence("ACGTACGTAC")
            .Qualities({30, 30, 30, 30, 30, 30, 30, 30, 30, 30});
        writer.Write(rec);
    }
    // Two records from "movieB" with same zmw 42
    for (int r{0}; r < 2; ++r) {
        BamRecord rec;
        rec.Name(std::format("movieB/{}/{}_{}", 42, r * 100, (r + 1) * 100))
            .Flag(0)
            .RefId(0)
            .Pos(100 + r)
            .MapQ(30)
            .Cigar({CigarOp{CigarOpType::M, 10}})
            .Sequence("ACGTACGTAC")
            .Qualities({30, 30, 30, 30, 30, 30, 30, 30, 30, 30});
        writer.Write(rec);
    }
    writer.Close();

    BamZmwReader zmwReader{MakeReader()};
    std::vector<BamRecord> group;

    // First group: movieA/42 (2 records)
    ASSERT_TRUE(zmwReader.GetNext(group));
    EXPECT_EQ(std::size(group), 2u);
    EXPECT_EQ(zmwReader.CurrentZmw().zmw, 42);
    for (const auto& rec : group) {
        EXPECT_NE(rec.Name().find("movieA/42/"), std::string_view::npos);
    }

    // Second group: movieB/42 (2 records)
    ASSERT_TRUE(zmwReader.GetNext(group));
    EXPECT_EQ(std::size(group), 2u);
    EXPECT_EQ(zmwReader.CurrentZmw().zmw, 42);
    for (const auto& rec : group) {
        EXPECT_NE(rec.Name().find("movieB/42/"), std::string_view::npos);
    }

    EXPECT_FALSE(zmwReader.GetNext(group));
}

}  // namespace Samoa
}  // namespace PacBio
