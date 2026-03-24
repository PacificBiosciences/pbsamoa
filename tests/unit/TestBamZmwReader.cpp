#include "TestTempDir.hpp"

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/io/BamRecordReader.hpp>
#include <pbsamoa/io/BamZmwReader.hpp>
#include <pbsamoa/io/ZmiBamWriter.hpp>

#include <gtest/gtest.h>

#include <chrono>
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
    tests::TempDirGuard tempDir_;
    std::filesystem::path tmpBamPath_;

    void SetUp() override
    {
        tempDir_.Reset("bam_zmw_reader");
        tmpBamPath_ = tempDir_.File("reads.bam");
    }

    static SamHeader MakeMinimalHeader()
    {
        SamHeader header;
        header.SetVersion("1.6");
        header.SetSortOrder("unknown");
        header.AddReferenceSequence(ReferenceSequence{"ref", 10000});
        return header;
    }

    static SamHeader MakeHeaderWithReadGroups()
    {
        SamHeader header{MakeMinimalHeader()};
        header.AddReadGroup(ReadGroup{"0000002a"});
        header.AddReadGroup(ReadGroup{"0000002b"});
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

TEST_F(BamZmwReaderTest, MetricsStartEmptyWithConfiguredCapacity)
{
    WriteBamWithZmws({{42, 2}, {99, 2}});
    BamZmwReader reader{MakeReader(), BamZmwReaderConfig{.PrefetchCapacityZmws = 2}};

    const ZmwReaderMetrics m = reader.GetMetrics();
    EXPECT_EQ(m.ConfiguredCapacity, 2u);
    EXPECT_EQ(m.GroupsConsumed, 0u);
}

TEST_F(BamZmwReaderTest, RejectsZeroPrefetchCapacity)
{
    WriteBamWithZmws({{42, 1}});
    EXPECT_THROW((BamZmwReader{MakeReader(), BamZmwReaderConfig{.PrefetchCapacityZmws = 0}}),
                 std::invalid_argument);
}

TEST_F(BamZmwReaderTest, MetricsShowBufferedGroupsBeforeConsumption)
{
    WriteBamWithZmws({{42, 2}, {99, 3}, {7, 1}});
    BamZmwReader reader{MakeReader(), BamZmwReaderConfig{.PrefetchCapacityZmws = 2}};

    std::this_thread::sleep_for(std::chrono::milliseconds{50});
    const ZmwReaderMetrics before = reader.GetMetrics();
    EXPECT_GT(before.QueueDepth, 0u);
    EXPECT_LE(before.QueueDepth, before.ConfiguredCapacity);
}

TEST_F(BamZmwReaderTest, GetNextPreservesZmwOrderWithPrefetch)
{
    WriteBamWithZmws({{42, 2}, {99, 3}, {7, 1}});
    BamZmwReader reader{MakeReader(), BamZmwReaderConfig{.PrefetchCapacityZmws = 2}};

    std::vector<BamRecord> group;
    ASSERT_TRUE(reader.GetNext(group));
    EXPECT_EQ(reader.CurrentZmw().zmw, 42);

    ASSERT_TRUE(reader.GetNext(group));
    EXPECT_EQ(reader.CurrentZmw().zmw, 99);

    ASSERT_TRUE(reader.GetNext(group));
    EXPECT_EQ(reader.CurrentZmw().zmw, 7);
}

TEST_F(BamZmwReaderTest, QueueDepthNeverExceedsConfiguredCapacity)
{
    WriteBamWithZmws({{42, 2}, {99, 2}, {7, 2}});
    BamZmwReader reader{MakeReader(), BamZmwReaderConfig{.PrefetchCapacityZmws = 1}};

    std::this_thread::sleep_for(std::chrono::milliseconds{50});
    const ZmwReaderMetrics m = reader.GetMetrics();
    EXPECT_LE(m.QueueDepth, 1u);
    EXPECT_LE(m.PeakQueueDepth, 1u);
}

TEST_F(BamZmwReaderTest, MetricsAdvanceAsGroupsAreConsumed)
{
    WriteBamWithZmws({{42, 1}, {99, 1}});
    BamZmwReader reader{MakeReader(), BamZmwReaderConfig{.PrefetchCapacityZmws = 2}};

    std::vector<BamRecord> group;
    ASSERT_TRUE(reader.GetNext(group));
    const ZmwReaderMetrics m = reader.GetMetrics();
    EXPECT_EQ(m.GroupsConsumed, 1u);
}

TEST_F(BamZmwReaderTest, ProducerStallsWhenQueueIsFull)
{
    WriteBamWithZmws({{42, 1}, {99, 1}, {7, 1}});
    BamZmwReader reader{MakeReader(), BamZmwReaderConfig{.PrefetchCapacityZmws = 1}};

    std::this_thread::sleep_for(std::chrono::milliseconds{50});
    const ZmwReaderMetrics m = reader.GetMetrics();
    EXPECT_GT(m.ProducerStalls, 0u);
}

TEST_F(BamZmwReaderTest, EarlyDestructionDoesNotHang)
{
    WriteBamWithZmws({{42, 2}, {99, 2}, {7, 2}, {8, 2}});
    {
        BamZmwReader reader{MakeReader(), BamZmwReaderConfig{.PrefetchCapacityZmws = 1}};
        std::this_thread::sleep_for(std::chrono::milliseconds{10});
    }
    SUCCEED();
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

TEST_F(BamZmwReaderTest, DifferentReadGroupsSameZmwAreSeparateGroups)
{
    // Same hole number, but different RG tags => different ZMW identities.
    const SamHeader header{MakeHeaderWithReadGroups()};
    ZmiBamWriter writer{tmpBamPath_, header};

    // Two records in RG 0x2a with zmw 42
    for (int r{0}; r < 2; ++r) {
        BamRecord rec;
        TagMap tags;
        tags.Set(TagKey{'R', 'G'}, std::string{"0000002a"});
        rec.Name(std::format("movieA/{}/{}_{}", 42, r * 100, (r + 1) * 100))
            .Flag(0)
            .RefId(0)
            .Pos(r)
            .MapQ(30)
            .Cigar({CigarOp{CigarOpType::M, 10}})
            .Sequence("ACGTACGTAC")
            .Qualities({30, 30, 30, 30, 30, 30, 30, 30, 30, 30})
            .Tags(std::move(tags));
        writer.Write(rec);
    }
    // Two records in RG 0x2b with same zmw 42
    for (int r{0}; r < 2; ++r) {
        BamRecord rec;
        TagMap tags;
        tags.Set(TagKey{'R', 'G'}, std::string{"0000002b"});
        rec.Name(std::format("movieA/{}/{}_{}", 42, r * 100, (r + 1) * 100))
            .Flag(0)
            .RefId(0)
            .Pos(100 + r)
            .MapQ(30)
            .Cigar({CigarOp{CigarOpType::M, 10}})
            .Sequence("ACGTACGTAC")
            .Qualities({30, 30, 30, 30, 30, 30, 30, 30, 30, 30})
            .Tags(std::move(tags));
        writer.Write(rec);
    }
    writer.Close();

    BamZmwReader zmwReader{MakeReader()};
    std::vector<BamRecord> group;

    // First group: RG 0x2a, zmw 42 (2 records)
    ASSERT_TRUE(zmwReader.GetNext(group));
    EXPECT_EQ(std::size(group), 2u);
    EXPECT_EQ(zmwReader.CurrentZmw().rgId, 42);
    EXPECT_EQ(zmwReader.CurrentZmw().zmw, 42);
    for (const auto& rec : group) {
        const auto* rg = rec.Tags().Get(TagKey{'R', 'G'});
        ASSERT_NE(rg, nullptr);
        EXPECT_EQ(std::get<std::string>(*rg), "0000002a");
    }

    // Second group: RG 0x2b, same zmw 42 (2 records)
    ASSERT_TRUE(zmwReader.GetNext(group));
    EXPECT_EQ(std::size(group), 2u);
    EXPECT_EQ(zmwReader.CurrentZmw().rgId, 43);
    EXPECT_EQ(zmwReader.CurrentZmw().zmw, 42);
    for (const auto& rec : group) {
        const auto* rg = rec.Tags().Get(TagKey{'R', 'G'});
        ASSERT_NE(rg, nullptr);
        EXPECT_EQ(std::get<std::string>(*rg), "0000002b");
    }

    EXPECT_FALSE(zmwReader.GetNext(group));
}

TEST_F(BamZmwReaderTest, ReadGroupIdParsedFromRgTag)
{
    const SamHeader header{MakeHeaderWithReadGroups()};
    ZmiBamWriter writer{tmpBamPath_, header};

    BamRecord rec;
    TagMap tags;
    tags.Set(TagKey{'R', 'G'}, std::string{"0000002a/0--0"});
    rec.Name("movie/42/0_100")
        .Flag(0)
        .RefId(0)
        .Pos(0)
        .MapQ(30)
        .Cigar({CigarOp{CigarOpType::M, 10}})
        .Sequence("ACGTACGTAC")
        .Qualities({30, 30, 30, 30, 30, 30, 30, 30, 30, 30})
        .Tags(std::move(tags));
    writer.Write(rec);
    writer.Close();

    BamZmwReader zmwReader{MakeReader()};
    std::vector<BamRecord> group;

    ASSERT_TRUE(zmwReader.GetNext(group));
    EXPECT_EQ(zmwReader.CurrentZmw().rgId, 42);
    EXPECT_EQ(zmwReader.CurrentZmw().zmw, 42);
}

}  // namespace Samoa
}  // namespace PacBio
