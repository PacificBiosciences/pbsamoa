#include "TestData.hpp"

#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamRecordReader.hpp>

#include <gtest/gtest.h>

#include <filesystem>
#include <string>
#include <vector>

#include <cstddef>

namespace PacBio {
namespace Samoa {

TEST(BamRecordReader, ReadRecordsMatchViewReaderToOwned)
{
    BamRawReader viewReader{tests::DataDir / "spec_example.bam"};
    BamRecordReader recordReader{tests::DataDir / "spec_example.bam",
                                 BamRecordReaderConfig{.DecodeWorkers = 0}};

    std::size_t count{0};
    for (auto viewRec = viewReader.ReadRecord(); viewRec.has_value();
         viewRec = viewReader.ReadRecord()) {
        const BamRecord owned{viewRec->ToOwned()};
        const auto decoded{recordReader.ReadRecord()};
        ASSERT_TRUE(decoded.has_value()) << "at record " << count;
        EXPECT_EQ(decoded->Name(), owned.Name()) << "at record " << count;
        EXPECT_EQ(decoded->Flag(), owned.Flag()) << "at record " << count;
        EXPECT_EQ(decoded->Pos(), owned.Pos()) << "at record " << count;
        EXPECT_EQ(decoded->RefId(), owned.RefId()) << "at record " << count;
        EXPECT_EQ(decoded->MapQ(), owned.MapQ()) << "at record " << count;
        EXPECT_EQ(decoded->Sequence(), owned.Sequence()) << "at record " << count;
        ++count;
    }
    EXPECT_FALSE(recordReader.ReadRecord().has_value());
    EXPECT_GT(count, 0U);
}

TEST(BamRecordReader, HeaderForwarded)
{
    const BamRecordReader reader{tests::DataDir / "spec_example.bam",
                                 BamRecordReaderConfig{.DecodeWorkers = 0}};
    EXPECT_EQ(reader.Header().Version(), "1.6");
}

TEST(BamRecordReader, HeaderOnlyBamProducesNoRecords)
{
    BamRecordReader reader{tests::DataDir / "header_only.bam",
                           BamRecordReaderConfig{.DecodeWorkers = 0}};
    EXPECT_FALSE(reader.ReadRecord().has_value());
}

TEST(BamRecordReader, RangeInterface)
{
    BamRecordReader reader{tests::DataDir / "spec_example.bam",
                           BamRecordReaderConfig{.DecodeWorkers = 0}};
    std::vector<std::string> names;
    for (const auto& rec : reader.Records()) {
        names.emplace_back(rec.Name());
    }
    ASSERT_EQ(std::size(names), 6u);
    EXPECT_EQ(names[0], "r001");
}

TEST(BamRecordReader, ParallelDecodeMatchesSerial)
{
    const auto path = tests::DataDir / "spec_example.bam";
    BamRecordReader serial{path, BamRecordReaderConfig{.DecodeWorkers = 0}};
    BamRecordReader parallel{path, BamRecordReaderConfig{
                                       .RawReaderConfig = {.BgzfWorkers = 2},
                                       .DecodeWorkers = 2,
                                   }};

    std::size_t count{0};
    while (true) {
        const auto s{serial.ReadRecord()};
        const auto p{parallel.ReadRecord()};
        ASSERT_EQ(s.has_value(), p.has_value()) << "at record " << count;
        if (!s.has_value()) {
            break;
        }
        EXPECT_EQ(s->Name(), p->Name()) << "at record " << count;
        EXPECT_EQ(s->Flag(), p->Flag()) << "at record " << count;
        EXPECT_EQ(s->Pos(), p->Pos()) << "at record " << count;
        ++count;
    }
    EXPECT_GT(count, 0U);
}

TEST(BamRecordReader, ManyRecordsParallel)
{
    const auto path = tests::DataDir / "many_records.bam";
    if (!std::filesystem::exists(path)) {
        GTEST_SKIP() << "many_records.bam not generated";
    }

    BamRecordReader reader{path, BamRecordReaderConfig{
                                     .RawReaderConfig = {.BgzfWorkers = 2},
                                     .DecodeWorkers = 2,
                                 }};

    std::size_t count{0};
    for (const auto& rec : reader.Records()) {
        EXPECT_FALSE(std::empty(rec.Name()));
        ++count;
    }
    EXPECT_EQ(count, 500u);
}

TEST(BamRecordReader, TagFilterDropTags)
{
    BamRecordReader noFilter{tests::DataDir / "spec_example.bam",
                             BamRecordReaderConfig{.DecodeWorkers = 0}};
    BamRecordReader withFilter{tests::DataDir / "spec_example.bam",
                               BamRecordReaderConfig{
                                   .DecodeWorkers = 0,
                                   .TagFilter = DropTags{TagKey{'N', 'M'}},
                               }};

    while (true) {
        const auto full{noFilter.ReadRecord()};
        const auto filtered{withFilter.ReadRecord()};
        ASSERT_EQ(full.has_value(), filtered.has_value());
        if (!full.has_value()) {
            break;
        }
        // Core fields match
        EXPECT_EQ(full->Name(), filtered->Name());
        EXPECT_EQ(full->Pos(), filtered->Pos());
        // Filtered record should not have the NM tag
        EXPECT_FALSE(filtered->Tags().Contains(TagKey{'N', 'M'}));
    }
}

TEST(BamRecordReader, EarlyDestruction)
{
    const auto path = tests::DataDir / "many_records.bam";
    if (!std::filesystem::exists(path)) {
        GTEST_SKIP() << "many_records.bam not generated";
    }

    {
        BamRecordReader reader{path, BamRecordReaderConfig{
                                         .RawReaderConfig = {.BgzfWorkers = 2},
                                         .DecodeWorkers = 2,
                                     }};
        // Read just one record, then destroy
        auto rec = reader.ReadRecord();
        EXPECT_TRUE(rec.has_value());
    }
    SUCCEED();
}

}  // namespace Samoa
}  // namespace PacBio
