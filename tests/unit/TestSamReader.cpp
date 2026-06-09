#include "TestData.hpp"
#include "TestTempDir.hpp"

#include <pbsamoa/io/SamReader.hpp>

#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {

// --- Core tests using spec_example.sam ---

TEST(SamReader, OpenSpecExample)
{
    const SamReader reader{tests::DataDir / "spec_example.sam"};
    const auto& header{reader.Header()};
    EXPECT_EQ(header.Version(), "1.6");
    EXPECT_EQ(header.SortOrder(), "coordinate");
    EXPECT_EQ(header.NumReferences(), 1);
    EXPECT_EQ(header.ReferenceSequences()[0].Name(), "ref");
    EXPECT_EQ(header.ReferenceSequences()[0].Length(), 45);
}

TEST(SamReader, ReadAllRecords)
{
    SamReader reader{tests::DataDir / "spec_example.sam"};
    std::vector<std::string> names;
    for (const auto& record : reader.Records()) {
        names.emplace_back(record.Name());
    }
    ASSERT_EQ(std::size(names), 6u);
    EXPECT_EQ(names.front(), "r001");
    EXPECT_EQ(names.back(), "r001");
}

TEST(SamReader, RecordFieldsParsedCorrectly)
{
    SamReader reader{tests::DataDir / "spec_example.sam"};
    const auto record{reader.ReadRecord()};
    ASSERT_TRUE(record);
    EXPECT_EQ(record->Name(), "r001");
    EXPECT_EQ(record->Flag(), 99u);
    EXPECT_EQ(record->RefId(), 0);
    EXPECT_EQ(record->Pos(), 6);  // 0-based (SAM POS=7)
    EXPECT_EQ(record->MapQ(), 30u);
    EXPECT_EQ(record->Tlen(), 39);
    EXPECT_EQ(record->Sequence(), "TTAGATAAAGGATACTG");
}

TEST(SamReader, CigarParsed)
{
    SamReader reader{tests::DataDir / "spec_example.sam"};
    const auto record{reader.ReadRecord()};
    ASSERT_TRUE(record);

    // 8M2I4M1D3M = 5 ops
    const auto& cigar = record->Cigar();
    ASSERT_EQ(std::size(cigar), 5u);
    EXPECT_EQ(cigar[0], CigarOp(CigarOpType::M, 8));
    EXPECT_EQ(cigar[1], CigarOp(CigarOpType::I, 2));
    EXPECT_EQ(cigar[2], CigarOp(CigarOpType::M, 4));
    EXPECT_EQ(cigar[3], CigarOp(CigarOpType::D, 1));
    EXPECT_EQ(cigar[4], CigarOp(CigarOpType::M, 3));
}

TEST(SamReader, TagsParsed)
{
    SamReader reader{tests::DataDir / "spec_example.sam"};

    // Skip r001, r002 to reach r003 (3rd record)
    reader.ReadRecord();
    reader.ReadRecord();
    const auto record{reader.ReadRecord()};
    ASSERT_TRUE(record);
    EXPECT_EQ(record->Name(), "r003");

    const TagKey saKey{'S', 'A'};
    EXPECT_TRUE(record->Tags().Contains(saKey));
    const auto* val{record->Tags().Get(saKey)};
    ASSERT_TRUE(val);
    ASSERT_TRUE(std::holds_alternative<std::string>(*val));
    EXPECT_EQ(std::get<std::string>(*val), "ref,29,-,6H5M,17,0;");
}

TEST(SamReader, UnmappedNextRef)
{
    SamReader reader{tests::DataDir / "spec_example.sam"};

    // r002 is the 2nd record, RNEXT=* → nextRefId=-1
    reader.ReadRecord();
    const auto record{reader.ReadRecord()};
    ASSERT_TRUE(record);
    EXPECT_EQ(record->Name(), "r002");
    EXPECT_EQ(record->NextRefId(), -1);
    EXPECT_EQ(record->NextPos(), -1);
}

TEST(SamReader, EqualSignNextRef)
{
    SamReader reader{tests::DataDir / "spec_example.sam"};
    const auto record{reader.ReadRecord()};
    ASSERT_TRUE(record);
    // r001 has RNEXT="=" → same as RefId
    EXPECT_EQ(record->Name(), "r001");
    EXPECT_EQ(record->NextRefId(), record->RefId());
    EXPECT_EQ(record->NextPos(), 36);  // 0-based (SAM PNEXT=37)
}

TEST(SamReader, NmTagParsed)
{
    SamReader reader{tests::DataDir / "spec_example.sam"};

    // Last record: r001 with NM:i:1
    std::optional<BamRecord> lastRecord;
    for (auto rec = reader.ReadRecord(); rec; rec = reader.ReadRecord()) {
        lastRecord = std::move(rec);
    }
    ASSERT_TRUE(lastRecord);
    EXPECT_EQ(lastRecord->Name(), "r001");

    const TagKey nmKey{'N', 'M'};
    EXPECT_TRUE(lastRecord->Tags().Contains(nmKey));
    const auto* val{lastRecord->Tags().Get(nmKey)};
    ASSERT_TRUE(val);
    ASSERT_TRUE(std::holds_alternative<std::int64_t>(*val));
    EXPECT_EQ(std::get<std::int64_t>(*val), 1);
}

TEST(SamReader, ThrowOnNonexistent)
{
    EXPECT_THROW(SamReader{"/nonexistent/file.sam"}, std::runtime_error);
}

// --- Edge case tests using temp files ---

class SamReaderTempFile : public ::testing::Test
{
protected:
    tests::TempDirGuard tempDir_;
    std::filesystem::path tmpDir_;

    void SetUp() override
    {
        tempDir_.Reset("sam_reader");
        tmpDir_ = tempDir_.Path();
    }

    std::filesystem::path WriteTempSam(std::string_view name, std::string_view content)
    {
        const auto path{tmpDir_ / name};
        std::ofstream out{path};
        out << content;
        return path;
    }
};

TEST_F(SamReaderTempFile, StarSequenceAndQuality)
{
    const auto path = WriteTempSam("star.sam",
                                   "@HD\tVN:1.6\n"
                                   "@SQ\tSN:ref\tLN:100\n"
                                   "read1\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n");
    SamReader reader{path};
    const auto record{reader.ReadRecord()};
    ASSERT_TRUE(record);
    EXPECT_TRUE(std::empty(record->Sequence()));
    EXPECT_TRUE(std::empty(record->Qualities()));
    EXPECT_EQ(record->RefId(), -1);
    EXPECT_EQ(record->Pos(), -1);
    EXPECT_TRUE(std::empty(record->Cigar()));
}

TEST_F(SamReaderTempFile, HeaderlessSamUnmappedParses)
{
    // SAM with no header lines: an unmapped record (RNAME '*') still parses with refId -1.
    const auto path =
        WriteTempSam("headerless.sam", "read1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\t*\n");
    SamReader reader{path};

    EXPECT_EQ(reader.Header().NumReferences(), 0);
    const auto record{reader.ReadRecord()};
    ASSERT_TRUE(record);
    EXPECT_EQ(record->Name(), "read1");
    EXPECT_EQ(record->RefId(), -1);
    EXPECT_EQ(record->Sequence(), "ACGTACGTAC");
}

TEST_F(SamReaderTempFile, HeaderlessSamWithMappedReadIsRejected)
{
    // htslib aborts when a non-'*' RNAME appears but the header has no @SQ lines.
    const auto path = WriteTempSam("headerless_mapped.sam",
                                   "read1\t0\tchr1\t100\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n");
    SamReader reader{path};
    EXPECT_THROW((void)reader.ReadRecord(), std::runtime_error);
}

TEST_F(SamReaderTempFile, RejectsInvalidQName)
{
    // SAMv1 §1.4 col 1: QNAME must match [!-?A-~]{1,254}; '@' (0x40) is excluded.
    const auto path = WriteTempSam("bad_qname.sam",
                                   "@HD\tVN:1.6\n"
                                   "@SQ\tSN:ref\tLN:100\n"
                                   "read@1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n");
    SamReader reader{path};
    EXPECT_THROW((void)reader.ReadRecord(), std::runtime_error);
}

TEST_F(SamReaderTempFile, EmptyFile)
{
    const auto path = WriteTempSam("empty.sam", "");
    SamReader reader{path};
    EXPECT_FALSE(reader.ReadRecord());
}

TEST_F(SamReaderTempFile, HeaderOnlyFile)
{
    const auto path = WriteTempSam("header_only.sam",
                                   "@HD\tVN:1.6\n"
                                   "@SQ\tSN:ref\tLN:100\n");
    SamReader reader{path};
    EXPECT_EQ(reader.Header().Version(), "1.6");
    EXPECT_FALSE(reader.ReadRecord());
}

TEST_F(SamReaderTempFile, MultipleTagsParsed)
{
    const auto path =
        WriteTempSam("multitag.sam",
                     "@HD\tVN:1.6\n"
                     "@SQ\tSN:ref\tLN:100\n"
                     "read1\t0\tref\t10\t30\t5M\t*\t0\t0\tACGTA\t*\tNM:i:2\tRG:Z:group1\n");
    SamReader reader{path};
    const auto record{reader.ReadRecord()};
    ASSERT_TRUE(record);

    const TagKey nmKey{'N', 'M'};
    const TagKey rgKey{'R', 'G'};

    EXPECT_TRUE(record->Tags().Contains(nmKey));
    EXPECT_TRUE(record->Tags().Contains(rgKey));
    const auto* nmVal{record->Tags().Get(nmKey)};
    ASSERT_TRUE(nmVal);
    ASSERT_TRUE(std::holds_alternative<std::int64_t>(*nmVal));
    EXPECT_EQ(std::get<std::int64_t>(*nmVal), 2);
    const auto* rgVal{record->Tags().Get(rgKey)};
    ASSERT_TRUE(rgVal);
    ASSERT_TRUE(std::holds_alternative<std::string>(*rgVal));
    EXPECT_EQ(std::get<std::string>(*rgVal), "group1");
}

TEST_F(SamReaderTempFile, QualitiesDecoded)
{
    // '!' = 33, so qual should be 0; 'I' = 73, so qual should be 40
    const auto path = WriteTempSam("qual.sam",
                                   "@HD\tVN:1.6\n"
                                   "@SQ\tSN:ref\tLN:100\n"
                                   "read1\t0\tref\t10\t30\t3M\t*\t0\t0\tACG\t!5I\n");
    SamReader reader{path};
    const auto record{reader.ReadRecord()};
    ASSERT_TRUE(record);
    const auto& quals = record->Qualities();
    ASSERT_EQ(std::size(quals), 3u);
    EXPECT_EQ(quals[0], 0u);   // '!' - 33 = 0
    EXPECT_EQ(quals[1], 20u);  // '5' - 33 = 20
    EXPECT_EQ(quals[2], 40u);  // 'I' - 33 = 40
}

TEST_F(SamReaderTempFile, NegativeTlen)
{
    const auto path = WriteTempSam("neg_tlen.sam",
                                   "@HD\tVN:1.6\n"
                                   "@SQ\tSN:ref\tLN:100\n"
                                   "read1\t147\tref\t37\t30\t9M\t=\t7\t-39\tCCCCCCCCC\t*\n");
    SamReader reader{path};
    const auto record{reader.ReadRecord()};
    ASSERT_TRUE(record);
    EXPECT_EQ(record->Tlen(), -39);
}

TEST_F(SamReaderTempFile, ThrowsOnMapqOutOfRange)
{
    const auto path = WriteTempSam("bad_mapq.sam",
                                   "@HD\tVN:1.6\n"
                                   "@SQ\tSN:ref\tLN:100\n"
                                   "read1\t0\tref\t10\t300\t3M\t*\t0\t0\tACG\t!!!\n");
    SamReader reader{path};
    EXPECT_THROW(reader.ReadRecord(), std::runtime_error);
}

TEST_F(SamReaderTempFile, ThrowsOnTrailingJunkInNumericField)
{
    const auto path = WriteTempSam("bad_flag.sam",
                                   "@HD\tVN:1.6\n"
                                   "@SQ\tSN:ref\tLN:100\n"
                                   "read1\t99x\tref\t10\t30\t3M\t*\t0\t0\tACG\t!!!\n");
    SamReader reader{path};
    EXPECT_THROW(reader.ReadRecord(), std::runtime_error);
}

TEST_F(SamReaderTempFile, ThrowsOnInvalidOptionalTag)
{
    const auto path = WriteTempSam("bad_tag.sam",
                                   "@HD\tVN:1.6\n"
                                   "@SQ\tSN:ref\tLN:100\n"
                                   "read1\t0\tref\t10\t30\t3M\t*\t0\t0\tACG\t!!!\tNM:i:5x\n");
    SamReader reader{path};
    EXPECT_THROW(reader.ReadRecord(), std::runtime_error);
}

TEST_F(SamReaderTempFile, ThrowsOnSeqQualLengthMismatch)
{
    const auto path = WriteTempSam("bad_seq_qual_len.sam",
                                   "@HD\tVN:1.6\n"
                                   "@SQ\tSN:ref\tLN:100\n"
                                   "read1\t0\tref\t10\t30\t4M\t*\t0\t0\tACGT\t!!!\n");
    SamReader reader{path};
    EXPECT_THROW(reader.ReadRecord(), std::runtime_error);
}

TEST_F(SamReaderTempFile, ThrowsOnQualityWithoutSequence)
{
    const auto path = WriteTempSam("bad_seq_missing_qual.sam",
                                   "@HD\tVN:1.6\n"
                                   "@SQ\tSN:ref\tLN:100\n"
                                   "read1\t4\t*\t0\t0\t*\t*\t0\t0\t*\t!!!\n");
    SamReader reader{path};
    EXPECT_THROW(reader.ReadRecord(), std::runtime_error);
}

TEST_F(SamReaderTempFile, ThrowsOnInvalidQualityCharacter)
{
    const auto path = WriteTempSam("bad_qual_char.sam",
                                   "@HD\tVN:1.6\n"
                                   "@SQ\tSN:ref\tLN:100\n"
                                   "read1\t0\tref\t10\t30\t3M\t*\t0\t0\tACG\t! !\n");
    SamReader reader{path};
    EXPECT_THROW(reader.ReadRecord(), std::runtime_error);
}

TEST_F(SamReaderTempFile, ThrowsOnUnknownRnameWithHeaderDictionary)
{
    const auto path = WriteTempSam("bad_rname.sam",
                                   "@HD\tVN:1.6\n"
                                   "@SQ\tSN:ref\tLN:100\n"
                                   "read1\t0\tchrX\t10\t30\t3M\t*\t0\t0\tACG\t!!!\n");
    SamReader reader{path};
    EXPECT_THROW(reader.ReadRecord(), std::runtime_error);
}

TEST_F(SamReaderTempFile, ThrowsOnUnknownRnextWithHeaderDictionary)
{
    const auto path = WriteTempSam("bad_rnext.sam",
                                   "@HD\tVN:1.6\n"
                                   "@SQ\tSN:ref\tLN:100\n"
                                   "read1\t0\tref\t10\t30\t3M\tchrY\t0\t0\tACG\t!!!\n");
    SamReader reader{path};
    EXPECT_THROW(reader.ReadRecord(), std::runtime_error);
}

TEST_F(SamReaderTempFile, RangeInterfaceEmptyFile)
{
    const auto path = WriteTempSam("empty2.sam", "");
    SamReader reader{path};
    std::size_t count{0};
    for ([[maybe_unused]] const auto& record : reader.Records()) {
        ++count;
    }
    EXPECT_EQ(count, 0u);
}

}  // namespace Samoa
}  // namespace PacBio
