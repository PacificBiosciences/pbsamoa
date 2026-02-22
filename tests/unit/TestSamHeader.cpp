#include "TestData.hpp"

#include <pbsamoa/core/SamHeader.hpp>

#include <pbsamoa/core/Bgzf.hpp>

#include <gtest/gtest.h>

#include <array>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {

TEST(ReferenceSequence, BasicConstruction)
{
    const ReferenceSequence ref{"chr1", 248956422};
    EXPECT_EQ(ref.Name(), "chr1");
    EXPECT_EQ(ref.Length(), 248956422);
    EXPECT_TRUE(std::empty(ref.CustomTags()));
}

TEST(ReferenceSequence, WithOptionalFields)
{
    ReferenceSequence ref{"chrM", 16569};
    ref.SetTag("TP", "circular");
    ref.SetTag("SP", "Homo sapiens");

    EXPECT_EQ(ref.Name(), "chrM");
    EXPECT_EQ(ref.Length(), 16569);
    const std::string* tp = ref.GetTag("TP");
    ASSERT_NE(tp, nullptr);
    EXPECT_EQ(*tp, "circular");
    const std::string* sp = ref.GetTag("SP");
    ASSERT_NE(sp, nullptr);
    EXPECT_EQ(*sp, "Homo sapiens");

    EXPECT_EQ(ref.GetTag("XX"), nullptr);
}

TEST(ReferenceSequence, AlternativeNames)
{
    ReferenceSequence ref{"MT", 16569};
    ref.SetTag("AN", "chrMT,M,chrM");
    const std::string* an = ref.GetTag("AN");
    ASSERT_NE(an, nullptr);
    EXPECT_EQ(*an, "chrMT,M,chrM");
}

TEST(ReadGroup, BasicConstruction)
{
    const ReadGroup rg{"movie1"};
    EXPECT_EQ(rg.Id(), "movie1");
    EXPECT_TRUE(std::empty(rg.CustomTags()));
}

TEST(ReadGroup, WithOptionalFields)
{
    ReadGroup rg{"flowcell1.lane2"};
    rg.SetTag("SM", "sample1");
    rg.SetTag("PL", "ILLUMINA");
    rg.SetTag("LB", "lib1");
    const std::string* sm = rg.GetTag("SM");
    ASSERT_NE(sm, nullptr);
    EXPECT_EQ(*sm, "sample1");
    const std::string* pl = rg.GetTag("PL");
    ASSERT_NE(pl, nullptr);
    EXPECT_EQ(*pl, "ILLUMINA");
}

TEST(ProgramRecord, BasicConstruction)
{
    const ProgramRecord pg{"bwa"};
    EXPECT_EQ(pg.Id(), "bwa");
    EXPECT_TRUE(std::empty(pg.CustomTags()));
}

TEST(ProgramRecord, WithChaining)
{
    ProgramRecord pg{"samtools"};
    pg.SetTag("PN", "samtools");
    pg.SetTag("VN", "1.17");
    pg.SetTag("PP", "bwa");
    pg.SetTag("CL", "samtools sort -o sorted.bam input.bam");
    const std::string* pp = pg.GetTag("PP");
    ASSERT_NE(pp, nullptr);
    EXPECT_EQ(*pp, "bwa");
    const std::string* cl = pg.GetTag("CL");
    ASSERT_NE(cl, nullptr);
    EXPECT_EQ(*cl, "samtools sort -o sorted.bam input.bam");
}

// --- SamHeader parsing tests ---

TEST(SamHeader, ParseMinimalHeader)
{
    const SamHeader header = SamHeader::FromText("@HD\tVN:1.6\n");

    EXPECT_EQ(header.Version(), "1.6");
    EXPECT_TRUE(std::empty(header.ReferenceSequences()));
    EXPECT_TRUE(std::empty(header.ReadGroups()));
    EXPECT_TRUE(std::empty(header.ProgramRecords()));
    EXPECT_TRUE(std::empty(header.Comments()));
}

TEST(SamHeader, ParseSpecExample)
{
    const std::string text = "@HD\tVN:1.6\tSO:coordinate\n"
                             "@SQ\tSN:ref\tLN:45\n";
    const SamHeader header = SamHeader::FromText(text);

    EXPECT_EQ(header.Version(), "1.6");
    EXPECT_EQ(header.SortOrder(), "coordinate");

    ASSERT_EQ(std::size(header.ReferenceSequences()), 1U);
    EXPECT_EQ(header.ReferenceSequences()[0].Name(), "ref");
    EXPECT_EQ(header.ReferenceSequences()[0].Length(), 45);
}

TEST(SamHeader, ParseMultipleReferenceSequences)
{
    const std::string text = "@HD\tVN:1.6\tSO:coordinate\n"
                             "@SQ\tSN:chr1\tLN:248956422\n"
                             "@SQ\tSN:chr2\tLN:242193529\n"
                             "@SQ\tSN:chrM\tLN:16569\tTP:circular\tSP:Homo sapiens\n";
    const SamHeader header = SamHeader::FromText(text);

    ASSERT_EQ(std::size(header.ReferenceSequences()), 3U);
    EXPECT_EQ(header.ReferenceSequences()[0].Name(), "chr1");
    EXPECT_EQ(header.ReferenceSequences()[0].Length(), 248956422);
    EXPECT_EQ(header.ReferenceSequences()[1].Name(), "chr2");
    EXPECT_EQ(header.ReferenceSequences()[1].Length(), 242193529);
    EXPECT_EQ(header.ReferenceSequences()[2].Name(), "chrM");
    EXPECT_EQ(header.ReferenceSequences()[2].Length(), 16569);
    const std::string* tp = header.ReferenceSequences()[2].GetTag("TP");
    ASSERT_NE(tp, nullptr);
    EXPECT_EQ(*tp, "circular");
}

TEST(SamHeader, ParseReadGroups)
{
    const std::string text = "@HD\tVN:1.6\n"
                             "@RG\tID:rg1\tSM:sample1\tPL:ILLUMINA\n"
                             "@RG\tID:rg2\tSM:sample2\tPL:PACBIO\tLB:lib2\n";
    const SamHeader header = SamHeader::FromText(text);

    ASSERT_EQ(std::size(header.ReadGroups()), 2U);
    EXPECT_EQ(header.ReadGroups()[0].Id(), "rg1");
    const std::string* pl = header.ReadGroups()[0].GetTag("PL");
    ASSERT_NE(pl, nullptr);
    EXPECT_EQ(*pl, "ILLUMINA");

    EXPECT_EQ(header.ReadGroups()[1].Id(), "rg2");
    const std::string* lb = header.ReadGroups()[1].GetTag("LB");
    ASSERT_NE(lb, nullptr);
    EXPECT_EQ(*lb, "lib2");
}

TEST(SamHeader, ParseProgramRecords)
{
    const std::string text = "@HD\tVN:1.6\n"
                             "@PG\tID:bwa\tPN:bwa\tVN:0.7.17\tCL:bwa mem ref.fa reads.fq\n"
                             "@PG\tID:samtools\tPN:samtools\tVN:1.17\tPP:bwa\tCL:samtools sort\n";
    const SamHeader header = SamHeader::FromText(text);

    ASSERT_EQ(std::size(header.ProgramRecords()), 2U);
    EXPECT_EQ(header.ProgramRecords()[0].Id(), "bwa");
    const std::string* cl = header.ProgramRecords()[0].GetTag("CL");
    ASSERT_NE(cl, nullptr);
    EXPECT_EQ(*cl, "bwa mem ref.fa reads.fq");

    EXPECT_EQ(header.ProgramRecords()[1].Id(), "samtools");
    const std::string* pp = header.ProgramRecords()[1].GetTag("PP");
    ASSERT_NE(pp, nullptr);
    EXPECT_EQ(*pp, "bwa");
}

TEST(SamHeader, ParseComments)
{
    const std::string text = "@HD\tVN:1.6\n"
                             "@CO\tThis is a comment\n"
                             "@CO\tAnother comment line\n";
    const SamHeader header = SamHeader::FromText(text);

    ASSERT_EQ(std::size(header.Comments()), 2U);
    EXPECT_EQ(header.Comments()[0], "This is a comment");
    EXPECT_EQ(header.Comments()[1], "Another comment line");
}

TEST(SamHeader, ParseNoHdLine)
{
    // Spec says @HD is optional
    const std::string text = "@SQ\tSN:ref\tLN:100\n";
    const SamHeader header = SamHeader::FromText(text);

    EXPECT_TRUE(std::empty(header.Version()));
    ASSERT_EQ(std::size(header.ReferenceSequences()), 1U);
    EXPECT_EQ(header.ReferenceSequences()[0].Name(), "ref");
}

TEST(SamHeader, ParseEmptyHeader)
{
    const SamHeader header = SamHeader::FromText("");

    EXPECT_TRUE(std::empty(header.Version()));
    EXPECT_TRUE(std::empty(header.ReferenceSequences()));
}

TEST(SamHeader, ParseGroupingOrder)
{
    const std::string text = "@HD\tVN:1.6\tGO:query\n";
    const SamHeader header = SamHeader::FromText(text);

    EXPECT_EQ(header.GroupOrder(), "query");
}

TEST(SamHeader, ParseSubSort)
{
    const std::string text = "@HD\tVN:1.6\tSO:coordinate\tSS:coordinate:queryname\n";
    const SamHeader header = SamHeader::FromText(text);

    EXPECT_EQ(header.SortOrder(), "coordinate");
    EXPECT_EQ(header.SubSort(), "coordinate:queryname");
}

TEST(SamHeader, ParseToleratesTrailingNewlines)
{
    const std::string text = "@HD\tVN:1.6\n@SQ\tSN:ref\tLN:45\n\n";
    const SamHeader header = SamHeader::FromText(text);

    EXPECT_EQ(header.Version(), "1.6");
    ASSERT_EQ(std::size(header.ReferenceSequences()), 1U);
}

TEST(SamHeader, ParseToleratesNoTrailingNewline)
{
    const std::string text = "@HD\tVN:1.6\n@SQ\tSN:ref\tLN:45";
    const SamHeader header = SamHeader::FromText(text);

    EXPECT_EQ(header.Version(), "1.6");
    ASSERT_EQ(std::size(header.ReferenceSequences()), 1U);
}

TEST(SamHeader, ParseSqWithAlternativeNames)
{
    const std::string text = "@SQ\tSN:MT\tLN:16569\tAN:chrMT,M,chrM\n";
    const SamHeader header = SamHeader::FromText(text);

    ASSERT_EQ(std::size(header.ReferenceSequences()), 1U);
    const std::string* an = header.ReferenceSequences()[0].GetTag("AN");
    ASSERT_NE(an, nullptr);
    EXPECT_EQ(*an, "chrMT,M,chrM");
}

TEST(SamHeader, RejectsInvalidSqMissingLn)
{
    const std::string text = "@HD\tVN:1.6\n@SQ\tSN:ref\n";
    EXPECT_THROW(SamHeader::FromText(text), std::runtime_error);
}

TEST(SamHeader, RejectsInvalidSqMissingSn)
{
    const std::string text = "@HD\tVN:1.6\n@SQ\tLN:100\n";
    EXPECT_THROW(SamHeader::FromText(text), std::runtime_error);
}

TEST(SamHeader, RejectsInvalidRgMissingId)
{
    const std::string text = "@HD\tVN:1.6\n@RG\tSM:sample\n";
    EXPECT_THROW(SamHeader::FromText(text), std::runtime_error);
}

TEST(SamHeader, RejectsInvalidPgMissingId)
{
    const std::string text = "@HD\tVN:1.6\n@PG\tPN:tool\n";
    EXPECT_THROW(SamHeader::FromText(text), std::runtime_error);
}

// --- SamHeader serialization tests ---

TEST(SamHeader, SerializeMinimalHeader)
{
    SamHeader header;
    header.SetVersion("1.6");
    const std::string text = header.ToText();
    EXPECT_EQ(text, "@HD\tVN:1.6\n");
}

TEST(SamHeader, SerializeWithSortOrder)
{
    SamHeader header;
    header.SetVersion("1.6");
    header.SetSortOrder("coordinate");
    const std::string text = header.ToText();
    EXPECT_EQ(text, "@HD\tVN:1.6\tSO:coordinate\n");
}

TEST(SamHeader, SerializeReferenceSequences)
{
    SamHeader header;
    header.SetVersion("1.6");
    header.AddReferenceSequence(ReferenceSequence{"chr1", 248956422});

    ReferenceSequence chrM{"chrM", 16569};
    chrM.SetTag("TP", "circular");
    header.AddReferenceSequence(std::move(chrM));
    const std::string text = header.ToText();
    const std::string expected = "@HD\tVN:1.6\n"
                                 "@SQ\tSN:chr1\tLN:248956422\n"
                                 "@SQ\tSN:chrM\tLN:16569\tTP:circular\n";
    EXPECT_EQ(text, expected);
}

TEST(SamHeader, SerializeReadGroups)
{
    SamHeader header;
    header.SetVersion("1.6");

    ReadGroup rg{"rg1"};
    rg.SetTag("SM", "sample1");
    rg.SetTag("PL", "ILLUMINA");
    header.AddReadGroup(std::move(rg));
    const std::string text = header.ToText();
    const std::string expected = "@HD\tVN:1.6\n"
                                 "@RG\tID:rg1\tSM:sample1\tPL:ILLUMINA\n";
    EXPECT_EQ(text, expected);
}

TEST(SamHeader, SerializeProgramRecords)
{
    SamHeader header;
    header.SetVersion("1.6");

    ProgramRecord pg{"bwa"};
    pg.SetTag("PN", "bwa");
    pg.SetTag("VN", "0.7.17");
    header.AddProgramRecord(std::move(pg));
    const std::string text = header.ToText();
    const std::string expected = "@HD\tVN:1.6\n"
                                 "@PG\tID:bwa\tPN:bwa\tVN:0.7.17\n";
    EXPECT_EQ(text, expected);
}

TEST(SamHeader, SerializeComments)
{
    SamHeader header;
    header.SetVersion("1.6");
    header.AddComment("This is a comment");
    const std::string text = header.ToText();
    const std::string expected = "@HD\tVN:1.6\n"
                                 "@CO\tThis is a comment\n";
    EXPECT_EQ(text, expected);
}

TEST(SamHeader, SerializeLineOrder)
{
    // Spec order: @HD, @SQ, @RG, @PG, @CO
    SamHeader header;
    header.SetVersion("1.6");
    header.AddComment("a comment");
    header.AddProgramRecord(ProgramRecord{"tool"});
    header.AddReadGroup(ReadGroup{"rg1"});
    header.AddReferenceSequence(ReferenceSequence{"ref", 100});
    const std::string text = header.ToText();
    // @HD always first, then @SQ, @RG, @PG, @CO
    const std::string expected = "@HD\tVN:1.6\n"
                                 "@SQ\tSN:ref\tLN:100\n"
                                 "@RG\tID:rg1\n"
                                 "@PG\tID:tool\n"
                                 "@CO\ta comment\n";
    EXPECT_EQ(text, expected);
}

TEST(SamHeader, RoundTripSpecExample)
{
    const std::string original = "@HD\tVN:1.6\tSO:coordinate\n"
                                 "@SQ\tSN:ref\tLN:45\n";
    const SamHeader header = SamHeader::FromText(original);
    const std::string serialized = header.ToText();
    EXPECT_EQ(serialized, original);
}

TEST(SamHeader, RoundTripFullHeader)
{
    const std::string original = "@HD\tVN:1.6\tSO:coordinate\n"
                                 "@SQ\tSN:chr1\tLN:248956422\n"
                                 "@SQ\tSN:chrM\tLN:16569\tTP:circular\n"
                                 "@RG\tID:rg1\tSM:sample1\tPL:ILLUMINA\n"
                                 "@PG\tID:bwa\tPN:bwa\tVN:0.7.17\n"
                                 "@PG\tID:samtools\tPN:samtools\tVN:1.17\tPP:bwa\n"
                                 "@CO\tGenerated by test\n";
    const SamHeader header = SamHeader::FromText(original);
    const std::string serialized = header.ToText();
    EXPECT_EQ(serialized, original);
}

TEST(SamHeader, RoundTripNoHdLine)
{
    // No @HD → serialization should still produce a valid header
    // Convention: omit @HD entirely if version is empty
    const std::string original = "@SQ\tSN:ref\tLN:100\n";
    const SamHeader header = SamHeader::FromText(original);
    const std::string serialized = header.ToText();
    EXPECT_EQ(serialized, original);
}

TEST(SamHeader, SerializeSubSort)
{
    SamHeader header;
    header.SetVersion("1.6");
    header.SetSortOrder("coordinate");
    header.SetSubSort("coordinate:queryname");
    const std::string text = header.ToText();
    EXPECT_EQ(text, "@HD\tVN:1.6\tSO:coordinate\tSS:coordinate:queryname\n");
}

// --- Reference name/ID lookup tests ---

TEST(SamHeader, ReferenceNameToId)
{
    const std::string text = "@HD\tVN:1.6\n"
                             "@SQ\tSN:chr1\tLN:248956422\n"
                             "@SQ\tSN:chr2\tLN:242193529\n"
                             "@SQ\tSN:chrM\tLN:16569\n";
    const SamHeader header = SamHeader::FromText(text);

    EXPECT_EQ(header.ReferenceId("chr1"), 0);
    EXPECT_EQ(header.ReferenceId("chr2"), 1);
    EXPECT_EQ(header.ReferenceId("chrM"), 2);
    EXPECT_EQ(header.ReferenceId("chrX"), -1);  // not found
    EXPECT_EQ(header.ReferenceId("*"), -1);     // unmapped sentinel
}

TEST(SamHeader, ReferenceIdToName)
{
    const std::string text = "@HD\tVN:1.6\n"
                             "@SQ\tSN:chr1\tLN:248956422\n"
                             "@SQ\tSN:chr2\tLN:242193529\n";
    const SamHeader header = SamHeader::FromText(text);

    EXPECT_EQ(header.ReferenceName(0), "chr1");
    EXPECT_EQ(header.ReferenceName(1), "chr2");
    EXPECT_EQ(header.ReferenceName(-1), "*");  // unmapped sentinel
}

TEST(SamHeader, ReferenceIdToNameBoundsCheck)
{
    const std::string text = "@SQ\tSN:ref\tLN:100\n";
    const SamHeader header = SamHeader::FromText(text);

    // Out of bounds should throw
    EXPECT_THROW(header.ReferenceName(1), std::out_of_range);
    EXPECT_THROW(header.ReferenceName(-2), std::out_of_range);
}

TEST(SamHeader, ReferenceCount)
{
    const std::string text = "@SQ\tSN:chr1\tLN:100\n"
                             "@SQ\tSN:chr2\tLN:200\n";
    const SamHeader header = SamHeader::FromText(text);
    EXPECT_EQ(header.NumReferences(), 2);
}

TEST(SamHeader, EmptyHeaderLookup)
{
    const SamHeader header;
    EXPECT_EQ(header.ReferenceId("anything"), -1);
    EXPECT_EQ(header.NumReferences(), 0);
}

TEST(SamHeader, LookupAfterMutation)
{
    SamHeader header;
    header.AddReferenceSequence(ReferenceSequence{"chr1", 100});
    header.AddReferenceSequence(ReferenceSequence{"chr2", 200});

    EXPECT_EQ(header.ReferenceId("chr1"), 0);
    EXPECT_EQ(header.ReferenceId("chr2"), 1);

    // Add a third reference — lookup should still work
    header.AddReferenceSequence(ReferenceSequence{"chr3", 300});
    EXPECT_EQ(header.ReferenceId("chr3"), 2);
}

// --- BAM binary header parsing tests ---

namespace {

/// \brief Build a minimal BAM header block in memory.
/// Layout: magic(4) + l_text(4) + text(l_text) + n_ref(4) + [l_name(4) + name(l_name) + l_ref(4)]*
std::vector<std::byte> BuildBamHeader(std::string_view headerText,
                                      const std::vector<std::pair<std::string, std::int32_t>>& refs)
{
    std::vector<std::byte> data;

    // magic: BAM\1
    data.push_back(std::byte{'B'});
    data.push_back(std::byte{'A'});
    data.push_back(std::byte{'M'});
    data.push_back(std::byte{1});

    // l_text (little-endian uint32)
    const std::uint32_t lText = std::size(headerText);
    const std::byte* lTextBytes = reinterpret_cast<const std::byte*>(&lText);
    data.insert(std::ranges::end(data), lTextBytes, lTextBytes + 4);

    // text
    for (char c : headerText) {
        data.push_back(static_cast<std::byte>(c));
    }

    // n_ref (little-endian uint32)
    const std::uint32_t nRef = std::size(refs);
    const std::byte* nRefBytes = reinterpret_cast<const std::byte*>(&nRef);
    data.insert(std::ranges::end(data), nRefBytes, nRefBytes + 4);

    for (const auto& [name, length] : refs) {
        // l_name = reference name length + 1 (NUL terminator)
        const std::uint32_t lName = std::size(name) + 1;
        const std::byte* lNameBytes = reinterpret_cast<const std::byte*>(&lName);
        data.insert(std::ranges::end(data), lNameBytes, lNameBytes + 4);

        // name (NUL-terminated)
        for (char c : name) {
            data.push_back(static_cast<std::byte>(c));
        }
        data.push_back(std::byte{0});

        // l_ref (little-endian int32)
        const std::byte* lRefBytes = reinterpret_cast<const std::byte*>(&length);
        data.insert(std::ranges::end(data), lRefBytes, lRefBytes + 4);
    }

    return data;
}

}  // namespace

TEST(SamHeader, ParseBamHeaderBlock)
{
    const std::string headerText = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:ref\tLN:45\n";
    const std::vector<std::pair<std::string, std::int32_t>> refs = {{"ref", 45}};
    const std::vector<std::byte> data = BuildBamHeader(headerText, refs);
    const SamHeader header = SamHeader::FromBamHeaderBlock(std::span<const std::byte>{data});

    EXPECT_EQ(header.Version(), "1.6");
    EXPECT_EQ(header.SortOrder(), "coordinate");
    ASSERT_EQ(header.NumReferences(), 1);
    EXPECT_EQ(header.ReferenceSequences()[0].Name(), "ref");
    EXPECT_EQ(header.ReferenceSequences()[0].Length(), 45);
}

TEST(SamHeader, ParseBamHeaderBlockMultipleRefs)
{
    const std::string headerText = "@HD\tVN:1.6\n"
                                   "@SQ\tSN:chr1\tLN:1000\n"
                                   "@SQ\tSN:chr2\tLN:2000\n";
    const std::vector<std::pair<std::string, std::int32_t>> refs = {
        {"chr1", 1000},
        {"chr2", 2000},
    };
    const std::vector<std::byte> data = BuildBamHeader(headerText, refs);
    const SamHeader header = SamHeader::FromBamHeaderBlock(std::span<const std::byte>{data});

    ASSERT_EQ(header.NumReferences(), 2);
    EXPECT_EQ(header.ReferenceId("chr1"), 0);
    EXPECT_EQ(header.ReferenceId("chr2"), 1);
}

TEST(SamHeader, ParseBamHeaderBlockNoSqInText)
{
    // Some BAM files have reference info only in the binary dict, not @SQ lines
    const std::string headerText = "@HD\tVN:1.6\n";
    const std::vector<std::pair<std::string, std::int32_t>> refs = {
        {"chr1", 1000},
    };
    const std::vector<std::byte> data = BuildBamHeader(headerText, refs);
    const SamHeader header = SamHeader::FromBamHeaderBlock(std::span<const std::byte>{data});

    // Should use binary dict as authoritative source
    ASSERT_EQ(header.NumReferences(), 1);
    EXPECT_EQ(header.ReferenceSequences()[0].Name(), "chr1");
    EXPECT_EQ(header.ReferenceSequences()[0].Length(), 1000);
}

TEST(SamHeader, ParseBamHeaderBlockBadMagic)
{
    std::vector<std::byte> data = {std::byte{'X'}, std::byte{'X'}, std::byte{'X'}, std::byte{1}};
    // Pad enough to read l_text
    data.resize(12, std::byte{0});

    EXPECT_THROW(SamHeader::FromBamHeaderBlock(std::span<const std::byte>{data}),
                 std::runtime_error);
}

TEST(SamHeader, ParseBamHeaderBlockTruncated)
{
    // Only magic, too short for anything else
    const std::vector<std::byte> data = {std::byte{'B'}, std::byte{'A'}, std::byte{'M'},
                                         std::byte{1}};

    EXPECT_THROW(SamHeader::FromBamHeaderBlock(std::span<const std::byte>{data}),
                 std::runtime_error);
}

TEST(SamHeader, ParseBamHeaderBlockEmptyHeaderText)
{
    const std::string headerText;
    const std::vector<std::pair<std::string, std::int32_t>> refs = {{"ref", 100}};
    const std::vector<std::byte> data = BuildBamHeader(headerText, refs);
    const SamHeader header = SamHeader::FromBamHeaderBlock(std::span<const std::byte>{data});

    ASSERT_EQ(header.NumReferences(), 1);
    EXPECT_EQ(header.ReferenceSequences()[0].Name(), "ref");
}

TEST(SamHeader, ParseBamHeaderBlockNullPaddedText)
{
    // Some tools NUL-pad the header text
    std::string headerText = "@HD\tVN:1.6\n";
    headerText.push_back('\0');
    headerText.push_back('\0');

    const std::vector<std::pair<std::string, std::int32_t>> refs;
    const std::vector<std::byte> data = BuildBamHeader(headerText, refs);
    const SamHeader header = SamHeader::FromBamHeaderBlock(std::span<const std::byte>{data});

    EXPECT_EQ(header.Version(), "1.6");
}

// --- Integration tests with real BAM fixtures ---

TEST(SamHeader, ParseFromRealBamFile)
{
    // Find the test fixture
    const std::filesystem::path bamPath{tests::DataDir / "spec_example.bam"};
    ASSERT_TRUE(std::filesystem::exists(bamPath)) << bamPath;

    // Read enough decompressed data to cover the header
    BgzfReader reader{bamPath};
    std::array<std::byte, 65536> buffer;
    const std::optional<std::size_t> bytesRead = reader.ReadBlock(buffer);
    ASSERT_TRUE(bytesRead.has_value());
    ASSERT_GT(*bytesRead, 0U);

    // Parse the BAM header from decompressed bytes
    const SamHeader header =
        SamHeader::FromBamHeaderBlock(std::span<const std::byte>{buffer}.first(*bytesRead));

    // The spec example BAM has: @HD VN:1.6 SO:coordinate, @SQ SN:ref LN:45
    EXPECT_EQ(header.Version(), "1.6");
    EXPECT_EQ(header.SortOrder(), "coordinate");
    ASSERT_GE(header.NumReferences(), 1);
    EXPECT_EQ(header.ReferenceSequences()[0].Name(), "ref");
    EXPECT_EQ(header.ReferenceSequences()[0].Length(), 45);
    EXPECT_EQ(header.ReferenceId("ref"), 0);
}

TEST(SamHeader, ParseHeaderOnlyBam)
{
    const std::filesystem::path bamPath{tests::DataDir / "header_only.bam"};
    ASSERT_TRUE(std::filesystem::exists(bamPath)) << bamPath;

    BgzfReader reader{bamPath};
    std::array<std::byte, 65536> buffer;
    const std::optional<std::size_t> bytesRead = reader.ReadBlock(buffer);
    ASSERT_TRUE(bytesRead.has_value());
    ASSERT_GT(*bytesRead, 0U);
    const SamHeader header =
        SamHeader::FromBamHeaderBlock(std::span<const std::byte>{buffer}.first(*bytesRead));

    // Header-only BAM should parse without error
    EXPECT_FALSE(std::empty(header.Version()));
}

// --- BAM binary serialization tests ---

TEST(SamHeader, BamBinaryRoundTrip)
{
    SamHeader header;
    header.SetVersion("1.6");
    header.SetSortOrder("coordinate");
    header.AddReferenceSequence(ReferenceSequence{"chr1", 248956422});
    header.AddReferenceSequence(ReferenceSequence{"chrM", 16569});

    ReadGroup rg{"rg1"};
    rg.SetTag("SM", "sample1");
    header.AddReadGroup(std::move(rg));
    const std::vector<std::byte> binary = header.ToBamHeaderBlock();

    // Parse back
    const SamHeader parsed = SamHeader::FromBamHeaderBlock(std::span<const std::byte>{binary});

    EXPECT_EQ(parsed.Version(), "1.6");
    EXPECT_EQ(parsed.SortOrder(), "coordinate");
    ASSERT_EQ(parsed.NumReferences(), 2);
    EXPECT_EQ(parsed.ReferenceSequences()[0].Name(), "chr1");
    EXPECT_EQ(parsed.ReferenceSequences()[0].Length(), 248956422);
    EXPECT_EQ(parsed.ReferenceSequences()[1].Name(), "chrM");
    EXPECT_EQ(parsed.ReferenceSequences()[1].Length(), 16569);
    ASSERT_EQ(std::size(parsed.ReadGroups()), 1U);
    EXPECT_EQ(parsed.ReadGroups()[0].Id(), "rg1");
}

TEST(SamHeader, BamBinaryRoundTripEmptyHeader)
{
    const SamHeader header;
    const std::vector<std::byte> binary = header.ToBamHeaderBlock();
    const SamHeader parsed = SamHeader::FromBamHeaderBlock(std::span<const std::byte>{binary});

    EXPECT_EQ(parsed.NumReferences(), 0);
    EXPECT_TRUE(std::empty(parsed.Version()));
}

TEST(SamHeader, BamBinaryMagicBytes)
{
    SamHeader header;
    header.SetVersion("1.6");
    const std::vector<std::byte> binary = header.ToBamHeaderBlock();

    ASSERT_GE(std::size(binary), 4U);
    EXPECT_EQ(binary[0], std::byte{'B'});
    EXPECT_EQ(binary[1], std::byte{'A'});
    EXPECT_EQ(binary[2], std::byte{'M'});
    EXPECT_EQ(binary[3], std::byte{1});
}

TEST(SamHeader, BamBinaryFromRealFileRoundTrip)
{
    const std::filesystem::path bamPath{tests::DataDir / "spec_example.bam"};
    ASSERT_TRUE(std::filesystem::exists(bamPath));

    BgzfReader reader{bamPath};
    std::array<std::byte, 65536> buffer;
    const std::optional<std::size_t> bytesRead = reader.ReadBlock(buffer);
    ASSERT_TRUE(bytesRead.has_value());
    const SamHeader original =
        SamHeader::FromBamHeaderBlock(std::span<const std::byte>{buffer}.first(*bytesRead));

    // Serialize and re-parse
    const std::vector<std::byte> binary = original.ToBamHeaderBlock();
    const SamHeader reparsed = SamHeader::FromBamHeaderBlock(std::span<const std::byte>{binary});

    EXPECT_EQ(reparsed.Version(), original.Version());
    EXPECT_EQ(reparsed.SortOrder(), original.SortOrder());
    ASSERT_EQ(reparsed.NumReferences(), original.NumReferences());
    for (std::int32_t i = 0; i < original.NumReferences(); ++i) {
        const std::size_t idx = static_cast<std::size_t>(i);
        EXPECT_EQ(reparsed.ReferenceSequences()[idx].Name(),
                  original.ReferenceSequences()[idx].Name());
        EXPECT_EQ(reparsed.ReferenceSequences()[idx].Length(),
                  original.ReferenceSequences()[idx].Length());
    }
}

}  // namespace Samoa
}  // namespace PacBio
