#include "TestTempDir.hpp"

#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/io/FastaReference.hpp>

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>

#include <cstdint>

namespace PacBio {
namespace Samoa {

namespace {

class FastaFileTest : public ::testing::Test
{
protected:
    tests::TempDirGuard tempDir_;

    void SetUp() override { tempDir_.Reset("fasta_reference"); }

    std::filesystem::path Write(std::string_view name, std::string_view text) const
    {
        const auto path{tempDir_.File(name)};
        std::ofstream out{path, std::ios::binary};
        out << text;
        return path;
    }

    static SamHeader Header(std::initializer_list<std::pair<std::string, std::int32_t>> sqs)
    {
        SamHeader header;
        for (const auto& [name, length] : sqs) {
            header.AddReferenceSequence(ReferenceSequence{name, length});
        }
        return header;
    }

    // Every window [start, start + length) of sequence refId equals the slice of expected.
    static void ExpectAllWindows(const CramReference& ref, std::int32_t refId,
                                 std::string_view expected)
    {
        ASSERT_EQ(ref.Length(refId), std::ssize(expected));
        std::string buffer;
        for (std::int64_t start = 0; start < std::ssize(expected); ++start) {
            for (std::int64_t length = 0; start + length <= std::ssize(expected); ++length) {
                EXPECT_EQ(ref.Bases(refId, start, length, buffer),
                          expected.substr(static_cast<std::size_t>(start),
                                          static_cast<std::size_t>(length)))
                    << "start " << start << " length " << length;
            }
        }
    }
};

}  // namespace

TEST_F(FastaFileTest, ScansRegularSequencesWithoutIndex)
{
    const auto fa{Write("a.fa", ">chr1 desc\nACGTA\ncgtac\nGT\n>chr2\nTTTT\nTT\n")};
    const FastaFile fasta{fa};
    EXPECT_EQ(fasta.NumSequences(), 2U);
    const auto ref{fasta.Reference(Header({{"chr2", 6}, {"chr1", 12}}))};
    ASSERT_EQ(ref->NumSequences(), 2);
    ExpectAllWindows(*ref, 0, "TTTTTT");
    ExpectAllWindows(*ref, 1, "ACGTAcgtacGT");
}

TEST_F(FastaFileTest, UsesTheFaiWhenPresent)
{
    // The .fai describes the file; with it, the offsets come from the index, not a scan.
    const auto fa{Write("b.fa", ">s1\nAAAC\nGG\n>s2\nCCCCA\n")};
    Write("b.fa.fai", "s1\t6\t4\t4\t5\ns2\t5\t16\t5\t6\n");
    const auto ref{FastaFile{fa}.Reference(Header({{"s1", 6}, {"s2", 5}}))};
    ExpectAllWindows(*ref, 0, "AAACGG");
    ExpectAllWindows(*ref, 1, "CCCCA");
}

TEST_F(FastaFileTest, RejectsAFaiThatDoesNotFitTheFile)
{
    const auto fa{Write("c.fa", ">s1\nACGT\n")};
    Write("c.fa.fai", "s1\t400\t4\t4\t5\n");
    EXPECT_THROW(FastaFile{fa}, std::runtime_error);
}

TEST_F(FastaFileTest, HandlesCrlfAndAMissingFinalNewline)
{
    const auto fa{Write("d.fa", ">s1\r\nACG\r\nTTA\r\nC\r\n>s2\nGGG\nCCC")};
    const auto ref{FastaFile{fa}.Reference(Header({{"s1", 7}, {"s2", 6}}))};
    ExpectAllWindows(*ref, 0, "ACGTTAC");
    ExpectAllWindows(*ref, 1, "GGGCCC");
}

TEST_F(FastaFileTest, ReadsUnevenLinesIntoMemory)
{
    const auto fa{Write("e.fa", ">s1\nAC\nGTA\nC\n\nGG\n>s2\nAC GT\n")};
    const auto ref{FastaFile{fa}.Reference(Header({{"s1", 8}, {"s2", 4}}))};
    ExpectAllWindows(*ref, 0, "ACGTACGG");
    ExpectAllWindows(*ref, 1, "ACGT");
}

TEST_F(FastaFileTest, MissingSequenceIsEmptyOrThrows)
{
    const auto fa{Write("f.fa", ">s1\nACGT\n")};
    const FastaFile fasta{fa};
    const auto header{Header({{"s1", 4}, {"other", 10}})};
    EXPECT_EQ(fasta.Reference(header)->Length(1), 0);
    EXPECT_THROW((void)fasta.Reference(header, MissingFastaSequence::THROW), std::runtime_error);
}

TEST_F(FastaFileTest, GzipAndEmptyFilesHaveNoSequences)
{
    EXPECT_EQ(FastaFile{Write("g.fa.gz", std::string{"\x1f\x8b\x08\x00rest", 8})}.NumSequences(),
              0U);
    EXPECT_EQ(FastaFile{Write("h.fa", "")}.NumSequences(), 0U);
    EXPECT_THROW(FastaFile{tempDir_.File("missing.fa")}, std::runtime_error);
}

TEST_F(FastaFileTest, ReferenceOutlivesTheFastaFile)
{
    std::shared_ptr<const CramReference> ref;
    {
        const FastaFile fasta{Write("i.fa", ">s1\nACGT\nAC\n")};
        ref = fasta.Reference(Header({{"s1", 6}}));
    }
    std::string buffer;
    EXPECT_EQ(ref->Bases(0, 2, 4, buffer), "GTAC");
}

}  // namespace Samoa
}  // namespace PacBio
