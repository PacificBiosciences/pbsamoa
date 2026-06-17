#include "TestData.hpp"
#include "TestTempDir.hpp"

#include <pbsamoa/core/GenomicInterval.hpp>
#include <pbsamoa/index/BaiIndex.hpp>
#include <pbsamoa/io/BamFile.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamRecordReader.hpp>

#include <gtest/gtest.h>

#include <filesystem>
#include <format>
#include <vector>

#include <unistd.h>

namespace PacBio {
namespace Samoa {
namespace {

using tests::TempFileGuard;

std::size_t CountRawQuery(const std::filesystem::path& bamPath, const GenomicInterval& interval)
{
    const BamFile file{bamPath};
    const BaiIndex index{BaiIndex::FromFile(file.StandardIndexFilename())};
    BamRawReader reader{bamPath};

    std::size_t count{0};
    for (const auto& record : reader.Query(index, file.Header().ReferenceId(interval.Name()),
                                           interval.Start(), interval.Stop())) {
        (void)record;
        ++count;
    }
    return count;
}

template <typename Query>
std::size_t CountRecords(Query& query)
{
    std::size_t count{0};
    for (const auto& record : query) {
        (void)record;
        ++count;
    }
    return count;
}

}  // namespace

TEST(BamRecordReaderQuery, CanReuseIntervalsAndHandleUnknownReferences)
{
    const BamFile file{tests::DataDir / "spec_example.bam"};
    BamRecordReader reader{file.Filename(), BamRecordReaderConfig{.DecodeWorkers = 0}};

    const GenomicInterval fullInterval{"ref", 0, 45};
    const auto fullExpected = CountRawQuery(file.Filename(), fullInterval);
    auto fullRange = reader.Query(fullInterval);
    EXPECT_EQ(CountRecords(fullRange), fullExpected);

    const GenomicInterval tailInterval{"ref", 35, 45};
    const auto tailExpected = CountRawQuery(file.Filename(), tailInterval);
    auto tailRange = reader.Query("ref", 35, 45);
    EXPECT_EQ(CountRecords(tailRange), tailExpected);

    EXPECT_THROW(reader.Query(GenomicInterval{"missing", 0, 45}), std::runtime_error);
}

TEST(BamRecordReaderQuery, ThrowsOnMissingBai)
{
    const TempFileGuard tmpBam{std::filesystem::temp_directory_path() /
                               std::format("pbsamoa_missing_bai_{}.bam", getpid())};
    const TempFileGuard tmpBai{std::filesystem::path{tmpBam.Path().string() + ".bai"}};
    std::filesystem::copy_file(tests::DataDir / "spec_example.bam", tmpBam.Path(),
                               std::filesystem::copy_options::overwrite_existing);
    std::filesystem::remove(tmpBai.Path());

    BamRecordReader reader{tmpBam.Path(), BamRecordReaderConfig{.DecodeWorkers = 0}};
    EXPECT_THROW(reader.Query("ref", 0, 45), std::runtime_error);
}

}  // namespace Samoa
}  // namespace PacBio
