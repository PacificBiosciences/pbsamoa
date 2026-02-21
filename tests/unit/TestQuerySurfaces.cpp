#include "TestData.hpp"

#include <pbsamoa/core/GenomicInterval.hpp>
#include <pbsamoa/index/BaiIndex.hpp>
#include <pbsamoa/io/BamFile.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamRecordReader.hpp>

#include <gtest/gtest.h>

#include <filesystem>
#include <format>
#include <system_error>
#include <utility>
#include <vector>

#include <unistd.h>

namespace PacBio {
namespace Samoa {
namespace {

class TempFileGuard
{
public:
    explicit TempFileGuard(std::filesystem::path path) : path_{std::move(path)} {}

    TempFileGuard(const TempFileGuard&) = delete;
    TempFileGuard& operator=(const TempFileGuard&) = delete;
    TempFileGuard(TempFileGuard&&) = delete;
    TempFileGuard& operator=(TempFileGuard&&) = delete;

    const std::filesystem::path& Path() const { return path_; }

    ~TempFileGuard()
    {
        std::error_code ec;
        std::filesystem::remove(path_, ec);
    }

private:
    std::filesystem::path path_;
};

std::size_t CountRawQuery(const std::filesystem::path& bamPath, const GenomicInterval& interval)
{
    const BamFile file{bamPath};
    const BaiIndex index{BaiIndex::FromFile(file.StandardIndexFilename())};
    BamRawReader reader{bamPath};

    std::size_t count{0};
    for (const auto& record : reader.Query(index, file.ReferenceId(interval.Name()),
                                           interval.Start(), interval.Stop())) {
        std::ignore = record;
        ++count;
    }
    return count;
}

template <typename Query>
std::size_t CountRecords(Query& query)
{
    std::size_t count{0};
    for (const auto& record : query) {
        std::ignore = record;
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
