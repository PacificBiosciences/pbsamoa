#ifndef PBSAMOA_TESTS_UNIT_TESTTEMPDIR_HPP
#define PBSAMOA_TESTS_UNIT_TESTTEMPDIR_HPP

#include <gtest/gtest.h>

#include <atomic>
#include <chrono>
#include <filesystem>
#include <format>
#include <string>
#include <string_view>
#include <system_error>

#include <cctype>
#include <cstdint>

namespace PacBio {
namespace Samoa {
namespace tests {

inline std::string SanitizePathComponent(std::string_view input)
{
    std::string result;
    result.reserve(input.size());
    for (const unsigned char c : input) {
        result.push_back(std::isalnum(c) ? static_cast<char>(c) : '_');
    }
    if (result.empty()) {
        result = "unnamed";
    }
    return result;
}

inline std::filesystem::path MakeUniqueTempDirPath(std::string_view tag)
{
    const auto* testInfo = ::testing::UnitTest::GetInstance()->current_test_info();
    const std::string suiteName{
        testInfo != nullptr ? SanitizePathComponent(testInfo->test_suite_name()) : "suite"};
    const std::string testName{testInfo != nullptr ? SanitizePathComponent(testInfo->name())
                                                   : "test"};

    static std::atomic<std::uint64_t> counter{0};
    const std::uint64_t nonce{counter.fetch_add(1, std::memory_order_relaxed)};
    const std::uint64_t nowNs{
        static_cast<std::uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(
                                       std::chrono::steady_clock::now().time_since_epoch())
                                       .count())};

    return std::filesystem::temp_directory_path() / "pbsamoa_tests" /
           std::format("{}_{}_{}_{}_{}", SanitizePathComponent(tag), suiteName, testName, nowNs,
                       nonce);
}

class TempDirGuard
{
public:
    TempDirGuard() = default;

    explicit TempDirGuard(std::string_view tag) { Reset(tag); }

    TempDirGuard(const TempDirGuard&) = delete;
    TempDirGuard& operator=(const TempDirGuard&) = delete;

    TempDirGuard(TempDirGuard&& other) noexcept : path_{std::move(other.path_)}
    {
        other.path_.clear();
    }

    TempDirGuard& operator=(TempDirGuard&& other) noexcept
    {
        if (this != &other) {
            Cleanup();
            path_ = std::move(other.path_);
            other.path_.clear();
        }
        return *this;
    }

    ~TempDirGuard() { Cleanup(); }

    void Reset(std::string_view tag)
    {
        Cleanup();
        path_ = MakeUniqueTempDirPath(tag);
        std::filesystem::create_directories(path_);
    }

    const std::filesystem::path& Path() const { return path_; }

    std::filesystem::path File(std::string_view filename) const
    {
        return path_ / std::filesystem::path{filename};
    }

private:
    void Cleanup() noexcept
    {
        if (path_.empty()) {
            return;
        }
        std::error_code ec;
        std::filesystem::remove_all(path_, ec);
        path_.clear();
    }

    std::filesystem::path path_;
};

}  // namespace tests
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TESTS_UNIT_TESTTEMPDIR_HPP
