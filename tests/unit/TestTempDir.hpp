#ifndef PBSAMOA_TESTS_UNIT_TESTTEMPDIR_HPP
#define PBSAMOA_TESTS_UNIT_TESTTEMPDIR_HPP

#include <gtest/gtest.h>

#include <atomic>
#include <cctype>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <format>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>

namespace PacBio {
namespace Samoa {
namespace tests {

inline char SanitizePathChar(unsigned char c)
{
    if (std::isalnum(c) != 0) {
        return static_cast<char>(c);
    }
    return '_';
}

inline std::string SanitizePathComponent(std::string_view input)
{
    std::string result;
    result.reserve(input.size());
    for (const unsigned char c : input) {
        result.push_back(SanitizePathChar(c));
    }
    if (result.empty()) {
        result = "unnamed";
    }
    return result;
}

inline std::filesystem::path MakeUniqueTempDirPath(std::string_view tag)
{
    const auto* testInfo = ::testing::UnitTest::GetInstance()->current_test_info();
    std::string suiteName{"suite"};
    std::string testName{"test"};
    if (testInfo) {
        suiteName = SanitizePathComponent(testInfo->test_suite_name());
        testName = SanitizePathComponent(testInfo->name());
    }

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

    TempDirGuard(TempDirGuard&& other) noexcept : path_{std::exchange(other.path_, {})} {}

    TempDirGuard& operator=(TempDirGuard&& other) noexcept
    {
        if (this != &other) {
            Cleanup();
            path_ = std::exchange(other.path_, {});
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

}  // namespace tests
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TESTS_UNIT_TESTTEMPDIR_HPP
