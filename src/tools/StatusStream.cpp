#include "StatusStream.hpp"

#include <array>
#include <exception>
#include <format>
#include <print>
#include <stdexcept>
#include <utility>

#include <cmath>

#include <sys/resource.h>
#include <unistd.h>

namespace PacBio {
namespace Samoa {
namespace Tools {
namespace {

// the reporter writes at most one progress line in this interval
constexpr std::chrono::seconds REPORT_INTERVAL{1};

std::int64_t WallMilliseconds()
{
    return std::chrono::duration_cast<std::chrono::milliseconds>(
               std::chrono::system_clock::now().time_since_epoch())
        .count();
}

// seconds, rounded to milliseconds: short lines, and no false precision
double RoundSeconds(const double seconds) { return std::round(seconds * 1000.0) / 1000.0; }

}  // namespace

StatusStream::Phase::Phase(StatusStream* const stream, std::string name,
                           const std::atomic<std::int64_t>* const done, std::string unit,
                           const std::optional<std::int64_t> total)
    : stream_{stream}, name_{std::move(name)}, done_{done}, unit_{std::move(unit)}, total_{total}
{
    if (stream_ == nullptr) {
        return;
    }
    const std::scoped_lock lock{stream_->mutex_};
    if (stream_->disabled_) {
        return;
    }
    try {
        JSON::Json line{{"name", name_}, {"state", "start"}};
        stream_->WriteLocked("phase", line);
        active_ = true;
        if (done_ != nullptr) {
            stream_->active_ = this;
            stream_->reported_ = -1;
        }
    } catch (...) {
        stream_->DisableLocked();
    }
}

StatusStream::Phase::~Phase()
{
    if (!active_) {
        return;
    }
    const std::scoped_lock lock{stream_->mutex_};
    if (stream_->active_ == this) {
        stream_->active_ = nullptr;
    }
}

void StatusStream::Phase::End(JSON::Json metrics) noexcept
{
    if (!active_) {
        return;
    }
    active_ = false;
    const std::scoped_lock lock{stream_->mutex_};
    try {
        // the order of the contract: stop the reporter for the phase, write the last progress
        // line, then the end line
        if (stream_->active_ == this) {
            stream_->active_ = nullptr;
        }
        if (done_ != nullptr) {
            // the end gives the exact count: done equals total, also if total was not known
            const std::int64_t done{done_->load(std::memory_order::relaxed)};
            stream_->WriteProgressLocked(*this, done, done);
        }
        JSON::Json line{{"name", name_}, {"state", "end"}};
        if (!metrics.empty()) {
            line["metrics"] = std::move(metrics);
        }
        stream_->WriteLocked("phase", line);
    } catch (...) {
        stream_->DisableLocked();
    }
}

StatusStream::StatusStream() noexcept = default;

StatusStream::~StatusStream() = default;

void StatusStream::Open(const std::filesystem::path& path, JSON::Json runInfo) noexcept
{
    if (std::empty(path)) {
        return;
    }
    try {
        const std::filesystem::path parent{path.parent_path()};
        if (!std::empty(parent)) {
            std::filesystem::create_directories(parent);
        }
        out_.open(path, std::ios::out | std::ios::trunc);
        if (!out_.is_open()) {
            throw std::runtime_error{std::format("cannot open {}", path.string())};
        }
        path_ = path;
        start_ = std::chrono::steady_clock::now();
        if (runInfo.contains("num_threads") && runInfo["num_threads"].is_number_integer()) {
            numThreads_ = runInfo["num_threads"].get<std::int64_t>();
        }
        {
            const std::scoped_lock lock{mutex_};
            disabled_ = false;
            WriteLocked("run_info", runInfo);
        }
        reporter_ = std::jthread{[this](const std::stop_token stop) { Report(stop); }};
    } catch (const std::exception& e) {
        const std::scoped_lock lock{mutex_};
        disabled_ = true;
        std::println(stderr, "pbsamoa: WARNING: the status stream is disabled: {}", e.what());
    }
}

bool StatusStream::Enabled() const noexcept
{
    const std::scoped_lock lock{mutex_};
    return !disabled_;
}

void StatusStream::Write(const std::string_view type, JSON::Json line) noexcept
{
    const std::scoped_lock lock{mutex_};
    try {
        WriteLocked(type, line);
    } catch (...) {
        DisableLocked();
    }
}

void StatusStream::Final(const std::string_view phase, const std::string_view cause,
                         JSON::Json metrics, JSON::Json outputs) noexcept
{
    try {
        // the reporter stops before the final line
        reporter_.request_stop();
        if (reporter_.joinable()) {
            reporter_.join();
        }
    } catch (...) {
        // a reporter that does not join cannot write after the final line: it needs mutex_
    }
    const std::scoped_lock lock{mutex_};
    if (disabled_ || finalWritten_) {
        return;
    }
    finalWritten_ = true;
    try {
        rusage usage{};
        getrusage(RUSAGE_SELF, &usage);
        const double cpuSeconds{usage.ru_utime.tv_sec + (usage.ru_utime.tv_usec / 1e6) +
                                usage.ru_stime.tv_sec + (usage.ru_stime.tv_usec / 1e6)};
#if defined(__APPLE__)
        const std::int64_t peakRss{usage.ru_maxrss};  // bytes
#else
        const std::int64_t peakRss{usage.ru_maxrss * std::int64_t{1024}};  // KiB
#endif
        const std::chrono::duration<double> wall{std::chrono::steady_clock::now() - start_};
        JSON::Json line{
            {"phase", phase},
            {"resources",
             {{"wall_s", RoundSeconds(wall.count())},
              {"cpu_s", RoundSeconds(cpuSeconds)},
              {"peak_rss_bytes", peakRss},
              {"num_threads", numThreads_}}},
            {"metrics", std::move(metrics)},
            {"outputs", std::move(outputs)},
        };
        if (!cause.empty()) {
            line["cause"] = cause;
        }
        WriteLocked("final", line);
        // the final line is the last line
        disabled_ = true;
    } catch (...) {
        DisableLocked();
    }
}

void StatusStream::WriteLocked(const std::string_view type, JSON::Json& line)
{
    if (disabled_) {
        return;
    }
    line["type"] = type;
    if (type != "run_info") {
        const std::chrono::duration<double> t{std::chrono::steady_clock::now() - start_};
        line["t"] = RoundSeconds(t.count());
        line["wall_ms"] = WallMilliseconds();
    }
    // replace bytes that are not UTF-8 (a file name, an exception message): dump() throws else
    out_ << line.dump(-1, ' ', false, JSON::Json::error_handler_t::replace) << '\n';
    out_.flush();
    if (!out_) {
        DisableLocked();
    }
}

void StatusStream::WriteProgressLocked(const Phase& phase, const std::int64_t done,
                                       const std::optional<std::int64_t> total)
{
    JSON::Json line{
        {"phase", phase.name_},
        {"done", done},
        {"total", total ? JSON::Json(*total) : JSON::Json(nullptr)},
        {"unit", phase.unit_},
    };
    WriteLocked("progress", line);
}

void StatusStream::DisableLocked() noexcept
{
    if (disabled_) {
        return;
    }
    disabled_ = true;
    try {
        std::println(stderr,
                     "pbsamoa: WARNING: the status stream is disabled after a write error on {}",
                     path_.string());
    } catch (...) {
        // best-effort by contract: the stream never stops the run
    }
}

void StatusStream::Report(const std::stop_token stop)
{
    std::unique_lock lock{mutex_};
    while (true) {
        // wakes at once on request_stop()
        wake_.wait_for(lock, stop, REPORT_INTERVAL, [] { return false; });
        if (stop.stop_requested()) {
            return;
        }
        if ((active_ == nullptr) || disabled_) {
            continue;
        }
        // the check of the active phase, the read and the write are under one lock: a progress
        // line never comes after the end line of its phase
        const std::int64_t done{active_->done_->load(std::memory_order::relaxed)};
        if (done == reported_) {
            continue;
        }
        reported_ = done;
        try {
            WriteProgressLocked(*active_, done, active_->total_);
        } catch (...) {
            DisableLocked();
        }
    }
}

JSON::Json MakeRunInfo(const std::string_view tool, const std::string_view command,
                       const std::string_view version, const std::string_view commandLine,
                       const std::int64_t numThreads, const std::vector<std::string>& phases)
{
    std::array<char, 256> host{};
    std::string hostname;
    if (gethostname(host.data(), host.size() - 1) == 0) {
        hostname = host.data();
    }
    return JSON::Json{
        {"schema_version", 1},
        {"tool", tool},
        {"command", command},
        {"version", version},
        {"commandline", commandLine},
        {"pid", getpid()},
        {"hostname", hostname},
        {"num_threads", numThreads},
        {"started_wall_ms", WallMilliseconds()},
        {"phases", phases},
    };
}

}  // namespace Tools
}  // namespace Samoa
}  // namespace PacBio
