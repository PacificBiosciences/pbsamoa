#ifndef PBSAMOA_TOOLS_STATUSSTREAM_HPP
#define PBSAMOA_TOOLS_STATUSSTREAM_HPP

#include <pbcopper/json/JSON.h>

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <filesystem>
#include <fstream>
#include <mutex>
#include <optional>
#include <stop_token>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {
namespace Tools {

/// \brief The live status stream of a command (--status-json): NDJSON, one JSON object per line.
///        A copy of the StatusStream of pbweave, with warnings on stderr in place of pbcopper
///        logging.
///
/// The first line is `run_info`. Each phase writes a `start` and an `end` line. While a phase with
/// a counter is active, a reporter thread writes a `progress` line each second if the counter
/// changed. The `final` line closes the stream. The contract is `tool-status-stream.md` of
/// pbhuman. Writing is best-effort: an I/O error logs one warning and disables the stream. The
/// run never fails because of the stream.
class StatusStream
{
public:
    /// \brief A phase of the command. The constructor writes the `start` line.
    ///
    /// The worker threads add to a counter of the caller with std::memory_order::relaxed. The
    /// reporter reads it. End() writes the last `progress` line and the `end` line. A phase that
    /// does not call End() (an exception) writes no more lines.
    class Phase
    {
    public:
        /// \param[in] stream  the stream; nullptr or a disabled stream makes all calls no-ops
        /// \param[in] name    the name of the phase, as in run_info.phases
        /// \param[in] done    the counter of the phase, or nullptr for a phase without progress;
        ///                    it must stay valid until End() or the destructor
        /// \param[in] unit    what the counter counts (`reads`, `records`, ...)
        /// \param[in] total   the value of the counter at the end, if it is known
        Phase(StatusStream* stream, std::string name,
              const std::atomic<std::int64_t>* done = nullptr, std::string unit = {},
              std::optional<std::int64_t> total = std::nullopt);

        // stops the progress lines of the phase if End() did not
        ~Phase();

        Phase(const Phase&) = delete;
        Phase& operator=(const Phase&) = delete;
        Phase(Phase&&) = delete;
        Phase& operator=(Phase&&) = delete;

        /// \brief Stops the progress lines of the phase, writes the last `progress` line (done
        ///        equal to total) and the `end` line.
        /// \param[in] metrics  the numbers that the log gives at the end of the phase (optional)
        void End(JSON::Json metrics = JSON::Json::object()) noexcept;

    private:
        friend class StatusStream;

        StatusStream* stream_;
        std::string name_;
        const std::atomic<std::int64_t>* done_;
        std::string unit_;
        std::optional<std::int64_t> total_;
        bool active_{false};
    };

    /// \brief A disabled stream: each call is a no-op until Open().
    StatusStream() noexcept;

    // stops the reporter thread
    ~StatusStream();

    StatusStream(const StatusStream&) = delete;
    StatusStream& operator=(const StatusStream&) = delete;
    StatusStream(StatusStream&&) = delete;
    StatusStream& operator=(StatusStream&&) = delete;

    /// \brief Creates (truncates) the file, writes the `run_info` line and starts the reporter.
    /// \param[in] path     the output file; the parent directory is created. Empty keeps the
    ///                     stream disabled.
    /// \param[in] runInfo  the fields of the first line (`type` is added here). Its
    ///                     `num_threads` also goes into the resources of the `final` line.
    void Open(const std::filesystem::path& path, JSON::Json runInfo) noexcept;

    /// \returns true while the stream writes lines
    [[nodiscard]] bool Enabled() const noexcept;

    /// \brief Appends one line with `type`, `t` (seconds since Open) and `wall_ms`. Thread-safe.
    ///        No-op after a write error.
    void Write(std::string_view type, JSON::Json line) noexcept;

    /// \brief Stops the reporter and writes the `final` line with the resources of the process.
    ///        After it, the stream writes no more lines.
    /// \param[in] phase    `FINISHED` or `FAILED`
    /// \param[in] cause    the message of the exception of a failure; empty for none
    /// \param[in] metrics  the main numbers of the run summary
    /// \param[in] outputs  the files of the command: objects with `path` and `kind`
    void Final(std::string_view phase, std::string_view cause = {},
               JSON::Json metrics = JSON::Json::object(),
               JSON::Json outputs = JSON::Json::array()) noexcept;

private:
    // writes one line; the caller holds mutex_
    void WriteLocked(std::string_view type, JSON::Json& line);

    // writes a `progress` line of the phase; the caller holds mutex_
    void WriteProgressLocked(const Phase& phase, std::int64_t done,
                             std::optional<std::int64_t> total);

    // disables the stream and logs one warning; the caller holds mutex_
    void DisableLocked() noexcept;

    // the loop of the reporter thread
    void Report(std::stop_token stop);

    std::filesystem::path path_;
    std::ofstream out_;
    std::chrono::steady_clock::time_point start_;
    std::int64_t numThreads_{0};
    bool disabled_{true};
    bool finalWritten_{false};
    // the phase with a counter that the reporter reports, and its last reported value
    const Phase* active_{nullptr};
    std::int64_t reported_{-1};
    // one mutex for the file and the active phase: a progress line of the reporter never comes
    // after the end line of its phase
    mutable std::mutex mutex_;
    std::condition_variable_any wake_;
    // declared last: it stops before the members that it reads are destroyed
    std::jthread reporter_;
};

/// \brief The fields of the `run_info` line, with `schema_version`, `pid`, `hostname` and
///        `started_wall_ms` of this process.
/// \param[in] phases  the phases of this run, in their order
[[nodiscard]] JSON::Json MakeRunInfo(std::string_view tool, std::string_view command,
                                     std::string_view version, std::string_view commandLine,
                                     std::int64_t numThreads,
                                     const std::vector<std::string>& phases);

}  // namespace Tools
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_STATUSSTREAM_HPP
