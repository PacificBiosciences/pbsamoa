#ifndef PBSAMOA_IO_SAMWRITER_HPP
#define PBSAMOA_IO_SAMWRITER_HPP

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <filesystem>
#include <memory>

namespace PacBio {
namespace Samoa {

struct SamWriterConfig
{
    /// \brief Write through a temporary file and atomically rename on Close().
    bool UseTempFile{false};
};

/// \brief Writes SAM text files.
///
/// Accepts both owned BamRecord and RawRecord. Fields are formatted
/// per SAM spec (1-based POS, "=" for RNEXT when same as RNAME, etc.).
class SamWriter
{
public:
    /// \param[in] path output SAM file path
    /// \param[in] header SAM header to write
    SamWriter(const std::filesystem::path& path, const SamHeader& header,
              const SamWriterConfig& config = SamWriterConfig{});
    ~SamWriter();

    SamWriter(const SamWriter&) = delete;
    SamWriter& operator=(const SamWriter&) = delete;
    SamWriter(SamWriter&&) noexcept;
    SamWriter& operator=(SamWriter&&) noexcept;

    /// \brief Write an owned record.
    void Write(const BamRecord& record);

    /// \brief Write an owning byte view.
    void Write(const RawRecord& view);

    /// \brief Write a batch of records.
    void WriteBatch(const RawRecordBatch& batch);

    /// \brief Flush and close. Called automatically by destructor.
    void Close();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_IO_SAMWRITER_HPP
