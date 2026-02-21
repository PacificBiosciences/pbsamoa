#ifndef PBSAMOA_IO_ZMIWRITER_HPP
#define PBSAMOA_IO_ZMIWRITER_HPP

#include <filesystem>
#include <memory>

#include <cstdint>

namespace PacBio {
namespace Samoa {

struct ZmiWriterConfig
{
    /// \brief Write through a temporary file and atomically rename on Close().
    bool UseTempFile{false};
};

/// \brief Streaming writer for ZMW index (.zmi) files.
///
/// Appends (rgId, zmw, virtualOffset) entries in BAM write order.
/// BGZF-compressed. Header is finalized with numRecords on close.
class ZmiWriter
{
public:
    explicit ZmiWriter(const std::filesystem::path& path,
                       const ZmiWriterConfig& config = ZmiWriterConfig{});
    ~ZmiWriter();

    ZmiWriter(const ZmiWriter&) = delete;
    ZmiWriter& operator=(const ZmiWriter&) = delete;
    ZmiWriter(ZmiWriter&&) noexcept;
    ZmiWriter& operator=(ZmiWriter&&) noexcept;

    /// \brief Append one entry to the index.
    /// \param[in] rgId numeric read-group ID
    /// \param[in] zmw ZMW hole number
    /// \param[in] virtualOffset BGZF virtual offset of the BAM record
    void AddRecord(std::int32_t rgId, std::int32_t zmw, std::int64_t virtualOffset);

    /// \brief Flush and close. Called automatically by destructor.
    void Close();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_IO_ZMIWRITER_HPP
