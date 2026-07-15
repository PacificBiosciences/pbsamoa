#ifndef PBSAMOA_IO_ZMIBAMWRITER_HPP
#define PBSAMOA_IO_ZMIBAMWRITER_HPP

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/io/BamWriter.hpp>

#include <filesystem>
#include <memory>

namespace PacBio {
namespace Samoa {

/// \brief Writes BAM + ZMI index files together.
///
/// Convenience wrapper that owns a BamWriter and ZmiWriter, wiring them
/// so every record written to the BAM also gets a ZMI index entry.
/// The .zmi file is created at `bamPath + ".zmi"`.
class ZmiBamWriter
{
public:
    /// \param[in] bamPath output BAM file path (.zmi created alongside)
    /// \param[in] header SAM header
    ZmiBamWriter(const std::filesystem::path& bamPath, const SamHeader& header,
                 const BamWriterConfig& config = BamWriterConfig{});
    ~ZmiBamWriter();

    ZmiBamWriter(const ZmiBamWriter&) = delete;
    ZmiBamWriter& operator=(const ZmiBamWriter&) = delete;
    ZmiBamWriter(ZmiBamWriter&&) noexcept;
    ZmiBamWriter& operator=(ZmiBamWriter&&) noexcept;

    /// \brief Write an owned record.
    void Write(const BamRecord& record);

    /// \brief Write an owning byte view.
    void Write(const RawRecord& view);

    /// \brief Write a batch of records.
    void WriteBatch(const RawRecordBatch& batch);

    /// \brief Flush and close both BAM and ZMI. Called automatically by
    /// destructor.
    void Close();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_IO_ZMIBAMWRITER_HPP
