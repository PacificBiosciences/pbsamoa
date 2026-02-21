#ifndef PBSAMOA_IO_BAMZMWREADER_HPP
#define PBSAMOA_IO_BAMZMWREADER_HPP

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/index/ZmwIndex.hpp>
#include <pbsamoa/io/BamCollection.hpp>
#include <pbsamoa/io/BamFile.hpp>
#include <pbsamoa/io/BamRecordReader.hpp>

#include <filesystem>
#include <memory>
#include <vector>

namespace PacBio {
namespace Samoa {

/// \brief A single ZMW's worth of owned records.
struct ZmwGroup
{
    ZmwIdentity zmw;
    std::vector<BamRecord> records;
};

/// \brief Reads BAM records grouped by ZMW identity.
///
/// Wraps a BamRecordReader and yields groups of consecutive records that share
/// the same `(rgId, zmw)`, where `rgId` is parsed from the `RG` tag's base
/// PacBio read-group ID when present and falls back to `0` otherwise.
/// Assumes the source yields records already grouped by this identity.
class BamZmwReader
{
public:
    /// \param[in] reader pre-decode record reader (takes ownership)
    explicit BamZmwReader(BamRecordReader reader);
    explicit BamZmwReader(BamCollection collection, BamRecordReaderConfig config = {});
    explicit BamZmwReader(std::vector<std::filesystem::path> paths,
                          BamRecordReaderConfig config = {});
    explicit BamZmwReader(std::vector<BamFile> files, BamRecordReaderConfig config = {});
    ~BamZmwReader();

    BamZmwReader(const BamZmwReader&) = delete;
    BamZmwReader& operator=(const BamZmwReader&) = delete;
    BamZmwReader(BamZmwReader&&) noexcept;
    BamZmwReader& operator=(BamZmwReader&&) noexcept;

    /// \brief Get the next group of records sharing the same ZMW.
    /// \param[out] records filled with owned BamRecords for one ZMW
    /// \returns true if a group was read, false at EOF
    bool GetNext(std::vector<BamRecord>& records);

    /// \brief The ZMW identity of the last returned group.
    ZmwIdentity CurrentZmw() const;

    /// \brief Access the header from the underlying source.
    const SamHeader& Header() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_IO_BAMZMWREADER_HPP
