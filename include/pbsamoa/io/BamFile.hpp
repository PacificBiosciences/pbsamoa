#ifndef PBSAMOA_IO_BAMFILE_HPP
#define PBSAMOA_IO_BAMFILE_HPP

#include <pbsamoa/core/SamHeader.hpp>

#include <filesystem>
#include <string>
#include <string_view>

#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief Convenience wrapper around a BAM path, header, and `.bai` sidecar.
class BamFile
{
public:
    explicit BamFile(std::filesystem::path fileName);

    const std::filesystem::path& Filename() const;

    void CreateStandardIndex() const;
    void EnsureStandardIndexExists() const;

    bool HasEOF() const;
    bool StandardIndexExists() const;
    std::filesystem::path StandardIndexFilename() const;
    bool StandardIndexIsNewer() const;

    const SamHeader& Header() const;
    std::string ReferenceName(std::int32_t id) const;
    std::uint32_t ReferenceLength(std::string_view name) const;
    std::uint32_t ReferenceLength(std::int32_t id) const;

private:
    std::filesystem::path fileName_;
    SamHeader header_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_IO_BAMFILE_HPP
