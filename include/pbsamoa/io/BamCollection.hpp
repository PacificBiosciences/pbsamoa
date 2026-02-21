#ifndef PBSAMOA_IO_BAMCOLLECTION_HPP
#define PBSAMOA_IO_BAMCOLLECTION_HPP

#include <pbsamoa/io/BamFile.hpp>

#include <filesystem>
#include <span>
#include <vector>

#include <cstddef>

namespace PacBio {
namespace Samoa {

/// \brief Ordered set of BAM inputs with a merged compatible header.
class BamCollection
{
public:
    explicit BamCollection(std::filesystem::path bam);
    explicit BamCollection(BamFile bamFile);
    explicit BamCollection(std::vector<std::filesystem::path> bams);
    explicit BamCollection(std::vector<BamFile> bamFiles);

    const SamHeader& Header() const;
    std::span<const BamFile> Files() const;
    std::size_t Size() const;

private:
    std::vector<BamFile> files_;
    SamHeader header_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_IO_BAMCOLLECTION_HPP
