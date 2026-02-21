#ifndef PBSAMOA_CORE_CIGARCLIPPING_HPP
#define PBSAMOA_CORE_CIGARCLIPPING_HPP

#include <pbsamoa/core/CigarOp.hpp>

#include <span>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief Result of a CIGAR clipping operation.
struct ClipResult
{
    std::vector<CigarOp> cigar;
    std::size_t clipOffset;  ///< query-base offset into original data
    std::size_t clipLength;  ///< query bases retained
    std::int32_t newPos;     ///< updated reference position
};

/// \brief Clip CIGAR to query (polymerase/ZMW) coordinates.
ClipResult ClipCigarToQuery(std::span<const CigarOp> cigar, std::int32_t queryStart,
                            std::int32_t queryEnd, std::int32_t origQueryStart,
                            std::int32_t origQueryEnd, std::int32_t refPos, bool isReverse);

/// \brief Clip CIGAR to reference (genomic) coordinates.
ClipResult ClipCigarToReference(std::span<const CigarOp> cigar, std::int32_t refStart,
                                std::int32_t refEnd, std::int32_t origRefPos, bool isReverse,
                                bool exciseFlankingInserts = false);

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_CORE_CIGARCLIPPING_HPP
