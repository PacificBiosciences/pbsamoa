#ifndef PBSAMOA_TOOLS_SAMOUTPUT_HPP
#define PBSAMOA_TOOLS_SAMOUTPUT_HPP

#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/core/Sequence.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <algorithm>
#include <print>
#include <span>
#include <string>

#include <cstdint>

namespace PacBio {
namespace Samoa {

inline std::int32_t OneBasedPosOrZero(std::int32_t pos)
{
    if (pos < 0) {
        return 0;
    }
    return pos + 1;
}

inline void WriteReferenceName(const SamHeader& header, std::int32_t refId)
{
    if (refId < 0) {
        std::print("*");
        return;
    }
    std::print("{}", header.ReferenceName(refId));
}

inline void WriteNextReferenceName(const SamHeader& header, std::int32_t refId,
                                   std::int32_t nextRefId)
{
    if (nextRefId < 0) {
        std::print("*");
        return;
    }
    if (nextRefId == refId) {
        std::print("=");
        return;
    }
    std::print("{}", header.ReferenceName(nextRefId));
}

inline bool IsQualityUnavailable(std::span<const std::uint8_t> qualities)
{
    return std::empty(qualities) ||
           std::ranges::all_of(qualities, [](std::uint8_t quality) { return quality == 0xFF; });
}

inline void WriteQualities(std::span<const std::uint8_t> qualities)
{
    if (IsQualityUnavailable(qualities)) {
        std::print("*");
        return;
    }

    std::string encodedQualities;
    encodedQualities.reserve(std::size(qualities));
    for (const std::uint8_t quality : qualities) {
        encodedQualities.push_back(static_cast<char>(quality + 33));
    }
    std::print("{}", encodedQualities);
}

inline void WriteViewAsSam(const SamHeader& header, const RawRecord& view)
{
    std::print("{}\t{}\t", view.Name(), view.Flag());

    const std::int32_t refId{view.RefId()};
    WriteReferenceName(header, refId);
    std::print("\t");

    const std::int32_t pos{view.Pos()};
    std::print("{}\t", OneBasedPosOrZero(pos));

    std::print("{}\t", view.MapQ());

    const CigarView cigar{view.CigarOps()};
    if (std::empty(cigar)) {
        std::print("*");
    } else {
        std::print("{}", CigarToString(cigar));
    }
    std::print("\t");

    const std::int32_t nextRefId{view.NextRefId()};
    WriteNextReferenceName(header, refId, nextRefId);
    std::print("\t");

    const std::int32_t nextPos{view.NextPos()};
    std::print("{}\t", OneBasedPosOrZero(nextPos));

    std::print("{}\t", view.Tlen());

    const std::string seq{view.Seq().ToString()};
    if (std::empty(seq)) {
        std::print("*");
    } else {
        std::print("{}", seq);
    }
    std::print("\t");

    WriteQualities(view.Qual());

    std::string tagBuf;
    SerializeRawTagsToSam(view.AuxData(), tagBuf);
    if (!std::empty(tagBuf)) {
        std::print("{}", tagBuf);
    }
    std::println();
}

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_SAMOUTPUT_HPP
