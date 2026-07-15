#ifndef PBSAMOA_TOOLS_SAMOUTPUT_HPP
#define PBSAMOA_TOOLS_SAMOUTPUT_HPP

#include "../BinaryUtils.hpp"
#include "../SamFieldUtils.hpp"

#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/core/Sequence.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <print>
#include <span>
#include <string>

#include <cstdint>

namespace PacBio {
namespace Samoa {

inline void WriteReferenceName(const SamHeader& header, std::int32_t refId)
{
    const std::string_view sentinel{RnameSentinel(refId)};
    std::print("{}", sentinel.empty() ? header.ReferenceName(refId) : sentinel);
}

inline void WriteNextReferenceName(const SamHeader& header, std::int32_t refId,
                                   std::int32_t nextRefId)
{
    const std::string_view sentinel{RnextSentinel(refId, nextRefId)};
    std::print("{}", sentinel.empty() ? header.ReferenceName(nextRefId) : sentinel);
}

inline void WriteQualities(std::span<const std::uint8_t> qualities)
{
    if (IsQualityUnavailable(qualities)) {
        std::print("*");
        return;
    }

    std::string encodedQualities(std::size(qualities), '\0');
    std::ranges::transform(qualities, std::begin(encodedQualities),
                           [](std::uint8_t q) { return static_cast<char>(q + 33); });
    std::print("{}", encodedQualities);
}

inline void WriteViewAsSam(const SamHeader& header, const RawRecord& view)
{
    std::print("{}\t{}\t", view.Name(), view.Flag());

    const std::int32_t refId{view.RefId()};
    WriteReferenceName(header, refId);
    std::print("\t");

    const std::int32_t pos{view.Pos()};
    std::print("{}\t", OneBasedPositionOrZero(pos));

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
    std::print("{}\t", OneBasedPositionOrZero(nextPos));

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
