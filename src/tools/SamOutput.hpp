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

inline void WriteViewAsSam(const SamHeader& header, const RawRecord& view)
{
    std::print("{}\t{}\t", view.Name(), view.Flag());

    const std::int32_t refId{view.RefId()};
    if (refId < 0) {
        std::print("*");
    } else {
        std::print("{}", header.ReferenceName(refId));
    }
    std::print("\t");

    const std::int32_t pos{view.Pos()};
    std::print("{}\t", pos >= 0 ? pos + 1 : 0);

    std::print("{}\t", view.MapQ());

    const CigarView cigar{view.CigarOps()};
    if (std::empty(cigar)) {
        std::print("*");
    } else {
        std::print("{}", CigarToString(cigar));
    }
    std::print("\t");

    const std::int32_t nextRefId{view.NextRefId()};
    if (nextRefId < 0) {
        std::print("*");
    } else if (nextRefId == refId) {
        std::print("=");
    } else {
        std::print("{}", header.ReferenceName(nextRefId));
    }
    std::print("\t");

    const std::int32_t nextPos{view.NextPos()};
    std::print("{}\t", nextPos >= 0 ? nextPos + 1 : 0);

    std::print("{}\t", view.Tlen());

    const std::string seq{view.Seq().ToString()};
    if (std::empty(seq)) {
        std::print("*");
    } else {
        std::print("{}", seq);
    }
    std::print("\t");

    const std::span<const std::uint8_t> qual{view.Qual()};
    const bool qualUnavailable{std::empty(qual) ||
                               std::ranges::all_of(qual, [](std::uint8_t q) { return q == 0xFF; })};
    if (qualUnavailable) {
        std::print("*");
    } else {
        for (const std::uint8_t q : qual) {
            std::print("{}", static_cast<char>(q + 33));
        }
    }

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
