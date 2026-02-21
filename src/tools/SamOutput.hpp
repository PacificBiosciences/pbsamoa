#ifndef PBSAMOA_TOOLS_SAMOUTPUT_HPP
#define PBSAMOA_TOOLS_SAMOUTPUT_HPP

#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/RawRecord.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/core/Sequence.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <algorithm>
#include <span>
#include <string>

#include <cstdint>
#include <cstdio>

namespace PacBio {
namespace Samoa {

inline void WriteViewAsSam(const SamHeader& header, const RawRecord& view)
{
    std::fputs(std::string{view.Name()}.c_str(), stdout);
    std::fputc('\t', stdout);

    std::fprintf(stdout, "%u", view.Flag());
    std::fputc('\t', stdout);

    const std::int32_t refId{view.RefId()};
    if (refId < 0) {
        std::fputc('*', stdout);
    } else {
        std::fputs(std::string{header.ReferenceName(refId)}.c_str(), stdout);
    }
    std::fputc('\t', stdout);

    const std::int32_t pos{view.Pos()};
    std::fprintf(stdout, "%d", pos >= 0 ? pos + 1 : 0);
    std::fputc('\t', stdout);

    std::fprintf(stdout, "%u", view.MapQ());
    std::fputc('\t', stdout);

    const CigarView cigar{view.CigarOps()};
    if (std::empty(cigar)) {
        std::fputc('*', stdout);
    } else {
        std::fputs(CigarToString(cigar).c_str(), stdout);
    }
    std::fputc('\t', stdout);

    const std::int32_t nextRefId{view.NextRefId()};
    if (nextRefId < 0) {
        std::fputc('*', stdout);
    } else if (nextRefId == refId) {
        std::fputc('=', stdout);
    } else {
        std::fputs(std::string{header.ReferenceName(nextRefId)}.c_str(), stdout);
    }
    std::fputc('\t', stdout);

    const std::int32_t nextPos{view.NextPos()};
    std::fprintf(stdout, "%d", nextPos >= 0 ? nextPos + 1 : 0);
    std::fputc('\t', stdout);

    std::fprintf(stdout, "%d", view.Tlen());
    std::fputc('\t', stdout);

    const std::string seq{view.Seq().ToString()};
    if (std::empty(seq)) {
        std::fputc('*', stdout);
    } else {
        std::fputs(seq.c_str(), stdout);
    }
    std::fputc('\t', stdout);

    const std::span<const std::uint8_t> qual{view.Qual()};
    const bool qualUnavailable{std::empty(qual) ||
                               std::ranges::all_of(qual, [](std::uint8_t q) { return q == 0xFF; })};
    if (qualUnavailable) {
        std::fputc('*', stdout);
    } else {
        for (const std::uint8_t q : qual) {
            std::fputc(static_cast<unsigned char>(q + 33), stdout);
        }
    }

    std::string tagBuf;
    SerializeRawTagsToSam(view.AuxData(), tagBuf);
    if (!std::empty(tagBuf)) {
        std::fputs(tagBuf.c_str(), stdout);
    }
    std::fputc('\n', stdout);
}

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TOOLS_SAMOUTPUT_HPP
