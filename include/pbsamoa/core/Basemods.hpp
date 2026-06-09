#ifndef PBSAMOA_CORE_BASEMODS_HPP
#define PBSAMOA_CORE_BASEMODS_HPP

#include <span>
#include <string>
#include <string_view>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief A single modification type parsed from an MM tag.
///
/// E.g. "C+m?,1,3,0" -> {Prefix = "C+m?", Skips = {1, 3, 0}}. Each skip is the
/// number of canonical bases to advance before the corresponding modified base.
struct BasemodRecord
{
    std::string Prefix;               ///< e.g. "C+m", "C+m?", "A+a."
    std::vector<std::int32_t> Skips;  ///< canonical-base skip counts
};

/// \brief Parse an MM tag string into ordered modification records.
///
/// "C+m,1,3,0;A+a,2;" -> [{"C+m", {1,3,0}}, {"A+a", {2}}].
std::vector<BasemodRecord> ParseBasemodString(std::string_view mm);

/// \brief Serialize modification records back into an MM tag string.
std::string WriteBasemodString(std::span<const BasemodRecord> records);

/// \brief Number of modification codes carried by a record's prefix: the count of
/// letter codes (e.g. "C+mh" -> 2) or 1 for a numeric ChEBI code (e.g. "C+76792").
///
/// This is the ML stride — ML stores this many quality values per modification site,
/// interleaved per site (spec SAMtags: "C+mh,5,12" -> ML m@0,h@0,m@1,h@1).
std::size_t ModCodeCount(std::string_view prefix);

/// \brief Partition of one record's modification sites against a query clip
/// window [clipOffset, clipOffset + clipLength).
struct BasemodClipWindow
{
    std::size_t FrontRemoved{0};              ///< sites entirely before the window
    std::size_t Retained{0};                  ///< sites within the window
    std::int32_t PrefixLost{0};               ///< canonical bases between the last front-removed
                                              ///< site and the clip start; must be added back to
                                              ///< the first retained skip to undo the clip
    std::vector<std::int32_t> RetainedSkips;  ///< rewritten skips for the retained sites
};

/// \brief Compute how a record's modification sites fall relative to a clip
/// window over the given (pre-clip) sequence.
///
/// Uses the same prefix-sum algorithm as pbbam's ClipBasemodsTag. The trailing
/// site count is `Skips.size() - FrontRemoved - Retained`; the ML quals for a
/// record are ordered [front | retained | trailing] in the original ML array.
///
/// MM/ML coordinates are always in original (5'->3') read orientation, whereas the
/// stored SEQ (and therefore clipOffset/clipLength) are in aligned orientation. For a
/// reverse-strand read (\p reverse true) the canonical base is complemented and the
/// clip window is mirrored onto the original orientation, matching htslib's seqi_rc walk.
BasemodClipWindow ClipBasemodRecord(const BasemodRecord& record, std::string_view sequence,
                                    std::size_t clipOffset, std::size_t clipLength,
                                    bool reverse = false);

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_CORE_BASEMODS_HPP
