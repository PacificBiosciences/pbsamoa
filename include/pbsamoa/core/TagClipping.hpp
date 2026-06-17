#ifndef PBSAMOA_CORE_TAGCLIPPING_HPP
#define PBSAMOA_CORE_TAGCLIPPING_HPP

#include <pbsamoa/core/Tags.hpp>

#include <functional>
#include <string>
#include <string_view>
#include <vector>

#include <cstddef>

namespace PacBio {
namespace Samoa {

/// \brief Context passed to tag clipping strategies for cross-tag access.
struct ClipContext
{
    const TagMap* tags{nullptr};

    /// \brief Snapshot of the original pulse-call string (the 'pc' tag)
    /// captured before the clipping loop. Used by PulseClipStrategy so
    /// that all pulse tags see the unmodified pulse-call data regardless
    /// of iteration order.
    std::string pulseCalls;

    /// \brief Original (pre-clip) sequence for base-modification clipping.
    /// Empty when not supplied (basemods clipping is skipped).
    std::string_view sequence;

    /// \brief Snapshot of the original MM tag value captured before the
    /// clipping loop. Used by BasemodClipStrategy when processing the ML
    /// tag so it can reconstruct the pre-clip modification layout.
    std::string basemodString;

    /// \brief True when the read is reverse-strand (FLAG 0x10). MM/ML
    /// coordinates are in original read orientation while the stored SEQ
    /// is reverse-complemented, so basemod clipping must mirror the window
    /// and complement the canonical base.
    bool isReverse{false};
};

/// \brief Abstract interface for a tag clipping strategy.
class TagClipStrategy
{
public:
    TagClipStrategy() = default;
    virtual ~TagClipStrategy() = default;
    TagClipStrategy(const TagClipStrategy&) = delete;
    TagClipStrategy& operator=(const TagClipStrategy&) = delete;
    TagClipStrategy(TagClipStrategy&&) = delete;
    TagClipStrategy& operator=(TagClipStrategy&&) = delete;

    /// \brief Clip a tag value in-place.
    /// \param[in,out] value      the tag value to clip
    /// \param[in] clipOffset     query-base offset into original data
    /// \param[in] clipLength     number of query bases retained
    /// \param[in] seqLength      original (pre-clip) sequence length
    /// \param[in] ctx            optional context for cross-tag access
    /// \returns false if the tag should be removed (malformed data)
    virtual bool Clip(TagValue& value, std::size_t clipOffset, std::size_t clipLength,
                      std::size_t seqLength, const ClipContext& ctx = {}) const = 0;
};

/// \brief Clips string and array tags by simple substring [offset, offset+length).
/// Used for: dq, iq, mq, sq, dt, st, ip, pw, fi, fp
class SubstringClipStrategy : public TagClipStrategy
{
public:
    bool Clip(TagValue& value, std::size_t clipOffset, std::size_t clipLength,
              std::size_t seqLength, const ClipContext& ctx = {}) const override;
};

/// \brief Clips reverse-orientation tags by mirrored offset.
/// reverseOffset = seqLength - (clipOffset + clipLength), then substring.
/// Used for: ri, rp
class ReverseSubstringClipStrategy : public TagClipStrategy
{
public:
    bool Clip(TagValue& value, std::size_t clipOffset, std::size_t clipLength,
              std::size_t seqLength, const ClipContext& ctx = {}) const override;
};

/// \brief Clips pulse-space tags via PulseBitset mapping.
/// Needs the 'pc' tag from ClipContext to build the bitset.
/// Used for: pc, pt, pq, pv, pg, pa, pm, ps, pi, pd, px, pe, sf
class PulseClipStrategy : public TagClipStrategy
{
public:
    bool Clip(TagValue& value, std::size_t clipOffset, std::size_t clipLength,
              std::size_t seqLength, const ClipContext& ctx = {}) const override;
};

/// \brief Clips MM and ML base-modification tags per SAMv1.6.
///
/// For MM (string): parses modification records, walks the original
/// sequence to find canonical-base positions, determines which
/// modification sites fall within the clip window, and rewrites the
/// skip counts.
///
/// For ML (B:C array): uses the snapshotted original MM string and
/// sequence to compute how many modification sites precede the clip
/// window and how many are retained, then slices the quality array.
///
/// Requires ClipContext::sequence to be non-empty; returns false
/// (tag removed) when the sequence is unavailable.
class BasemodClipStrategy : public TagClipStrategy
{
public:
    bool Clip(TagValue& value, std::size_t clipOffset, std::size_t clipLength,
              std::size_t seqLength, const ClipContext& ctx = {}) const override;
};

/// \brief Clips the run-length-encoded sa (subread pileup coverage) tag.
///
/// The sa tag is a B:C array with alternating pairs: [runLen1, cov1, runLen2, cov2, ...].
/// The sum of all run lengths equals the sequence length. This strategy
/// clips the RLE structure to the query window [clipOffset, clipOffset+clipLength),
/// splitting runs at the boundaries as needed.
class PileupClipStrategy : public TagClipStrategy
{
public:
    bool Clip(TagValue& value, std::size_t clipOffset, std::size_t clipLength,
              std::size_t seqLength, const ClipContext& ctx = {}) const override;
};

/// \brief Registry mapping tag keys to clipping strategies.
class TagClipper
{
public:
    /// \brief Register a strategy for a set of tags.
    void Register(std::initializer_list<TagKey> tags, const TagClipStrategy& strategy);

    /// \brief Clip all registered tags in the TagMap.
    /// Tags not in the registry are left untouched.
    /// \param[in] sequence  original (pre-clip) sequence for basemods;
    ///                      empty by default (basemods clipping skipped)
    /// \param[in] isReverse true if the read is reverse-strand (FLAG 0x10), so
    ///                      MM/ML clipping mirrors the window and complements the base
    void ClipTags(TagMap& tags, std::size_t clipOffset, std::size_t clipLength,
                  std::size_t seqLength, std::string_view sequence = {},
                  bool isReverse = false) const;

    /// \brief Create a TagClipper with all PacBio-standard strategies registered.
    static TagClipper PacBioDefault();

private:
    struct Registration
    {
        TagKey key;
        std::reference_wrapper<const TagClipStrategy> strategy;
    };

    std::vector<Registration> registrations_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_CORE_TAGCLIPPING_HPP
