#ifndef PBSAMOA_CORE_BAMRECORD_HPP
#define PBSAMOA_CORE_BAMRECORD_HPP

#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <string>
#include <string_view>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {

class TagClipper;

/// \brief Type of coordinate clipping.
enum class ClipType
{
    CLIP_TO_QUERY,      ///< clip to polymerase/ZMW coordinates
    CLIP_TO_REFERENCE,  ///< clip to genomic/reference coordinates
};

/// \brief Owned, mutable BAM alignment record with structured fields.
///
/// Fields are stored as high-level types (strings, vectors) — not BAM binary
/// layout. Mutations are plain field assignments. Serialization to BAM binary
/// happens once in the writer via SerializeToBam().
class BamRecord
{
public:
    // --- accessors ---

    std::string_view Name() const;
    std::uint16_t Flag() const;
    std::int32_t RefId() const;
    std::int32_t Pos() const;
    std::uint8_t MapQ() const;
    const std::vector<CigarOp>& Cigar() const;
    std::int32_t NextRefId() const;
    std::int32_t NextPos() const;
    std::int32_t Tlen() const;
    std::string_view Sequence() const;
    const std::vector<std::uint8_t>& Qualities() const;
    const TagMap& Tags() const;

    // --- mutators (fluent interface) ---

    BamRecord& Name(std::string name);
    BamRecord& Flag(std::uint16_t flag);
    BamRecord& RefId(std::int32_t refId);
    BamRecord& Pos(std::int32_t pos);
    BamRecord& MapQ(std::uint8_t mapq);
    BamRecord& Cigar(std::vector<CigarOp> cigar);
    BamRecord& NextRefId(std::int32_t nextRefId);
    BamRecord& NextPos(std::int32_t nextPos);
    BamRecord& Tlen(std::int32_t tlen);
    BamRecord& Sequence(std::string sequence);
    BamRecord& Qualities(std::vector<std::uint8_t> qualities);
    BamRecord& Tags(TagMap tags);

    // --- derived fields ---

    bool IsMapped() const;
    bool IsReverseStrand() const;
    bool IsPrimary() const;
    std::int32_t ReferenceEnd() const;

    /// \brief Serialize to BAM binary layout (excludes block_size prefix).
    std::vector<std::byte> SerializeToBam() const;

    // --- clipping ---

    /// \brief Clip record in-place to the given coordinate range.
    BamRecord& Clip(ClipType type, std::int32_t start, std::int32_t end,
                    bool exciseFlankingInserts = false);
    BamRecord& Clip(ClipType type, std::int32_t start, std::int32_t end, const TagClipper& clipper,
                    bool exciseFlankingInserts = false);

    /// \brief Return a clipped copy of this record.
    BamRecord Clipped(ClipType type, std::int32_t start, std::int32_t end,
                      bool exciseFlankingInserts = false) const;
    BamRecord Clipped(ClipType type, std::int32_t start, std::int32_t end,
                      const TagClipper& clipper, bool exciseFlankingInserts = false) const;

    /// \brief Mutable access to the tag map.
    TagMap& MutableTags();

private:
    std::string name_;
    std::uint16_t flag_{0};
    std::int32_t refId_{-1};
    std::int32_t pos_{-1};
    std::uint8_t mapQ_{0};
    std::vector<CigarOp> cigar_;
    std::int32_t nextRefId_{-1};
    std::int32_t nextPos_{-1};
    std::int32_t tlen_{0};
    std::string sequence_;
    std::vector<std::uint8_t> qualities_;
    TagMap tags_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_CORE_BAMRECORD_HPP
