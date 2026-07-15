#ifndef PBSAMOA_CORE_BAMRECORD_HPP
#define PBSAMOA_CORE_BAMRECORD_HPP

#include <pbsamoa/core/CigarOp.hpp>
#include <pbsamoa/core/TagClipping.hpp>
#include <pbsamoa/core/Tags.hpp>

#include <pbcopper/data/Accuracy.h>
#include <pbcopper/data/Frames.h>
#include <pbcopper/data/LocalContextFlags.h>
#include <pbcopper/data/SNR.h>

#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include <cstdint>

namespace PacBio {
namespace Samoa {

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
    CigarView Cigar() const;
    std::int32_t NextRefId() const;
    std::int32_t NextPos() const;
    std::int32_t Tlen() const;
    std::string_view Sequence() const;
    std::span<const std::uint8_t> Qualities() const;
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
    bool IsSecondary() const;
    bool IsSupplementary() const;
    std::int32_t ReferenceEnd() const;

    /// \brief First query position consumed by aligned ops, in polymerase
    /// coordinates (QueryStart + strand-appropriate leading soft-clip). Returns
    /// -1 if a hard-clip indicates the query interval is not fully recoverable.
    std::int32_t AlignedStart() const;

    /// \brief One-past-last query position consumed by aligned ops, in polymerase
    /// coordinates (QueryEnd - strand-appropriate trailing soft-clip). Returns -1
    /// if a hard-clip indicates the query interval is not fully recoverable.
    std::int32_t AlignedEnd() const;

    /// \brief Turn an unmapped record into a mapped one.
    ///
    /// Clears the unmapped bit, sets reference id/position/mapq/CIGAR, and sets
    /// the reverse-strand bit per \p reverse. On the reverse strand the SEQ is
    /// reverse-complemented and QUAL is reversed in place, matching the BAM
    /// convention that aligned records store SEQ/QUAL in alignment orientation.
    /// Per-base tags (ip, pw, ML, …) are left in native orientation, matching
    /// pbbam. Precondition: the record currently stores forward-strand SEQ/QUAL
    /// (i.e. it is not already reverse-mapped).
    /// \param[in] refId    reference sequence index
    /// \param[in] pos      0-based reference start position
    /// \param[in] reverse  true if the alignment is on the reverse strand
    /// \param[in] cigar    alignment CIGAR (reference orientation)
    /// \param[in] mapQ     mapping quality
    BamRecord& Map(std::int32_t refId, std::int32_t pos, bool reverse, std::vector<CigarOp> cigar,
                   std::uint8_t mapQ);

    /// \brief Serialize to BAM binary layout (excludes block_size prefix).
    [[nodiscard]] std::vector<std::byte> SerializeToBam() const;

    // --- clipping ---

    /// \brief Clip record in-place to the given coordinate range.
    BamRecord& Clip(ClipType type, std::int32_t start, std::int32_t end,
                    bool exciseFlankingInserts = false);
    BamRecord& Clip(ClipType type, std::int32_t start, std::int32_t end, const TagClipper& clipper,
                    bool exciseFlankingInserts = false);

    /// \brief Return a clipped copy of this record.
    [[nodiscard]] BamRecord Clipped(ClipType type, std::int32_t start, std::int32_t end,
                                    bool exciseFlankingInserts = false) const;
    [[nodiscard]] BamRecord Clipped(ClipType type, std::int32_t start, std::int32_t end,
                                    const TagClipper& clipper,
                                    bool exciseFlankingInserts = false) const;

    /// \brief Mutable access to the tag map.
    TagMap& MutableTags();

    // --- PacBio BAM accessors ---

    /// \brief Full read name (e.g. "movie/123/0_1000").
    std::string FullName() const;

    /// \brief Movie name extracted from the read name prefix before the first '/'.
    std::string MovieName() const;

    /// \brief ZMW hole number extracted from the read name.
    std::int32_t HoleNumber() const;

    /// \brief Query start position. Uses 'qs' tag if present, otherwise parses from read name.
    std::int32_t QueryStart() const;

    /// \brief Query end position. Uses 'qe' tag if present, otherwise parses from read name.
    std::int32_t QueryEnd() const;

    /// \brief Read group ID from the 'RG' auxiliary tag.
    std::string ReadGroupId() const;

    /// \brief Local context flags from the 'cx' tag. Returns nullopt if absent.
    [[nodiscard]] std::optional<Data::LocalContextFlags> LocalContextFlags() const;

    /// \brief Set local context flags in the 'cx' tag.
    void LocalContextFlags(Data::LocalContextFlags flags);

    /// \brief Signal-to-noise ratios from the 'sn' tag (4-element float array).
    Data::SNR SignalToNoise() const;

    /// \brief Read accuracy from the 'rq' tag.
    Data::Accuracy ReadAccuracy() const;

    /// \brief Pulse widths from the 'pw' tag. Returns nullopt if absent.
    [[nodiscard]] std::optional<Data::Frames> PulseWidth() const;

    /// \brief Inter-pulse durations from the 'ip' tag. Returns nullopt if absent.
    [[nodiscard]] std::optional<Data::Frames> IPD() const;

    /// \brief Wall-clock start time from the 'ws' tag. Returns nullopt if absent.
    [[nodiscard]] std::optional<std::int32_t> WallStart() const;

    /// \brief Wall-clock end time from the 'we' tag. Returns nullopt if absent.
    [[nodiscard]] std::optional<std::int32_t> WallEnd() const;

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
