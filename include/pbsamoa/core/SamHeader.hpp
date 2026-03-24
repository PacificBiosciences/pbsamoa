#ifndef PBSAMOA_CORE_SAMHEADER_HPP
#define PBSAMOA_CORE_SAMHEADER_HPP

#include <pbcopper/data/Strand.h>

#include <expected>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief A single @SQ reference sequence entry.
///
/// Required fields (SN, LN) are stored directly. All optional fields
/// (AS, M5, SP, TP, UR, AN, AH, DS) stored in a generic tag map.
class ReferenceSequence
{
public:
    ReferenceSequence(std::string name, std::int32_t length);

    std::string_view Name() const;
    std::int32_t Length() const;

    const std::string* GetTag(std::string_view key) const;
    void SetTag(std::string key, std::string value);
    std::span<const std::pair<std::string, std::string>> CustomTags() const;

private:
    std::string name_;
    std::int32_t length_;
    std::vector<std::pair<std::string, std::string>> customTags_;
};

/// \brief Parsed key=value fields from a PacBio @RG DS tag.
using DsTagFields = std::vector<std::pair<std::string, std::string>>;

/// \brief A single @RG read group entry.
///
/// Required field (ID) stored directly. All optional fields
/// (SM, PL, LB, PU, CN, DS, DT, FO, KS, PG, PI, PM, BC) in tag map.
/// PacBio-specific DS tag fields are accessible via named accessors.
class ReadGroup
{
public:
    explicit ReadGroup(std::string id);

    std::string_view Id() const;

    // --- Generic tag access ---

    const std::string* GetTag(std::string_view key) const;
    void SetTag(std::string key, std::string value);
    std::span<const std::pair<std::string, std::string>> CustomTags() const;

    // --- PacBio convenience accessors ---

    /// \brief Movie name from the PU tag.
    const std::string* MovieName() const;

    // --- PacBio DS tag field accessors ---

    /// \brief Access a DS field by key. Returns nullptr if DS tag absent or key not found.
    const std::string* DsField(std::string_view key) const;

    /// \brief Set a DS field, creating or updating the DS tag as needed.
    void SetDsField(std::string key, std::string value);

    /// \brief Parse the full DS tag into key-value fields.
    DsTagFields ParsedDsFields() const;

    /// \brief Replace the entire DS tag from a DsTagFields collection.
    void SetDsFields(std::span<const std::pair<std::string, std::string>> fields);

    /// \brief Read type (e.g. "SUBREAD", "CCS").
    const std::string* ReadType() const;
    /// \brief Binding kit part number.
    const std::string* BindingKit() const;
    /// \brief Sequencing kit part number.
    const std::string* SequencingKit() const;
    /// \brief Basecaller version string.
    const std::string* BasecallerVersion() const;

private:
    std::string id_;
    std::vector<std::pair<std::string, std::string>> customTags_;

    /// \brief Lazily parsed DS tag fields cache.
    mutable std::optional<DsTagFields> parsedDs_;

    const DsTagFields& EnsureParsedDs() const;
    void InvalidateDsCache();
};

// --- PacBio read group ID generation ---

/// \brief Generate a PacBio read group ID from movie name, read type, and optional strand.
///
/// Computes MD5("{movieName}//{readType}[//fwd|//rev]") and returns the first
/// 8 hex characters.
std::string MakeReadGroupId(std::string_view movieName, std::string_view readType,
                            std::optional<Data::Strand> strand = std::nullopt);

/// \brief Generate a barcoded PacBio read group ID.
///
/// Returns MakeReadGroupId(movieName, readType, strand) + "/" + barcodeSuffix.
std::string MakeReadGroupId(std::string_view movieName, std::string_view readType,
                            std::string_view barcodeSuffix,
                            std::optional<Data::Strand> strand = std::nullopt);

/// \brief Strip barcode suffix from a read group ID, returning the base hash.
std::string_view ReadGroupBaseId(std::string_view readGroupId);

/// \brief A single @PG program record entry.
///
/// Required field (ID) stored directly. Optional fields
/// (PN, VN, CL, PP, DS) in tag map. PP enables chaining.
class ProgramRecord
{
public:
    explicit ProgramRecord(std::string id);

    std::string_view Id() const;

    const std::string* GetTag(std::string_view key) const;
    void SetTag(std::string key, std::string value);
    std::span<const std::pair<std::string, std::string>> CustomTags() const;

private:
    std::string id_;
    std::vector<std::pair<std::string, std::string>> customTags_;
};

/// \brief Parsed SAM/BAM header.
///
/// Holds @HD metadata, @SQ reference sequences (in file order, which defines
/// refID mapping), @RG read groups, @PG program records, and @CO comments.
class SamHeader
{
public:
    SamHeader();

    /// \brief Parse header from SAM text (one or more lines starting with @).
    static std::expected<SamHeader, std::string> FromText(std::string_view text);

    /// \brief Serialize header to SAM text.
    /// Line order: @HD, @SQ (in insertion order), @RG, @PG, @CO.
    std::string ToText() const;

    /// \brief Parse header from BAM binary header block.
    ///
    /// The span covers everything from the BAM magic through the end of the
    /// reference dictionary: magic(4) + l_text(4) + text(l_text) + n_ref(4)
    /// + [l_name(4) + name(l_name) + l_ref(4)] * n_ref.
    ///
    /// The binary reference dictionary is authoritative. If @SQ lines are
    /// present in header text, they are parsed for optional tags only.
    static std::expected<SamHeader, std::string> FromBamHeaderBlock(
        std::span<const std::byte> data);

    /// \brief Serialize to BAM binary header block.
    ///
    /// Produces: magic(4) + l_text(4) + text(l_text) + n_ref(4) +
    ///           [l_name(4) + name(l_name) + l_ref(4)] * n_ref
    std::vector<std::byte> ToBamHeaderBlock() const;

    // --- @HD fields ---
    std::string_view Version() const;
    std::string_view SortOrder() const;
    std::string_view GroupOrder() const;
    std::string_view SubSort() const;

    // --- @HD mutators ---
    void SetVersion(std::string version);
    void SetSortOrder(std::string sortOrder);
    void SetGroupOrder(std::string groupOrder);
    void SetSubSort(std::string subSort);

    // --- @SQ ---
    std::span<const ReferenceSequence> ReferenceSequences() const;
    std::vector<ReferenceSequence>& ReferenceSequences();
    void AddReferenceSequence(ReferenceSequence seq);

    // --- @RG ---
    std::span<const ReadGroup> ReadGroups() const;
    std::vector<ReadGroup>& ReadGroups();
    void AddReadGroup(ReadGroup rg);

    // --- @PG ---
    std::span<const ProgramRecord> ProgramRecords() const;
    std::vector<ProgramRecord>& ProgramRecords();
    void AddProgramRecord(ProgramRecord pg);

    // --- @CO ---
    std::span<const std::string> Comments() const;
    void AddComment(std::string comment);

    // --- Reference name/ID lookup ---

    /// \brief Look up reference sequence ID by name. Returns -1 if not found.
    std::int32_t ReferenceId(std::string_view name) const;

    /// \brief Look up reference sequence name by ID. Returns "*" for -1.
    /// \throws std::out_of_range for invalid IDs
    std::string_view ReferenceName(std::int32_t refId) const;

    /// \brief Number of reference sequences.
    std::int32_t NumReferences() const;

private:
    // @HD
    std::string version_;
    std::string sortOrder_;
    std::string groupOrder_;
    std::string subSort_;

    std::vector<ReferenceSequence> references_;
    std::vector<ReadGroup> readGroups_;
    std::vector<ProgramRecord> programRecords_;
    std::vector<std::string> comments_;

    void BuildNameIndex() const;

    struct TransparentStringHash
    {
        using is_transparent = void;

        std::size_t operator()(std::string_view sv) const noexcept
        {
            return std::hash<std::string_view>{}(sv);
        }
    };

    mutable std::unordered_map<std::string, std::int32_t, TransparentStringHash, std::equal_to<>>
        nameToId_;
    mutable std::size_t nameIndexSize_{0};  // tracks reference count at last rebuild
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_CORE_SAMHEADER_HPP
