#ifndef PBSAMOA_CORE_SAMHEADER_HPP
#define PBSAMOA_CORE_SAMHEADER_HPP

#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
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
    const std::vector<std::pair<std::string, std::string>>& CustomTags() const;

private:
    std::string name_;
    std::int32_t length_;
    std::vector<std::pair<std::string, std::string>> customTags_;
};

/// \brief A single @RG read group entry.
///
/// Required field (ID) stored directly. All optional fields
/// (SM, PL, LB, PU, CN, DS, DT, FO, KS, PG, PI, PM, BC) in tag map.
class ReadGroup
{
public:
    explicit ReadGroup(std::string id);

    std::string_view Id() const;

    const std::string* GetTag(std::string_view key) const;
    void SetTag(std::string key, std::string value);
    const std::vector<std::pair<std::string, std::string>>& CustomTags() const;

private:
    std::string id_;
    std::vector<std::pair<std::string, std::string>> customTags_;
};

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
    const std::vector<std::pair<std::string, std::string>>& CustomTags() const;

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
    /// \throws std::runtime_error on required fields missing from @SQ/@RG/@PG
    static SamHeader FromText(std::string_view text);

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
    ///
    /// \throws std::runtime_error on bad magic, truncated data
    static SamHeader FromBamHeaderBlock(std::span<const std::byte> data);

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
    const std::vector<ReferenceSequence>& ReferenceSequences() const;
    std::vector<ReferenceSequence>& ReferenceSequences();
    void AddReferenceSequence(ReferenceSequence seq);

    // --- @RG ---
    const std::vector<ReadGroup>& ReadGroups() const;
    std::vector<ReadGroup>& ReadGroups();
    void AddReadGroup(ReadGroup rg);

    // --- @PG ---
    const std::vector<ProgramRecord>& ProgramRecords() const;
    std::vector<ProgramRecord>& ProgramRecords();
    void AddProgramRecord(ProgramRecord pg);

    // --- @CO ---
    const std::vector<std::string>& Comments() const;
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
