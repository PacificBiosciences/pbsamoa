#include <pbsamoa/core/SamHeader.hpp>

#include "BinaryUtils.hpp"
#include "Md5.hpp"

#include <algorithm>
#include <charconv>
#include <expected>
#include <format>
#include <stdexcept>
#include <unordered_set>
#include <utility>

#include <cctype>
#include <cstddef>

namespace {

/// \brief Split string_view by a delimiter character.
/// Returns vector of non-owning views. Skips empty trailing entries.
std::vector<std::string_view> Split(std::string_view sv, char delim)
{
    std::vector<std::string_view> result;
    while (!std::empty(sv)) {
        const std::size_t pos{sv.find(delim)};
        if (pos == std::string_view::npos) {
            result.push_back(sv);
            break;
        }
        result.push_back(sv.substr(0, pos));
        sv.remove_prefix(pos + 1);
    }
    return result;
}

/// \brief Parse a TAG:VALUE pair. Returns {tag, value}.
/// The tag is the first two characters before the colon.
/// The value is everything after the first colon.
std::pair<std::string_view, std::string_view> ParseTagValue(std::string_view field)
{
    const std::size_t colon{field.find(':')};
    if (colon == std::string_view::npos) {
        return {field, {}};
    }
    return {field.substr(0, colon), field.substr(colon + 1)};
}

template <typename Record>
void MoveCustomTagsToRecord(std::vector<std::pair<std::string, std::string>>& tags, Record& record)
{
    for (auto& [key, value] : tags) {
        record.SetTag(std::move(key), std::move(value));
    }
}

void AppendSamField(std::string& result, std::string_view key, std::string_view value)
{
    result += '\t';
    result += key;
    result += ':';
    result += value;
}

void AppendSamField(std::string& result, std::string_view key, std::int32_t value)
{
    result += '\t';
    result += key;
    result += ':';
    std::format_to(std::back_inserter(result), "{}", value);
}

void AppendCustomTags(std::string& result,
                      std::span<const std::pair<std::string, std::string>> customTags)
{
    for (const auto& [key, value] : customTags) {
        AppendSamField(result, key, value);
    }
}

const std::string* FindValue(std::span<const std::pair<std::string, std::string>> fields,
                             std::string_view key)
{
    const auto it{std::ranges::find(fields, key, &std::pair<std::string, std::string>::first)};
    if (it == std::end(fields)) {
        return nullptr;
    }
    return &it->second;
}

void UpsertValue(std::vector<std::pair<std::string, std::string>>& fields, std::string key,
                 std::string value)
{
    auto it{std::ranges::find(fields, key, &std::pair<std::string, std::string>::first)};
    if (it != std::end(fields)) {
        it->second = std::move(value);
        return;
    }
    fields.emplace_back(std::move(key), std::move(value));
}

std::expected<std::int32_t, std::string> ParseInt32Field(std::string_view value,
                                                         std::string_view errorMessage)
{
    std::int32_t parsed{0};
    const char* const end{std::data(value) + std::size(value)};
    const auto [ptr, ec]{std::from_chars(std::data(value), end, parsed)};
    if ((ec != std::errc{}) || (ptr != end)) {
        return std::unexpected{std::string{errorMessage}};
    }
    return parsed;
}

template <typename Record>
std::expected<Record, std::string> ParseTaggedRecord(std::span<const std::string_view> fields,
                                                     std::string_view requiredTag,
                                                     std::string_view missingFieldMessage)
{
    std::string requiredValue;
    std::vector<std::pair<std::string, std::string>> tags;

    for (const std::string_view field : fields.subspan(1)) {
        const auto [tag, value] = ParseTagValue(field);
        if (tag == requiredTag) {
            requiredValue = std::string{value};
        } else {
            tags.emplace_back(std::string{tag}, std::string{value});
        }
    }

    if (std::empty(requiredValue)) {
        return std::unexpected{std::string{missingFieldMessage}};
    }

    Record record{std::move(requiredValue)};
    MoveCustomTagsToRecord(tags, record);
    return record;
}

std::expected<PacBio::Samoa::ReferenceSequence, std::string> ParseReferenceSequenceRecord(
    std::span<const std::string_view> fields)
{
    std::string name;
    std::optional<std::int32_t> length;
    std::vector<std::pair<std::string, std::string>> tags;

    for (const std::string_view field : fields.subspan(1)) {
        const auto [tag, value] = ParseTagValue(field);
        if (tag == "SN") {
            name = std::string{value};
        } else if (tag == "LN") {
            auto parsedLength{ParseInt32Field(value, "Invalid LN value in @SQ line")};
            if (!parsedLength) {
                return std::unexpected{std::move(parsedLength.error())};
            }
            // htslib (header.c:158-187) tolerates a repeated LN with the same value but
            // rejects conflicting values; a duplicate same-value LN is harmless.
            if (length && (*length != *parsedLength)) {
                return std::unexpected{"@SQ line has multiple LN tags with different values"};
            }
            length = *parsedLength;
        } else {
            tags.emplace_back(std::string{tag}, std::string{value});
        }
    }

    if (std::empty(name)) {
        return std::unexpected{"Missing required SN field in @SQ line"};
    }
    if (!length || (*length < 1)) {
        return std::unexpected{"Missing required LN field in @SQ line"};
    }

    PacBio::Samoa::ReferenceSequence ref{std::move(name), *length};
    MoveCustomTagsToRecord(tags, ref);
    return ref;
}

}  // namespace

namespace PacBio {
namespace Samoa {
// --- ReferenceSequence ---

ReferenceSequence::ReferenceSequence(std::string name, std::int32_t length)
    : name_{std::move(name)}, length_{length}
{
}

std::string_view ReferenceSequence::Name() const { return name_; }

std::int32_t ReferenceSequence::Length() const { return length_; }

const std::string* ReferenceSequence::GetTag(std::string_view key) const
{
    return FindValue(customTags_, key);
}

void ReferenceSequence::SetTag(std::string key, std::string value)
{
    UpsertValue(customTags_, std::move(key), std::move(value));
}

std::span<const std::pair<std::string, std::string>> ReferenceSequence::CustomTags() const
{
    return customTags_;
}

// --- ReadGroup ---

ReadGroup::ReadGroup(std::string id) : id_{std::move(id)} {}

std::string_view ReadGroup::Id() const { return id_; }

const std::string* ReadGroup::GetTag(std::string_view key) const
{
    return FindValue(customTags_, key);
}

void ReadGroup::SetTag(std::string key, std::string value)
{
    if (key == "DS") {
        InvalidateDsCache();
    }
    UpsertValue(customTags_, std::move(key), std::move(value));
}

std::span<const std::pair<std::string, std::string>> ReadGroup::CustomTags() const
{
    return customTags_;
}

// --- PacBio convenience accessors ---

const std::string* ReadGroup::MovieName() const { return GetTag("PU"); }

// --- PacBio DS tag field accessors ---

namespace {

bool IsSpace(char ch) { return std::isspace(static_cast<unsigned char>(ch)) != 0; }

std::string TrimCopy(std::string_view text)
{
    const auto first{
        static_cast<std::size_t>(std::ranges::find_if_not(text, IsSpace) - text.begin())};
    std::size_t last{std::size(text)};
    while ((last > first) && IsSpace(text[last - 1])) {
        --last;
    }
    return std::string{text.substr(first, last - first)};
}

void AppendDsField(DsTagFields& fields, std::string_view token)
{
    const std::size_t equals{token.find('=')};
    const std::string key{TrimCopy(token.substr(0, equals))};
    std::string value;
    if (equals != std::string_view::npos) {
        value = TrimCopy(token.substr(equals + 1));
    }
    fields.emplace_back(key, value);
}

DsTagFields ParseDsTagImpl(std::string_view dsTag)
{
    DsTagFields fields;
    for (const std::string_view token : Split(dsTag, ';')) {
        if (!token.empty()) {
            AppendDsField(fields, token);
        }
    }
    return fields;
}

std::string BuildDsTagImpl(std::span<const std::pair<std::string, std::string>> fields)
{
    std::string result;
    for (std::size_t i{0}; i < std::size(fields); ++i) {
        if (i > 0) {
            result.push_back(';');
        }
        result += fields[i].first;
        result.push_back('=');
        result += fields[i].second;
    }
    return result;
}

}  // namespace

const DsTagFields& ReadGroup::EnsureParsedDs() const
{
    if (!parsedDs_) {
        const std::string* dsTag{GetTag("DS")};
        if (dsTag) {
            parsedDs_.emplace(ParseDsTagImpl(*dsTag));
        } else {
            parsedDs_.emplace();
        }
    }
    return *parsedDs_;
}

void ReadGroup::InvalidateDsCache() { parsedDs_.reset(); }

const std::string* ReadGroup::DsField(std::string_view key) const
{
    return FindValue(EnsureParsedDs(), key);
}

void ReadGroup::SetDsField(std::string key, std::string value)
{
    DsTagFields fields{EnsureParsedDs()};
    UpsertValue(fields, std::move(key), std::move(value));
    SetTag("DS", BuildDsTagImpl(fields));
    parsedDs_ = std::move(fields);
}

DsTagFields ReadGroup::ParsedDsFields() const { return EnsureParsedDs(); }

void ReadGroup::SetDsFields(std::span<const std::pair<std::string, std::string>> fields)
{
    SetTag("DS", BuildDsTagImpl(fields));
}

const std::string* ReadGroup::ReadType() const { return DsField("READTYPE"); }

const std::string* ReadGroup::BindingKit() const { return DsField("BINDINGKIT"); }

const std::string* ReadGroup::SequencingKit() const { return DsField("SEQUENCINGKIT"); }

const std::string* ReadGroup::BasecallerVersion() const { return DsField("BASECALLERVERSION"); }

// --- PacBio read group ID generation ---

namespace {

constexpr std::string_view ForwardSuffix{"//fwd"};
constexpr std::string_view ReverseSuffix{"//rev"};

}  // namespace

std::string MakeReadGroupId(std::string_view movieName, std::string_view readType,
                            std::optional<Data::Strand> strand)
{
    std::string content{std::format("{}//{}", movieName, readType)};
    const Data::Strand resolvedStrand{strand.value_or(Data::Strand::UNMAPPED)};
    if (resolvedStrand == Data::Strand::FORWARD) {
        content += ForwardSuffix;
    } else if (resolvedStrand == Data::Strand::REVERSE) {
        content += ReverseSuffix;
    }

    return detail::Md5Hex(content).substr(0, 8);
}

std::string MakeReadGroupId(std::string_view movieName, std::string_view readType,
                            std::string_view barcodeSuffix, std::optional<Data::Strand> strand)
{
    return MakeReadGroupId(movieName, readType, strand) + "/" + std::string{barcodeSuffix};
}

std::string_view ReadGroupBaseId(std::string_view readGroupId)
{
    return readGroupId.substr(0, readGroupId.find('/'));
}

// --- ProgramRecord ---

ProgramRecord::ProgramRecord(std::string id) : id_{std::move(id)} {}

std::string_view ProgramRecord::Id() const { return id_; }

const std::string* ProgramRecord::GetTag(std::string_view key) const
{
    return FindValue(customTags_, key);
}

void ProgramRecord::SetTag(std::string key, std::string value)
{
    UpsertValue(customTags_, std::move(key), std::move(value));
}

std::span<const std::pair<std::string, std::string>> ProgramRecord::CustomTags() const
{
    return customTags_;
}

// --- SamHeader parsing ---

SamHeader::SamHeader() = default;

std::expected<SamHeader, std::string> SamHeader::FromText(std::string_view text)
{
    SamHeader header;
    // The spec (SAMv1 §1.3) requires @SQ SN names to be distinct; htslib treats a repeat
    // as a fatal parse error (header.c:236-240). Track seen names to enforce the same.
    std::unordered_set<std::string> seenRefNames;
    std::unordered_set<std::string> seenRgIds;
    bool seenHd{false};
    const std::vector<std::string_view> lines{Split(text, '\n')};
    for (const std::string_view line : lines) {
        if (std::empty(line)) {
            continue;
        }
        const std::vector<std::string_view> fields{Split(line, '\t')};
        if (std::empty(fields)) {
            continue;
        }
        const std::string_view recordType{fields[0]};

        if (recordType == "@HD") {
            // SAMv1 §1.3: if present, there must be only one @HD line.
            if (seenHd) {
                return std::unexpected{"Multiple @HD lines in sam header"};
            }
            seenHd = true;
            for (const std::string_view field : std::span{fields}.subspan(1)) {
                const auto [tag, value] = ParseTagValue(field);
                if (tag == "VN") {
                    header.version_ = std::string{value};
                } else if (tag == "SO") {
                    header.sortOrder_ = std::string{value};
                } else if (tag == "GO") {
                    header.groupOrder_ = std::string{value};
                } else if (tag == "SS") {
                    header.subSort_ = std::string{value};
                } else {
                    // Preserve unknown/custom @HD tags so they round-trip, matching htslib.
                    header.hdCustomTags_.emplace_back(std::string{tag}, std::string{value});
                }
            }
            continue;
        }
        if (recordType == "@SQ") {
            auto parsed{ParseReferenceSequenceRecord(fields)};
            if (!parsed) {
                return std::unexpected{std::move(parsed.error())};
            }
            if (!seenRefNames.insert(std::string{parsed->Name()}).second) {
                return std::unexpected{"Duplicate @SQ SN \"" + std::string{parsed->Name()} +
                                       "\" in sam header"};
            }
            header.references_.push_back(std::move(*parsed));
            continue;
        }
        if (recordType == "@RG") {
            // SAMv1 §1.3: each @RG ID must be unique.
            auto parsed{ParseTaggedRecord<ReadGroup>(fields, "ID",
                                                     "Missing required ID field in @RG line")};
            if (!parsed) {
                return std::unexpected{std::move(parsed.error())};
            }
            if (!seenRgIds.insert(std::string{parsed->Id()}).second) {
                return std::unexpected{"Duplicate @RG ID \"" + std::string{parsed->Id()} +
                                       "\" in sam header"};
            }
            header.readGroups_.push_back(std::move(*parsed));
            continue;
        }
        if (recordType == "@PG") {
            // SAMv1 §1.3 requires unique @PG IDs, but a duplicate does not stop a read.
            // @PG IDs are opaque; only PP refers to them, and pbsamoa never resolves PP.
            // Other tools and pbsamoa builds before the UniqueProgramId fix emit
            // duplicates, so keep every line instead of rejecting the whole file. htslib
            // keeps them too (header.c). pbsamoa's own writers keep new IDs unique.
            auto parsed{ParseTaggedRecord<ProgramRecord>(fields, "ID",
                                                         "Missing required ID field in @PG line")};
            if (!parsed) {
                return std::unexpected{std::move(parsed.error())};
            }
            header.programRecords_.push_back(std::move(*parsed));
            continue;
        }
        if (recordType == "@CO") {
            // @CO line: everything after the first tab is the comment. A @CO without a tab
            // is malformed; reject it (htslib header.c:805-809) rather than dropping it.
            const std::size_t tabPos{line.find('\t')};
            if (tabPos == std::string_view::npos) {
                return std::unexpected{"Missing tab in @CO line"};
            }
            header.comments_.emplace_back(line.substr(tabPos + 1));
            continue;
        }
        // Ignore unknown header line types (forward compatibility)
    }

    return header;
}

std::expected<SamHeader, std::string> SamHeader::FromBamHeaderBlock(std::span<const std::byte> data)
{
    // Minimum: magic(4) + l_text(4) + n_ref(4) = 12 bytes
    if (std::size(data) < 12) {
        return std::unexpected{"BAM header block too short"};
    }

    // Validate magic: BAM\1
    if ((data[0] != std::byte{'B'}) || (data[1] != std::byte{'A'}) || (data[2] != std::byte{'M'}) ||
        (data[3] != std::byte{1})) {
        return std::unexpected{"Invalid BAM magic bytes"};
    }
    const std::uint32_t lText{PacBio::Samoa::ReadU32LE(std::data(data) + 4)};
    if (std::size(data) < 8 + lText + 4) {
        return std::unexpected{"BAM header block truncated in header text"};
    }

    std::string_view headerText{reinterpret_cast<const char*>(std::data(data) + 8), lText};
    while (!std::empty(headerText) && (headerText.back() == '\0')) {
        headerText.remove_suffix(1);
    }
    std::expected<SamHeader, std::string> headerResult{FromText(headerText)};
    if (!headerResult) {
        return std::unexpected{std::move(headerResult.error())};
    }
    SamHeader header{std::move(*headerResult)};

    // Parse binary reference dictionary
    std::size_t offset{8 + lText};
    const std::int32_t nRefSigned{PacBio::Samoa::ReadI32LE(std::data(data) + offset)};
    offset += 4;
    if (nRefSigned < 0) {
        return std::unexpected{"Invalid BAM binary header: negative n_ref"};
    }
    const std::uint32_t nRef{static_cast<std::uint32_t>(nRefSigned)};

    // Build references from binary dict (authoritative for name/length). Bound the
    // reservation by the bytes actually available (each ref entry is >= 8 bytes: 4-byte
    // l_name + name + 4-byte l_ref) so a hostile n_ref cannot trigger std::bad_alloc.
    std::vector<ReferenceSequence> binaryRefs;
    binaryRefs.reserve(
        std::min<std::size_t>(static_cast<std::size_t>(nRef), (std::size(data) - offset) / 8U));

    for (std::uint32_t i = 0; i < nRef; ++i) {
        if (offset + 4 > std::size(data)) {
            return std::unexpected{"BAM header block truncated in reference dictionary"};
        }
        const std::uint32_t lName{PacBio::Samoa::ReadU32LE(std::data(data) + offset)};
        offset += 4;
        if (lName == 0) {
            // l_name includes the NUL terminator, so it must be >= 1 (htslib sam.c:334).
            return std::unexpected{"Invalid BAM binary header: zero-length reference name"};
        }

        if (offset + lName + 4 > std::size(data)) {
            return std::unexpected{"BAM header block truncated in reference name"};
        }

        // Name is NUL-terminated, lName includes the NUL
        const std::size_t nameLength{static_cast<std::size_t>(lName - 1)};
        std::string name{reinterpret_cast<const char*>(std::data(data) + offset), nameLength};
        offset += lName;
        const std::int32_t lRef{PacBio::Samoa::ReadI32LE(std::data(data) + offset)};
        offset += 4;

        binaryRefs.emplace_back(std::move(name), lRef);
    }

    // If @SQ lines were in text, transfer optional tags to binary refs
    if ((!std::empty(header.references_)) &&
        (std::size(header.references_) == std::size(binaryRefs))) {
        for (std::size_t i = 0; i < std::size(binaryRefs); ++i) {
            for (const auto& [key, value] : header.references_[i].CustomTags()) {
                binaryRefs[i].SetTag(key, value);
            }
        }
    }

    // Replace text-parsed refs with authoritative binary dict
    header.references_ = std::move(binaryRefs);
    header.nameIndexSize_ = 0;  // invalidate name index

    return header;
}

std::vector<std::byte> SamHeader::ToBamHeaderBlock() const
{
    std::vector<std::byte> result;

    // magic: BAM\1
    result.push_back(std::byte{'B'});
    result.push_back(std::byte{'A'});
    result.push_back(std::byte{'M'});
    result.push_back(std::byte{1});

    // Header text
    const std::string text{ToText()};
    const std::uint32_t textLen{static_cast<std::uint32_t>(std::size(text))};
    AppendLE(result, textLen);
    for (char c : text) {
        result.push_back(static_cast<std::byte>(c));
    }

    // Reference dictionary
    const std::uint32_t numRefs{static_cast<std::uint32_t>(std::size(references_))};
    AppendLE(result, numRefs);

    for (const ReferenceSequence& ref : references_) {
        // l_name includes NUL terminator
        const std::uint32_t lName{static_cast<std::uint32_t>(std::size(ref.Name()) + 1)};
        AppendLE(result, lName);

        for (char c : ref.Name()) {
            result.push_back(static_cast<std::byte>(c));
        }
        result.push_back(std::byte{0});  // NUL terminator

        AppendLE(result, ref.Length());
    }

    return result;
}

// --- SamHeader serialization ---

std::string SamHeader::ToText() const
{
    std::string result;

    // @HD line: emit if any @HD field is set, not only VN, so a VN-less @HD carrying
    // sort-order metadata round-trips like htslib (build_header_line, header.c:743-756).
    if (!std::empty(version_) || !std::empty(sortOrder_) || !std::empty(groupOrder_) ||
        !std::empty(subSort_) || !std::empty(hdCustomTags_)) {
        result += "@HD";
        if (!std::empty(version_)) {
            AppendSamField(result, "VN", version_);
        }
        if (!std::empty(sortOrder_)) {
            AppendSamField(result, "SO", sortOrder_);
        }
        if (!std::empty(groupOrder_)) {
            AppendSamField(result, "GO", groupOrder_);
        }
        if (!std::empty(subSort_)) {
            AppendSamField(result, "SS", subSort_);
        }
        AppendCustomTags(result, hdCustomTags_);
        result += '\n';
    }

    // @SQ lines
    for (const ReferenceSequence& ref : references_) {
        result += "@SQ";
        AppendSamField(result, "SN", ref.Name());
        AppendSamField(result, "LN", ref.Length());
        AppendCustomTags(result, ref.CustomTags());
        result += '\n';
    }

    // @RG lines
    for (const ReadGroup& rg : readGroups_) {
        result += "@RG";
        AppendSamField(result, "ID", rg.Id());
        AppendCustomTags(result, rg.CustomTags());
        result += '\n';
    }

    // @PG lines
    for (const ProgramRecord& pg : programRecords_) {
        result += "@PG";
        AppendSamField(result, "ID", pg.Id());
        AppendCustomTags(result, pg.CustomTags());
        result += '\n';
    }

    // @CO lines
    for (const std::string_view comment : comments_) {
        result += "@CO\t";
        result += comment;
        result += '\n';
    }

    return result;
}

// --- SamHeader accessors ---

std::string_view SamHeader::Version() const { return version_; }

std::string_view SamHeader::SortOrder() const { return sortOrder_; }

std::string_view SamHeader::GroupOrder() const { return groupOrder_; }

std::string_view SamHeader::SubSort() const { return subSort_; }

void SamHeader::SetVersion(std::string version) { version_ = std::move(version); }

void SamHeader::SetSortOrder(std::string sortOrder) { sortOrder_ = std::move(sortOrder); }

void SamHeader::SetGroupOrder(std::string groupOrder) { groupOrder_ = std::move(groupOrder); }

void SamHeader::SetSubSort(std::string subSort) { subSort_ = std::move(subSort); }

std::span<const ReferenceSequence> SamHeader::ReferenceSequences() const { return references_; }

std::vector<ReferenceSequence>& SamHeader::ReferenceSequences() { return references_; }

void SamHeader::AddReferenceSequence(ReferenceSequence seq)
{
    references_.push_back(std::move(seq));
}

std::span<const ReadGroup> SamHeader::ReadGroups() const { return readGroups_; }

std::vector<ReadGroup>& SamHeader::ReadGroups() { return readGroups_; }

void SamHeader::AddReadGroup(ReadGroup rg) { readGroups_.push_back(std::move(rg)); }

std::span<const ProgramRecord> SamHeader::ProgramRecords() const { return programRecords_; }

std::vector<ProgramRecord>& SamHeader::ProgramRecords() { return programRecords_; }

void SamHeader::AddProgramRecord(ProgramRecord pg) { programRecords_.push_back(std::move(pg)); }

std::span<const std::string> SamHeader::Comments() const { return comments_; }

void SamHeader::AddComment(std::string comment) { comments_.push_back(std::move(comment)); }

// --- Reference name/ID lookup ---

void SamHeader::BuildNameIndex() const
{
    if (nameIndexSize_ == std::size(references_)) {
        return;
    }
    nameToId_.clear();
    nameToId_.reserve(std::size(references_));
    for (std::size_t i{0}; i < std::size(references_); ++i) {
        nameToId_.emplace(std::string{references_[i].Name()}, static_cast<std::int32_t>(i));
    }
    nameIndexSize_ = std::size(references_);
}

std::int32_t SamHeader::ReferenceId(std::string_view name) const
{
    if (name == "*") {
        return -1;
    }
    BuildNameIndex();
    const auto it = nameToId_.find(name);
    if (it == std::ranges::end(nameToId_)) {
        return -1;
    }
    return it->second;
}

std::string_view SamHeader::ReferenceName(std::int32_t refId) const
{
    if (refId == -1) {
        return "*";
    }
    if ((refId < 0) || (refId >= std::ssize(references_))) {
        throw std::out_of_range{std::format("Invalid reference ID: {}", refId)};
    }
    return references_[refId].Name();
}

std::int32_t SamHeader::NumReferences() const { return std::ssize(references_); }

}  // namespace Samoa
}  // namespace PacBio
