#include <pbsamoa/core/SamHeader.hpp>

#include "BinaryUtils.hpp"
#include "Md5.hpp"

#include <algorithm>
#include <charconv>
#include <expected>
#include <format>
#include <stdexcept>
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

std::expected<std::int32_t, std::string> ParseInt32Field(std::string_view value,
                                                         std::string_view errorMessage)
{
    std::int32_t parsed{0};
    const std::from_chars_result result{
        std::from_chars(std::data(value), std::data(value) + std::size(value), parsed)};
    if (result.ec != std::errc{}) {
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
            if (!parsedLength.has_value()) {
                return std::unexpected{std::move(parsedLength.error())};
            }
            length = *parsedLength;
        } else {
            tags.emplace_back(std::string{tag}, std::string{value});
        }
    }

    if (std::empty(name)) {
        return std::unexpected{"Missing required SN field in @SQ line"};
    }
    if (!length.has_value() || (*length < 0)) {
        return std::unexpected{"Missing required LN field in @SQ line"};
    }

    PacBio::Samoa::ReferenceSequence ref{std::move(name), *length};
    MoveCustomTagsToRecord(tags, ref);
    return ref;
}

}  // namespace

namespace PacBio {
namespace Samoa {
namespace detail {

const std::string* FindTag(std::span<const std::pair<std::string, std::string>> tags,
                           std::string_view key)
{
    const auto it{std::ranges::find(tags, key, &std::pair<std::string, std::string>::first)};
    if (it == std::ranges::end(tags)) {
        return nullptr;
    }
    return &it->second;
}

void UpsertTag(std::vector<std::pair<std::string, std::string>>& tags, std::string key,
               std::string value)
{
    auto it{std::ranges::find(tags, key, &std::pair<std::string, std::string>::first)};
    if (it != std::ranges::end(tags)) {
        it->second = std::move(value);
    } else {
        tags.emplace_back(std::move(key), std::move(value));
    }
}

}  // namespace detail

// --- ReferenceSequence ---

ReferenceSequence::ReferenceSequence(std::string name, std::int32_t length)
    : name_{std::move(name)}, length_{length}
{
}

std::string_view ReferenceSequence::Name() const { return name_; }

std::int32_t ReferenceSequence::Length() const { return length_; }

const std::string* ReferenceSequence::GetTag(std::string_view key) const
{
    return detail::FindTag(customTags_, key);
}

void ReferenceSequence::SetTag(std::string key, std::string value)
{
    detail::UpsertTag(customTags_, std::move(key), std::move(value));
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
    return detail::FindTag(customTags_, key);
}

void ReadGroup::SetTag(std::string key, std::string value)
{
    if (key == "DS") {
        InvalidateDsCache();
    }
    detail::UpsertTag(customTags_, std::move(key), std::move(value));
}

std::span<const std::pair<std::string, std::string>> ReadGroup::CustomTags() const
{
    return customTags_;
}

// --- PacBio convenience accessors ---

const std::string* ReadGroup::MovieName() const { return GetTag("PU"); }

// --- PacBio DS tag field accessors ---

namespace {

std::string TrimCopy(std::string_view text)
{
    const auto first{static_cast<std::size_t>(
        std::ranges::find_if(text, [](unsigned char ch) { return std::isspace(ch) == 0; }) -
        text.begin())};
    std::size_t last{std::size(text)};
    while ((last > first) && (std::isspace(static_cast<unsigned char>(text[last - 1])) != 0)) {
        --last;
    }
    return std::string{text.substr(first, last - first)};
}

DsTagFields ParseDsTagImpl(std::string_view dsTag)
{
    DsTagFields fields;

    std::size_t start{0};
    while (start <= std::size(dsTag)) {
        const std::size_t end{dsTag.find(';', start)};
        const std::string_view token{(end == std::string_view::npos)
                                         ? dsTag.substr(start)
                                         : dsTag.substr(start, end - start)};
        if (!token.empty()) {
            const std::size_t equals{token.find('=')};
            const std::string key{TrimCopy(token.substr(0, equals))};
            const std::string value{(equals == std::string_view::npos)
                                        ? std::string{}
                                        : TrimCopy(token.substr(equals + 1))};
            fields.emplace_back(key, value);
        }

        if (end == std::string_view::npos) {
            break;
        }
        start = end + 1;
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

const std::string* FindInDsFields(const DsTagFields& fields, std::string_view key)
{
    const auto it{std::ranges::find(fields, key, &std::pair<std::string, std::string>::first)};
    return (it == fields.end()) ? nullptr : &it->second;
}

void UpsertDsField(DsTagFields& fields, std::string key, std::string value)
{
    auto it{std::ranges::find(fields, key, &std::pair<std::string, std::string>::first)};
    if (it != fields.end()) {
        it->second = std::move(value);
    } else {
        fields.emplace_back(std::move(key), std::move(value));
    }
}

}  // namespace

const DsTagFields& ReadGroup::EnsureParsedDs() const
{
    if (!parsedDs_.has_value()) {
        const std::string* dsTag{GetTag("DS")};
        parsedDs_ = dsTag ? ParseDsTagImpl(*dsTag) : DsTagFields{};
    }
    return *parsedDs_;
}

void ReadGroup::InvalidateDsCache() { parsedDs_.reset(); }

const std::string* ReadGroup::DsField(std::string_view key) const
{
    return FindInDsFields(EnsureParsedDs(), key);
}

void ReadGroup::SetDsField(std::string key, std::string value)
{
    auto fields{EnsureParsedDs()};  // copy
    UpsertDsField(fields, std::move(key), std::move(value));
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
    std::string content{std::string{movieName} + "//" + std::string{readType}};
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

// ReadGroupBaseId is defined in ZmwUtils.cpp (same namespace, same signature).

// --- ProgramRecord ---

ProgramRecord::ProgramRecord(std::string id) : id_{std::move(id)} {}

std::string_view ProgramRecord::Id() const { return id_; }

const std::string* ProgramRecord::GetTag(std::string_view key) const
{
    return detail::FindTag(customTags_, key);
}

void ProgramRecord::SetTag(std::string key, std::string value)
{
    detail::UpsertTag(customTags_, std::move(key), std::move(value));
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
            for (const std::string_view field :
                 std::span<const std::string_view>{fields}.subspan(1)) {
                const auto [tag, value] = ParseTagValue(field);
                if (tag == "VN") {
                    header.version_ = std::string{value};
                } else if (tag == "SO") {
                    header.sortOrder_ = std::string{value};
                } else if (tag == "GO") {
                    header.groupOrder_ = std::string{value};
                } else if (tag == "SS") {
                    header.subSort_ = std::string{value};
                }
            }
        } else if (recordType == "@SQ") {
            auto reference{ParseReferenceSequenceRecord(fields)};
            if (!reference.has_value()) {
                return std::unexpected{std::move(reference.error())};
            }
            header.references_.push_back(std::move(*reference));

        } else if (recordType == "@RG") {
            auto readGroup{ParseTaggedRecord<ReadGroup>(fields, "ID",
                                                        "Missing required ID field in @RG line")};
            if (!readGroup.has_value()) {
                return std::unexpected{std::move(readGroup.error())};
            }
            header.readGroups_.push_back(std::move(*readGroup));

        } else if (recordType == "@PG") {
            auto programRecord{ParseTaggedRecord<ProgramRecord>(
                fields, "ID", "Missing required ID field in @PG line")};
            if (!programRecord.has_value()) {
                return std::unexpected{std::move(programRecord.error())};
            }
            header.programRecords_.push_back(std::move(*programRecord));

        } else if (recordType == "@CO") {
            // @CO line: everything after first tab is the comment
            const std::size_t tabPos{line.find('\t')};
            if (tabPos != std::string_view::npos) {
                header.comments_.emplace_back(line.substr(tabPos + 1));
            }
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

    // Parse SAM header text (may be NUL-padded)
    std::string_view headerText{reinterpret_cast<const char*>(std::data(data) + 8), lText};

    // Strip trailing NULs
    while (!std::empty(headerText) && headerText.back() == '\0') {
        headerText.remove_suffix(1);
    }

    std::expected<SamHeader, std::string> headerResult{FromText(headerText)};
    if (!headerResult.has_value()) {
        return std::unexpected{std::move(headerResult.error())};
    }
    SamHeader header{std::move(*headerResult)};

    // Parse binary reference dictionary
    std::size_t offset{8 + lText};
    const std::uint32_t nRef{PacBio::Samoa::ReadU32LE(std::data(data) + offset)};
    offset += 4;

    // Build references from binary dict (authoritative for name/length)
    std::vector<ReferenceSequence> binaryRefs;
    binaryRefs.reserve(nRef);

    for (std::uint32_t i = 0; i < nRef; ++i) {
        if (offset + 4 > std::size(data)) {
            return std::unexpected{"BAM header block truncated in reference dictionary"};
        }
        const std::uint32_t lName{PacBio::Samoa::ReadU32LE(std::data(data) + offset)};
        offset += 4;

        if (offset + lName + 4 > std::size(data)) {
            return std::unexpected{"BAM header block truncated in reference name"};
        }

        // Name is NUL-terminated, lName includes the NUL
        std::string name{reinterpret_cast<const char*>(std::data(data) + offset),
                         lName > 0 ? lName - 1 : 0};
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
    WriteU32LE(result, textLen);
    for (char c : text) {
        result.push_back(static_cast<std::byte>(c));
    }

    // Reference dictionary
    const std::uint32_t numRefs{static_cast<std::uint32_t>(std::size(references_))};
    WriteU32LE(result, numRefs);

    for (const ReferenceSequence& ref : references_) {
        // l_name includes NUL terminator
        const std::uint32_t lName{static_cast<std::uint32_t>(std::size(ref.Name()) + 1)};
        WriteU32LE(result, lName);

        for (char c : ref.Name()) {
            result.push_back(static_cast<std::byte>(c));
        }
        result.push_back(std::byte{0});  // NUL terminator

        WriteI32LE(result, ref.Length());
    }

    return result;
}

// --- SamHeader serialization ---

std::string SamHeader::ToText() const
{
    std::string result;

    // @HD line (only if version is set)
    if (!std::empty(version_)) {
        result += "@HD";
        AppendSamField(result, "VN", version_);
        if (!std::empty(sortOrder_)) {
            AppendSamField(result, "SO", sortOrder_);
        }
        if (!std::empty(groupOrder_)) {
            AppendSamField(result, "GO", groupOrder_);
        }
        if (!std::empty(subSort_)) {
            AppendSamField(result, "SS", subSort_);
        }
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
