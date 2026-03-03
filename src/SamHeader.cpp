#include <pbsamoa/core/SamHeader.hpp>

#include "BinaryUtils.hpp"

#include <algorithm>
#include <charconv>
#include <expected>
#include <format>
#include <stdexcept>
#include <utility>

namespace {

/// \brief Split string_view by a delimiter character.
/// Returns vector of non-owning views. Skips empty trailing entries.
std::vector<std::string_view> Split(std::string_view sv, char delim)
{
    std::vector<std::string_view> result;
    while (!std::empty(sv)) {
        const std::size_t pos = sv.find(delim);
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
    const std::size_t colon = field.find(':');
    if (colon == std::string_view::npos) {
        return {field, {}};
    }
    return {field.substr(0, colon), field.substr(colon + 1)};
}

}  // namespace

namespace PacBio {
namespace Samoa {
namespace detail {

const std::string* FindTag(std::span<const std::pair<std::string, std::string>> tags,
                           std::string_view key)
{
    const auto it = std::ranges::find(tags, key, &std::pair<std::string, std::string>::first);
    if (it == std::ranges::end(tags)) {
        return nullptr;
    }
    return &it->second;
}

void UpsertTag(std::vector<std::pair<std::string, std::string>>& tags, std::string key,
               std::string value)
{
    auto it = std::ranges::find(tags, key, &std::pair<std::string, std::string>::first);
    if (it \!= std::ranges::end(tags)) {
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
    auto it = std::ranges::find_if(customTags_,
                                   [&key](const std::pair<std::string, std::string>& tagPair) {
                                       return tagPair.first == key;
                                   });
    if (it != std::ranges::end(customTags_)) {
        it->second = std::move(value);
    } else {
        customTags_.emplace_back(std::move(key), std::move(value));
    }
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
    auto it = std::ranges::find_if(customTags_,
                                   [&key](const std::pair<std::string, std::string>& tagPair) {
                                       return tagPair.first == key;
                                   });
    if (it != std::ranges::end(customTags_)) {
        it->second = std::move(value);
    } else {
        customTags_.emplace_back(std::move(key), std::move(value));
    }
}

std::span<const std::pair<std::string, std::string>> ReadGroup::CustomTags() const
{
    return customTags_;
}

// --- ProgramRecord ---

ProgramRecord::ProgramRecord(std::string id) : id_{std::move(id)} {}

std::string_view ProgramRecord::Id() const { return id_; }

const std::string* ProgramRecord::GetTag(std::string_view key) const
{
    return detail::FindTag(customTags_, key);
}

void ProgramRecord::SetTag(std::string key, std::string value)
{
    auto it = std::ranges::find_if(customTags_,
                                   [&key](const std::pair<std::string, std::string>& tagPair) {
                                       return tagPair.first == key;
                                   });
    if (it != std::ranges::end(customTags_)) {
        it->second = std::move(value);
    } else {
        customTags_.emplace_back(std::move(key), std::move(value));
    }
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
    const std::vector<std::string_view> lines = Split(text, '\n');
    for (const std::string_view line : lines) {
        if (std::empty(line)) {
            continue;
        }
        const std::vector<std::string_view> fields = Split(line, '\t');
        if (std::empty(fields)) {
            continue;
        }
        const std::string_view recordType = fields[0];

        if (recordType == "@HD") {
            for (std::size_t i = 1; i < std::size(fields); ++i) {
                const auto [tag, value] = ParseTagValue(fields[i]);
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
            std::string name;
            std::int32_t length{-1};
            std::vector<std::pair<std::string, std::string>> tags;

            for (std::size_t i = 1; i < std::size(fields); ++i) {
                const auto [tag, value] = ParseTagValue(fields[i]);
                if (tag == "SN") {
                    name = std::string{value};
                } else if (tag == "LN") {
                    std::int32_t parsed{0};
                    const std::from_chars_result parseResult{std::from_chars(
                        std::data(value), std::data(value) + std::size(value), parsed)};
                    const std::errc ec = parseResult.ec;
                    if (ec != std::errc{}) {
                        return std::unexpected{"Invalid LN value in @SQ line"};
                    }
                    length = parsed;
                } else {
                    tags.emplace_back(std::string{tag}, std::string{value});
                }
            }

            if (std::empty(name)) {
                return std::unexpected{"Missing required SN field in @SQ line"};
            }
            if (length < 0) {
                return std::unexpected{"Missing required LN field in @SQ line"};
            }

            ReferenceSequence ref{std::move(name), length};
            for (auto& [k, v] : tags) {
                ref.SetTag(std::move(k), std::move(v));
            }
            header.references_.push_back(std::move(ref));

        } else if (recordType == "@RG") {
            std::string id;
            std::vector<std::pair<std::string, std::string>> tags;

            for (std::size_t i = 1; i < std::size(fields); ++i) {
                const auto [tag, value] = ParseTagValue(fields[i]);
                if (tag == "ID") {
                    id = std::string{value};
                } else {
                    tags.emplace_back(std::string{tag}, std::string{value});
                }
            }

            if (std::empty(id)) {
                return std::unexpected{"Missing required ID field in @RG line"};
            }

            ReadGroup rg{std::move(id)};
            for (auto& [k, v] : tags) {
                rg.SetTag(std::move(k), std::move(v));
            }
            header.readGroups_.push_back(std::move(rg));

        } else if (recordType == "@PG") {
            std::string id;
            std::vector<std::pair<std::string, std::string>> tags;

            for (std::size_t i = 1; i < std::size(fields); ++i) {
                const auto [tag, value] = ParseTagValue(fields[i]);
                if (tag == "ID") {
                    id = std::string{value};
                } else {
                    tags.emplace_back(std::string{tag}, std::string{value});
                }
            }

            if (std::empty(id)) {
                return std::unexpected{"Missing required ID field in @PG line"};
            }

            ProgramRecord pg{std::move(id)};
            for (auto& [k, v] : tags) {
                pg.SetTag(std::move(k), std::move(v));
            }
            header.programRecords_.push_back(std::move(pg));

        } else if (recordType == "@CO") {
            // @CO line: everything after first tab is the comment
            const std::size_t tabPos = line.find('\t');
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
    const std::uint32_t lText = PacBio::Samoa::ReadU32LE(std::data(data) + 4);
    if (std::size(data) < 8 + lText + 4) {
        return std::unexpected{"BAM header block truncated in header text"};
    }

    // Parse SAM header text (may be NUL-padded)
    std::string_view headerText{reinterpret_cast<const char*>(std::data(data) + 8), lText};

    // Strip trailing NULs
    while (!std::empty(headerText) && headerText.back() == '\0') {
        headerText.remove_suffix(1);
    }

    auto headerResult{FromText(headerText)};
    if (!headerResult.has_value()) {
        return std::unexpected{std::move(headerResult.error())};
    }
    SamHeader header{std::move(*headerResult)};

    // Parse binary reference dictionary
    std::size_t offset = 8 + lText;
    const std::uint32_t nRef = PacBio::Samoa::ReadU32LE(std::data(data) + offset);
    offset += 4;

    // Build references from binary dict (authoritative for name/length)
    std::vector<ReferenceSequence> binaryRefs;
    binaryRefs.reserve(nRef);

    for (std::uint32_t i = 0; i < nRef; ++i) {
        if (offset + 4 > std::size(data)) {
            return std::unexpected{"BAM header block truncated in reference dictionary"};
        }
        const std::uint32_t lName = PacBio::Samoa::ReadU32LE(std::data(data) + offset);
        offset += 4;

        if (offset + lName + 4 > std::size(data)) {
            return std::unexpected{"BAM header block truncated in reference name"};
        }

        // Name is NUL-terminated, lName includes the NUL
        std::string name{reinterpret_cast<const char*>(std::data(data) + offset),
                         lName > 0 ? lName - 1 : 0};
        offset += lName;
        const std::int32_t lRef = PacBio::Samoa::ReadI32LE(std::data(data) + offset);
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
    const std::string text = ToText();
    const std::uint32_t textLen = std::size(text);
    WriteU32LE(result, textLen);
    for (char c : text) {
        result.push_back(static_cast<std::byte>(c));
    }

    // Reference dictionary
    const std::uint32_t numRefs = std::size(references_);
    WriteU32LE(result, numRefs);

    for (const ReferenceSequence& ref : references_) {
        // l_name includes NUL terminator
        const std::uint32_t lName = std::size(ref.Name()) + 1;
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
        result += "@HD\tVN:";
        result += version_;
        if (!std::empty(sortOrder_)) {
            result += "\tSO:";
            result += sortOrder_;
        }
        if (!std::empty(groupOrder_)) {
            result += "\tGO:";
            result += groupOrder_;
        }
        if (!std::empty(subSort_)) {
            result += "\tSS:";
            result += subSort_;
        }
        result += '\n';
    }

    // @SQ lines
    for (const ReferenceSequence& ref : references_) {
        result += "@SQ\tSN:";
        result += ref.Name();
        result += "\tLN:";
        result += std::format("{}", ref.Length());
        for (const auto& [key, value] : ref.CustomTags()) {
            result += '\t';
            result += key;
            result += ':';
            result += value;
        }
        result += '\n';
    }

    // @RG lines
    for (const ReadGroup& rg : readGroups_) {
        result += "@RG\tID:";
        result += rg.Id();
        for (const auto& [key, value] : rg.CustomTags()) {
            result += '\t';
            result += key;
            result += ':';
            result += value;
        }
        result += '\n';
    }

    // @PG lines
    for (const ProgramRecord& pg : programRecords_) {
        result += "@PG\tID:";
        result += pg.Id();
        for (const auto& [key, value] : pg.CustomTags()) {
            result += '\t';
            result += key;
            result += ':';
            result += value;
        }
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
