#include <pbsamoa/io/BamCollection.hpp>

#include <algorithm>
#include <format>
#include <functional>
#include <ranges>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace PacBio {
namespace Samoa {
namespace {

template <typename Record>
void MergeCustomTags(Record& target, const Record& incoming, std::string_view label)
{
    for (const auto& [key, value] : incoming.CustomTags()) {
        const std::string* existing{target.GetTag(key)};
        if (existing) {
            if (*existing != value) {
                throw std::runtime_error{
                    std::format("BamCollection: conflicting {} tag '{}'", label, key)};
            }
        } else {
            target.SetTag(key, value);
        }
    }
}

std::string MergeHdField(std::string_view lhs, std::string_view rhs, std::string_view label)
{
    if (lhs.empty()) {
        return std::string{rhs};
    }
    if (rhs.empty() || (lhs == rhs)) {
        return std::string{lhs};
    }
    throw std::runtime_error{std::format("BamCollection: conflicting @HD {}", label)};
}

template <typename Record, typename KeyFn>
void MergeNamedRecords(std::vector<Record>& target, const Record& incoming, KeyFn keyFn,
                       std::string_view label)
{
    const std::string_view key{std::invoke(keyFn, incoming)};
    const auto it{std::ranges::find(target, key, keyFn)};
    if (it == std::ranges::end(target)) {
        target.push_back(incoming);
        return;
    }
    MergeCustomTags(*it, incoming, std::format("{} '{}'", label, key));
}

void MergeReferences(SamHeader& merged, const SamHeader& incoming)
{
    auto& mergedRefs{merged.ReferenceSequences()};
    const auto incomingRefs{incoming.ReferenceSequences()};
    const std::size_t numRefs{std::size(mergedRefs)};
    if (numRefs != std::size(incomingRefs)) {
        throw std::runtime_error{"BamCollection: conflicting reference counts"};
    }

    for (std::size_t i{0}; i < numRefs; ++i) {
        if ((mergedRefs[i].Name() != incomingRefs[i].Name()) ||
            (mergedRefs[i].Length() != incomingRefs[i].Length())) {
            throw std::runtime_error{
                std::format("BamCollection: conflicting reference dictionary at index {}", i)};
        }
        MergeCustomTags(mergedRefs[i], incomingRefs[i],
                        std::format("reference '{}'", mergedRefs[i].Name()));
    }
}

void MergeComments(SamHeader& merged, const SamHeader& incoming)
{
    for (const std::string_view comment : incoming.Comments()) {
        // AddComment may reallocate the underlying comment storage, so re-fetch
        // the span on each iteration instead of caching it across insertions.
        if (!std::ranges::contains(merged.Comments(), comment)) {
            merged.AddComment(std::string{comment});
        }
    }
}

void MergeHeader(SamHeader& merged, const SamHeader& incoming)
{
    merged.SetVersion(MergeHdField(merged.Version(), incoming.Version(), "VN"));
    merged.SetSortOrder(MergeHdField(merged.SortOrder(), incoming.SortOrder(), "SO"));
    merged.SetGroupOrder(MergeHdField(merged.GroupOrder(), incoming.GroupOrder(), "GO"));
    merged.SetSubSort(MergeHdField(merged.SubSort(), incoming.SubSort(), "SS"));

    MergeReferences(merged, incoming);
    for (const ReadGroup& readGroup : incoming.ReadGroups()) {
        MergeNamedRecords(merged.ReadGroups(), readGroup, &ReadGroup::Id, "read group");
    }
    for (const ProgramRecord& program : incoming.ProgramRecords()) {
        MergeNamedRecords(merged.ProgramRecords(), program, &ProgramRecord::Id, "program");
    }
    MergeComments(merged, incoming);
}

SamHeader MergeHeaders(std::span<const BamFile> files)
{
    if (std::empty(files)) {
        throw std::invalid_argument{"BamCollection: at least one BAM file is required"};
    }

    SamHeader merged{files.front().Header()};
    for (const BamFile& file : files.subspan(1)) {
        MergeHeader(merged, file.Header());
    }

    if (std::size(files) > 1) {
        // A concatenated multi-file stream should not claim a stronger global
        // ordering guarantee than the collection layer can prove.
        merged.SetSortOrder("unknown");
        merged.SetGroupOrder({});
        merged.SetSubSort({});
    }

    return merged;
}

std::vector<BamFile> ToBamFiles(std::vector<std::filesystem::path> bams)
{
    std::vector<BamFile> bamFiles;
    bamFiles.reserve(std::size(bams));
    for (std::filesystem::path& bam : bams) {
        bamFiles.emplace_back(std::move(bam));
    }
    return bamFiles;
}

}  // namespace

BamCollection::BamCollection(std::filesystem::path bam) : BamCollection{BamFile{std::move(bam)}} {}

BamCollection::BamCollection(BamFile bamFile)
    : BamCollection{std::vector<BamFile>{std::move(bamFile)}}
{
}

BamCollection::BamCollection(std::vector<std::filesystem::path> bams)
    : BamCollection{ToBamFiles(std::move(bams))}
{
}

BamCollection::BamCollection(std::vector<BamFile> bamFiles)
    : files_{std::move(bamFiles)}, header_{MergeHeaders(files_)}
{
}

const SamHeader& BamCollection::Header() const { return header_; }

std::span<const BamFile> BamCollection::Files() const { return files_; }

std::size_t BamCollection::Size() const { return std::size(files_); }

}  // namespace Samoa
}  // namespace PacBio
