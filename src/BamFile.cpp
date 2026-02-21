#include <pbsamoa/io/BamFile.hpp>

#include <pbsamoa/core/Bgzf.hpp>
#include <pbsamoa/index/BaiIndex.hpp>
#include <pbsamoa/io/BamRawReader.hpp>

#include <filesystem>
#include <stdexcept>
#include <string>
#include <utility>

#include <cstdint>

namespace PacBio {
namespace Samoa {
namespace {

std::filesystem::path IndexPathFor(const std::filesystem::path& bamPath)
{
    std::filesystem::path indexPath{bamPath};
    indexPath += ".bai";
    return indexPath;
}

bool StandardIndexNeedsRefresh(const std::filesystem::path& bamPath)
{
    const std::filesystem::path indexPath{IndexPathFor(bamPath)};
    if (!std::filesystem::exists(indexPath) || !std::filesystem::exists(bamPath)) {
        return true;
    }
    return std::filesystem::last_write_time(indexPath) < std::filesystem::last_write_time(bamPath);
}

}  // namespace

BamFile::BamFile(std::filesystem::path fileName) : fileName_{std::move(fileName)}
{
    BamRawReader reader{fileName_};
    header_ = reader.Header();
}

const std::filesystem::path& BamFile::Filename() const { return fileName_; }

void BamFile::CreateStandardIndex() const
{
    const BaiIndex index{BaiIndex::Build(fileName_)};
    index.ToFile(StandardIndexFilename());
}

void BamFile::EnsureStandardIndexExists() const
{
    if (StandardIndexNeedsRefresh(fileName_)) {
        CreateStandardIndex();
    }
}

bool BamFile::HasEOF() const
{
    const BgzfReader reader{fileName_};
    return reader.HasEofMarker();
}

bool BamFile::StandardIndexExists() const
{
    return std::filesystem::exists(IndexPathFor(fileName_));
}

std::filesystem::path BamFile::StandardIndexFilename() const { return IndexPathFor(fileName_); }

bool BamFile::StandardIndexIsNewer() const { return !StandardIndexNeedsRefresh(fileName_); }

bool BamFile::HasReference(std::string_view name) const { return ReferenceId(name) >= 0; }

const SamHeader& BamFile::Header() const { return header_; }

std::int32_t BamFile::ReferenceId(std::string_view name) const { return header_.ReferenceId(name); }

std::string BamFile::ReferenceName(std::int32_t id) const
{
    try {
        return std::string{header_.ReferenceName(id)};
    } catch (const std::out_of_range&) {
        return {};
    }
}

std::uint32_t BamFile::ReferenceLength(std::string_view name) const
{
    const std::int32_t id{ReferenceId(name)};
    return ReferenceLength(id);
}

std::uint32_t BamFile::ReferenceLength(std::int32_t id) const
{
    const auto refs{header_.ReferenceSequences()};
    if ((id < 0) || (id >= std::ssize(refs))) {
        return 0;
    }
    return static_cast<std::uint32_t>(refs[id].Length());
}

}  // namespace Samoa
}  // namespace PacBio
