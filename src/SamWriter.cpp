#include <pbsamoa/io/SamWriter.hpp>

#include "WriterUtils.hpp"

#include <algorithm>
#include <format>
#include <fstream>
#include <string>

#include <cstdint>

namespace PacBio {
namespace Samoa {
namespace {

constexpr std::uint8_t MISSING_BAM_QUALITY{0xFF};

std::int32_t ToSamPosition(const std::int32_t zeroBasedPos)
{
    if (zeroBasedPos < 0) {
        return 0;
    }
    return zeroBasedPos + 1;
}

bool QualitiesUnavailable(std::span<const std::uint8_t> qual)
{
    return std::empty(qual) || std::ranges::all_of(qual, [](const std::uint8_t quality) {
               return quality == MISSING_BAM_QUALITY;
           });
}

void FormatFields(std::string& buf, const SamHeader& header, std::string_view name,
                  std::uint16_t flag, std::int32_t refId, std::int32_t pos, std::uint8_t mapq,
                  CigarView cigar, std::int32_t nextRefId, std::int32_t nextPos, std::int32_t tlen,
                  std::string_view seq, std::span<const std::uint8_t> qual)
{
    buf.clear();

    // QNAME
    buf.append(name);
    buf += '\t';

    // FLAG
    std::format_to(std::back_inserter(buf), "{}", flag);
    buf += '\t';

    // RNAME
    if (refId < 0) {
        buf += '*';
    } else {
        buf.append(header.ReferenceName(refId));
    }
    buf += '\t';

    // POS (0-based -> 1-based)
    std::format_to(std::back_inserter(buf), "{}", ToSamPosition(pos));
    buf += '\t';

    // MAPQ
    std::format_to(std::back_inserter(buf), "{}", static_cast<unsigned>(mapq));
    buf += '\t';

    // CIGAR
    if (std::empty(cigar)) {
        buf += '*';
    } else {
        buf.append(CigarToString(cigar));
    }
    buf += '\t';

    // RNEXT
    if (nextRefId < 0) {
        buf += '*';
    } else if (nextRefId == refId) {
        buf += '=';
    } else {
        buf.append(header.ReferenceName(nextRefId));
    }
    buf += '\t';

    // PNEXT (0-based -> 1-based)
    std::format_to(std::back_inserter(buf), "{}", ToSamPosition(nextPos));
    buf += '\t';

    // TLEN
    std::format_to(std::back_inserter(buf), "{}", tlen);
    buf += '\t';

    // SEQ
    if (std::empty(seq)) {
        buf += '*';
    } else {
        buf.append(seq);
    }
    buf += '\t';

    // QUAL - 0xFF means unavailable in BAM; output '*'
    if (QualitiesUnavailable(qual)) {
        buf += '*';
    } else {
        for (const std::uint8_t quality : qual) {
            buf += static_cast<char>(quality + 33);
        }
    }
}

void FlushLine(std::string& buf, std::ofstream& stream)
{
    buf += '\n';
    stream.write(std::data(buf), static_cast<std::streamsize>(std::size(buf)));
}

}  // namespace

struct SamWriter::Impl
{
    std::filesystem::path path;
    std::filesystem::path writePath;
    SamWriterConfig config;
    std::ofstream stream;
    SamHeader header;
    std::string buf;  // reusable line buffer
    bool closed{false};

    Impl(const std::filesystem::path& outputPath, const SamHeader& hdr, const SamWriterConfig& cfg)
        : path{outputPath}, config{cfg}, header{hdr}
    {
        writePath = detail::ResolveWritePath(path, config.UseTempFile);
        stream.open(writePath, std::ios::out | std::ios::trunc);
        if (!stream) {
            throw std::runtime_error{std::format("SamWriter: cannot open {}", writePath.string())};
        }
        stream << header.ToText();
    }
};

SamWriter::SamWriter(const std::filesystem::path& path, const SamHeader& header,
                     const SamWriterConfig& config)
    : impl_{std::make_unique<Impl>(path, header, config)}
{
}

SamWriter::~SamWriter()
{
    if ((impl_ != nullptr) && (!impl_->closed)) {
        Close();
    }
}

SamWriter::SamWriter(SamWriter&&) noexcept = default;
SamWriter& SamWriter::operator=(SamWriter&&) noexcept = default;

void SamWriter::Write(const BamRecord& record)
{
    FormatFields(impl_->buf, impl_->header, record.Name(), record.Flag(), record.RefId(),
                 record.Pos(), record.MapQ(), record.Cigar(), record.NextRefId(), record.NextPos(),
                 record.Tlen(), record.Sequence(), record.Qualities());

    for (const auto& [key, value] : record.Tags().Entries()) {
        impl_->buf += '\t';
        impl_->buf.append(SerializeTagToSam(key, value));
    }

    FlushLine(impl_->buf, impl_->stream);
}

void SamWriter::Write(const RawRecord& view)
{
    FormatFields(impl_->buf, impl_->header, view.Name(), view.Flag(), view.RefId(), view.Pos(),
                 view.MapQ(), view.CigarOps(), view.NextRefId(), view.NextPos(), view.Tlen(),
                 view.Seq().ToString(), view.Qual());

    SerializeRawTagsToSam(view.AuxData(), impl_->buf);

    FlushLine(impl_->buf, impl_->stream);
}

void SamWriter::WriteBatch(const RawRecordBatch& batch)
{
    const std::size_t recordCount{batch.RecordCount()};
    for (std::size_t i{0}; i < recordCount; ++i) {
        const RawRecord view{batch.RecordData(i)};
        Write(view);
    }
}

void SamWriter::Close()
{
    if (impl_->closed) {
        return;
    }
    impl_->stream.flush();
    if (!impl_->stream) {
        throw std::runtime_error{"SamWriter: flush failed"};
    }
    impl_->stream.close();
    if (!impl_->stream) {
        throw std::runtime_error{
            std::format("SamWriter: failed to close output {}", impl_->writePath.string())};
    }

    if (impl_->config.UseTempFile) {
        detail::AtomicRename(impl_->writePath, impl_->path, "SamWriter");
    }

    impl_->closed = true;
}

}  // namespace Samoa
}  // namespace PacBio
