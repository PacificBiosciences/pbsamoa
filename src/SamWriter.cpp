#include <pbsamoa/io/SamWriter.hpp>

#include <algorithm>
#include <format>
#include <fstream>
#include <string>

#include <cstdint>

namespace PacBio {
namespace Samoa {
namespace {

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
    std::format_to(std::back_inserter(buf), "{}", pos >= 0 ? pos + 1 : 0);
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
    std::format_to(std::back_inserter(buf), "{}", nextPos >= 0 ? nextPos + 1 : 0);
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
    const bool qualUnavailable{std::empty(qual) ||
                               std::ranges::all_of(qual, [](std::uint8_t q) { return q == 0xFF; })};
    if (qualUnavailable) {
        buf += '*';
    } else {
        for (const std::uint8_t q : qual) {
            buf += static_cast<char>(q + 33);
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
    std::ofstream stream;
    SamHeader header;
    std::string buf;  // reusable line buffer
    bool closed{false};

    Impl(const std::filesystem::path& path, const SamHeader& hdr) : stream{path}, header{hdr}
    {
        if (!stream) {
            throw std::runtime_error{std::format("SamWriter: cannot open {}", path.string())};
        }
        stream << header.ToText();
    }
};

SamWriter::SamWriter(const std::filesystem::path& path, const SamHeader& header)
    : impl_{std::make_unique<Impl>(path, header)}
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
    for (std::size_t i{0}; i < batch.RecordCount(); ++i) {
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
    impl_->stream.close();
    impl_->closed = true;
}

}  // namespace Samoa
}  // namespace PacBio
