#include "Undo.hpp"

#include <pbsamoa/PbSamoaLibraryInfo.hpp>
#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/LosslessClipping.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/core/Tags.hpp>
#include <pbsamoa/io/BamRecordReader.hpp>
#include <pbsamoa/io/BamWriter.hpp>

#include <pbcopper/cli2/Interface.h>
#include <pbcopper/cli2/PositionalArgument.h>
#include <pbcopper/cli2/Results.h>

#include <filesystem>
#include <optional>
#include <print>
#include <ranges>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include <cstdint>
#include <cstdlib>

namespace PacBio {
namespace Samoa {
namespace Undo {
namespace {

const CLI_v2::PositionalArgument Input{
    R"({
    "name" : "input",
    "description" : "Demultiplexed BAM (lima or horizon).",
    "type" : "file"
})"};

const CLI_v2::PositionalArgument Output{
    R"({
    "name" : "output",
    "description" : "Output BAM with barcode clipping undone.",
    "type" : "file"
})"};

// Collapse barcoded read groups ("<base>/<label>") onto their base ID and drop the
// barcode bookkeeping (BC tag, DS Barcode* fields) like lima-undo does.
ReadGroup UnbarcodedReadGroup(const ReadGroup& group)
{
    ReadGroup base{std::string{ReadGroupBaseId(group.Id())}};
    for (const auto& [key, value] : group.CustomTags()) {
        if ((key != "BC") && (key != "DS")) {
            base.SetTag(key, value);
        }
    }
    DsTagFields fields;
    for (const auto& field : group.ParsedDsFields()) {
        if (!field.first.starts_with("Barcode")) {
            fields.push_back(field);
        }
    }
    if (!std::empty(fields)) {
        base.SetDsFields(fields);
    }
    return base;
}

SamHeader UndoneHeader(const SamHeader& input, std::string_view commandLine)
{
    SamHeader header{input};
    std::vector<ReadGroup> groups;
    for (const ReadGroup& group : input.ReadGroups()) {
        const std::string_view baseId{ReadGroupBaseId(group.Id())};
        const bool seen{
            std::ranges::any_of(groups, [baseId](const ReadGroup& g) { return g.Id() == baseId; })};
        if (!seen) {
            groups.push_back(UnbarcodedReadGroup(group));
        }
    }
    header.ReadGroups() = std::move(groups);

    // Drop lima's @PG lines. A surviving @PG whose PP names a dropped record is
    // reparented to that record's own PP, so the chain stays unbroken (SAMv1 §1.3).
    std::vector<ProgramRecord>& programs{header.ProgramRecords()};
    for (auto it{std::begin(programs)}; it != std::end(programs);) {
        if (!it->Id().starts_with("lima")) {
            ++it;
            continue;
        }
        const std::string removedId{it->Id()};
        const std::string* removedPp{it->GetTag("PP")};
        const std::optional<std::string> parent{removedPp ? std::optional{*removedPp}
                                                          : std::nullopt};
        it = programs.erase(it);
        for (ProgramRecord& pg : programs) {
            const std::string* pp{pg.GetTag("PP")};
            if ((pp != nullptr) && (*pp == removedId)) {
                if (parent) {
                    pg.SetTag("PP", *parent);
                } else {
                    pg.RemoveTag("PP");
                }
            }
        }
    }
    AppendProgramRecord(header, "pbsamoa.undo", std::string{commandLine});
    return header;
}

}  // namespace

CLI_v2::Interface CreateInterface()
{
    CLI_v2::Interface interface{"pbsamoa undo",
                                "Undo barcode clipping and demultiplexing tags (lima-compatible).",
                                LibraryFormattedVersion()};
    interface.DisableNumThreadsOption();
    interface.AddPositionalArguments({Input, Output});
    return interface;
}

int Runner(const CLI_v2::Results& results)
{
    const std::vector<std::string>& positional{results.PositionalArguments()};
    if (std::size(positional) != 2) {
        throw std::runtime_error{"undo requires exactly two arguments: <input.bam> <output.bam>"};
    }

    BamRecordReader reader{std::filesystem::path{positional[0]}};
    const std::string commandLine{"pbsamoa " + results.InputCommandLine()};
    BamWriter writer{std::filesystem::path{positional[1]},
                     UndoneHeader(reader.Header(), commandLine)};

    std::int64_t numRecords{0};
    std::int64_t numRestored{0};
    for (BamRecord record : reader.Records()) {
        ++numRecords;
        if (RestoreFromLossless(record)) {
            ++numRestored;
        } else if (const TagValue* rg{record.Tags().Get(RG_TAG)}; rg != nullptr) {
            // No ls stack (unclipped read): still point RG at the collapsed base group.
            if (const std::string* id{std::get_if<std::string>(rg)}; id != nullptr) {
                record.MutableTags().Set(RG_TAG, std::string{ReadGroupBaseId(*id)});
            }
        }
        writer.Write(record);
    }
    writer.Close();
    std::println(stderr, "pbsamoa undo: {} records, {} restored", numRecords, numRestored);
    return EXIT_SUCCESS;
}

}  // namespace Undo
}  // namespace Samoa
}  // namespace PacBio
