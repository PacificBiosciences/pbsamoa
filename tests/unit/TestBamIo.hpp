#ifndef PBSAMOA_TESTS_UNIT_TESTBAMIO_HPP
#define PBSAMOA_TESTS_UNIT_TESTBAMIO_HPP

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/io/BamRawReader.hpp>
#include <pbsamoa/io/BamWriter.hpp>

#include <filesystem>
#include <string>
#include <vector>

namespace PacBio {
namespace Samoa {
namespace tests {

inline void WriteBam(const std::filesystem::path& path, const SamHeader& header,
                     const std::vector<BamRecord>& records)
{
    BamWriter writer{path, header};
    for (const BamRecord& record : records) {
        writer.Write(record);
    }
    writer.Close();
}

inline std::vector<std::string> ReadNames(const std::filesystem::path& path)
{
    BamRawReader reader{path};
    std::vector<std::string> names;
    for (const auto& view : reader.Records()) {
        names.emplace_back(view.Name());
    }
    return names;
}

}  // namespace tests
}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_TESTS_UNIT_TESTBAMIO_HPP
