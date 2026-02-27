#ifndef PBSAMOA_IO_SAMREADER_HPP
#define PBSAMOA_IO_SAMREADER_HPP

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/SamHeader.hpp>

#include <filesystem>
#include <memory>
#include <optional>
#include <ranges>

namespace PacBio {
namespace Samoa {

class SamReader
{
public:
    explicit SamReader(const std::filesystem::path& path);
    ~SamReader();

    SamReader(const SamReader&) = delete;
    SamReader& operator=(const SamReader&) = delete;
    SamReader(SamReader&&) noexcept;
    SamReader& operator=(SamReader&&) noexcept;

    const SamHeader& Header() const;

    std::optional<BamRecord> ReadRecord();

    class RecordRange : public std::ranges::view_interface<RecordRange>
    {
    public:
        class Iterator
        {
        public:
            Iterator();
            explicit Iterator(SamReader* reader);

            const BamRecord& operator*() const;
            const BamRecord* operator->() const;
            Iterator& operator++();
            void operator++(int);
            bool operator==(const Iterator& other) const;

        private:
            SamReader* reader_{nullptr};
            std::optional<BamRecord> current_;
        };

        explicit RecordRange(SamReader* reader);
        Iterator begin();
        Iterator end();

    private:
        SamReader* reader_;
    };

    RecordRange Records();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;

    BamRecord ParseAlignmentLine(std::string_view line);
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_IO_SAMREADER_HPP
