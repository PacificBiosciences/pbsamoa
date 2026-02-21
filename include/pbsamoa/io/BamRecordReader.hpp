#ifndef PBSAMOA_IO_BAMRECORDREADER_HPP
#define PBSAMOA_IO_BAMRECORDREADER_HPP

#include <pbsamoa/core/BamRecord.hpp>
#include <pbsamoa/core/Metrics.hpp>
#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/core/Tags.hpp>
#include <pbsamoa/io/BamRawReader.hpp>

#include <filesystem>
#include <memory>
#include <optional>
#include <variant>

#include <cstddef>

namespace PacBio {
namespace Samoa {

struct BamRecordReaderConfig
{
    BamRawReaderConfig ViewConfig{};
    std::size_t DecodeWorkers{4};
    ByteLimit BatchBudget{ByteLimit{4 * 1024 * 1024}};
    std::size_t OutputCapacity{4096};
    std::variant<std::monostate, DropTags, KeepTags> TagFilter{};
};

class BamRecordReader
{
public:
    explicit BamRecordReader(const std::filesystem::path& path, BamRecordReaderConfig config = {});
    ~BamRecordReader();

    BamRecordReader(const BamRecordReader&) = delete;
    BamRecordReader& operator=(const BamRecordReader&) = delete;
    BamRecordReader(BamRecordReader&&) noexcept;
    BamRecordReader& operator=(BamRecordReader&&) noexcept;

    const SamHeader& Header() const;

    /// \brief Read next pre-decoded record. Returns nullopt at EOF.
    /// \throws std::runtime_error if the producer thread encountered an error.
    std::optional<BamRecord> ReadRecord();

    /// \brief Range interface: for (const auto& rec : reader.Records()) { ... }
    class RecordRange
    {
    public:
        class Iterator
        {
        public:
            Iterator();
            explicit Iterator(BamRecordReader* reader);

            const BamRecord& operator*() const;
            const BamRecord* operator->() const;
            Iterator& operator++();
            void operator++(int);
            bool operator==(const Iterator& other) const;

        private:
            BamRecordReader* reader_{nullptr};
            std::optional<BamRecord> current_;
        };

        explicit RecordRange(BamRecordReader* reader);
        Iterator begin();
        Iterator end();

    private:
        BamRecordReader* reader_;
    };

    RecordRange Records();

    /// \brief Snapshot of all reader metrics (BGZF + decode layers). Thread-safe, lock-free.
    ReaderMetrics GetMetrics() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_IO_BAMRECORDREADER_HPP
