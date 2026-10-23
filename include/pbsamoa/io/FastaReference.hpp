#ifndef PBSAMOA_IO_FASTAREFERENCE_HPP
#define PBSAMOA_IO_FASTAREFERENCE_HPP

#include <pbsamoa/core/SamHeader.hpp>
#include <pbsamoa/cram/CramReference.hpp>

#include <filesystem>
#include <memory>

#include <cstddef>

namespace PacBio {
namespace Samoa {

/// \brief What FastaFile::Reference does with an @SQ line whose sequence is not in the FASTA.
enum class MissingFastaSequence
{
    EMPTY,  ///< the sequence has length 0
    THROW,  ///< std::runtime_error
};

/// \brief A memory-mapped FASTA file, read in windows: the text of the whole genome is never
/// copied.
///
/// The sequence offsets come from `<fasta>.fai` when it exists, else from one scan of the file.
/// A sequence with lines of uneven length, which an index cannot describe, is read into memory
/// instead. A gzip-compressed file has no sequences. Copies share the mapping.
class FastaFile
{
public:
    /// \throws std::runtime_error if the file cannot be opened or mapped, or the .fai does not
    ///         fit the file.
    explicit FastaFile(const std::filesystem::path& path);

    [[nodiscard]] std::size_t NumSequences() const;

    /// \brief The reference of the @SQ lines of header, in header order, by name. It keeps the
    /// file mapped.
    ///
    /// \throws std::runtime_error with MissingFastaSequence::THROW if an @SQ sequence is not
    ///         in the file.
    [[nodiscard]] std::shared_ptr<const CramReference> Reference(
        const SamHeader& header, MissingFastaSequence missing = MissingFastaSequence::EMPTY) const;

    struct Data;

private:
    std::shared_ptr<const Data> data_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_IO_FASTAREFERENCE_HPP
