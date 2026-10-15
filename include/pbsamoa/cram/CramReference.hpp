#ifndef PBSAMOA_CRAM_CRAMREFERENCE_HPP
#define PBSAMOA_CRAM_CRAMREFERENCE_HPP

#include <string>
#include <string_view>

#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief The reference of a reference-based CRAM file, given in windows: for example from a
/// packed form, so the text of the whole reference need not be in memory.
///
/// CramWriter and CramReader use it. Several threads call Bases at the same time.
class CramReference
{
public:
    virtual ~CramReference() = default;

    /// \brief The number of reference sequences: one per @SQ line, in header order.
    [[nodiscard]] virtual std::int32_t NumSequences() const = 0;

    /// \brief The length of reference sequence refId.
    [[nodiscard]] virtual std::int64_t Length(std::int32_t refId) const = 0;

    /// \brief The bases [start, start + length) of reference sequence refId (inside it).
    ///
    /// The result can be a view of buffer, which the call can change. The case does not
    /// matter.
    [[nodiscard]] virtual std::string_view Bases(std::int32_t refId, std::int64_t start,
                                                 std::int64_t length,
                                                 std::string& buffer) const = 0;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_CRAM_CRAMREFERENCE_HPP
