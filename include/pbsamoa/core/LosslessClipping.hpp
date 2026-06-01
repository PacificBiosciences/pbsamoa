#ifndef PBSAMOA_CORE_LOSSLESSCLIPPING_HPP
#define PBSAMOA_CORE_LOSSLESSCLIPPING_HPP

#include <pbsamoa/core/BamRecord.hpp>

#include <cstdint>

namespace PacBio {
namespace Samoa {

/// \brief Clip a record to the query range [clipLeft, clipRight) while preserving
/// the removed flanks so the operation can be undone.
///
/// Before clipping, the leading [0, clipLeft) and trailing [clipRight, seqLen)
/// portions of the sequence, qualities, and per-base kinetics tags are captured
/// into a msgpack-encoded `ls` tag, chaining any pre-existing `ls` under a
/// `nested` key (so repeated clips remain reversible). Forward tags (pw/ip/fp/fi)
/// are stored in read orientation; reverse tags (rp/ri) are stored mirrored. The
/// record is then clipped with the standard query clipper. The `ls` format matches
/// lima, so `lima-undo` can reconstruct the read.
///
/// Sequence, qualities, per-base kinetics, and base modifications (MM/ML) are all
/// captured for undo. Intended for unmapped reads (lima's use case): the record is
/// clipped to query and restore rebuilds the query fields, but it does not reverse
/// CIGAR/POS, so a mapped record's alignment is not reconstructed on undo.
///
/// No-op if the range is empty, out of bounds, or covers the whole sequence.
/// \throws std::runtime_error if the record carries raw-subread pulse (pc/pt/pq/...)
///         or subread-pileup (sf/sm/sx/sa) per-base tags: the tag clipper trims them
///         but capturing the removed values into `ls` for undo is not yet
///         implemented, so they are rejected rather than silently dropped.
void ClipToQueryLossless(BamRecord& record, std::int32_t clipLeft, std::int32_t clipRight);

/// \brief Inverse of ClipToQueryLossless: reconstruct the pre-clip record from its
/// `ls` tag, unwinding the full `nested` chain, then remove the demultiplexing
/// tags (bc/bq/bx/bl/bt/ql/qt/ls).
/// \returns true if an `ls` tag was present and applied.
[[nodiscard]] bool RestoreFromLossless(BamRecord& record);

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_CORE_LOSSLESSCLIPPING_HPP
