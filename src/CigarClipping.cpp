#include <pbsamoa/core/CigarClipping.hpp>

#include <pbsamoa/core/CigarOp.hpp>

#include <algorithm>
#include <utility>
#include <vector>

#include <cstddef>
#include <cstdint>

namespace {

PacBio::Samoa::ClipResult MakeQueryClipResult(std::vector<PacBio::Samoa::CigarOp>&& cigar,
                                              std::int32_t queryOffset, std::int32_t queryStart,
                                              std::int32_t queryEnd, std::int32_t newPos)
{
    return PacBio::Samoa::ClipResult{
        std::move(cigar),
        static_cast<std::size_t>(queryOffset),
        static_cast<std::size_t>(queryEnd - queryStart),
        newPos,
    };
}

std::vector<PacBio::Samoa::CigarOp> RemainingOps(std::span<const PacBio::Samoa::CigarOp> cigar,
                                                 std::size_t first)
{
    return {
        std::begin(cigar) + static_cast<std::ptrdiff_t>(first),
        std::end(cigar),
    };
}

void TrimBack(std::vector<PacBio::Samoa::CigarOp>& ops, std::int32_t count)
{
    while (!std::empty(ops) && (count > 0)) {
        const PacBio::Samoa::CigarOp op{ops.back()};
        const PacBio::Samoa::CigarOpType type{op.Type()};
        const std::uint32_t len{op.Length()};

        if (!ConsumesQuery(type)) {
            ops.pop_back();
            continue;
        }

        if (len <= static_cast<std::uint32_t>(count)) {
            count -= static_cast<std::int32_t>(len);
            ops.pop_back();
        } else {
            const std::uint32_t shrunk{len - static_cast<std::uint32_t>(count)};
            ops.back() = PacBio::Samoa::CigarOp{type, shrunk};
            count = 0;
        }
    }
}

}  // namespace

namespace PacBio {
namespace Samoa {

ClipResult ClipCigarToQuery(std::span<const CigarOp> cigar, std::int32_t queryStart,
                            std::int32_t queryEnd, std::int32_t origQueryStart,
                            std::int32_t origQueryEnd, std::int32_t refPos, bool isReverse)
{
    std::int32_t frontRemove{std::max(std::int32_t{0}, queryStart - origQueryStart)};
    std::int32_t backRemove{std::max(std::int32_t{0}, origQueryEnd - queryEnd)};

    if (isReverse) {
        std::swap(frontRemove, backRemove);
    }

    // --- Front trim ---
    std::int32_t refAdvance{0};
    std::int32_t queryOffset{0};
    std::size_t idx{0};

    while ((idx < std::size(cigar)) && (frontRemove > 0)) {
        const CigarOp op{cigar[idx]};
        const CigarOpType type{op.Type()};
        const std::uint32_t len{op.Length()};

        if (!ConsumesQuery(type)) {
            // Non-query-consuming ops (D, N, H, P) are skipped entirely
            // and their reference consumption counts toward refAdvance.
            if (ConsumesReference(type)) {
                refAdvance += static_cast<std::int32_t>(len);
            }
            ++idx;
            continue;
        }

        // Query-consuming op
        if (len <= static_cast<std::uint32_t>(frontRemove)) {
            // Consume entire op
            frontRemove -= static_cast<std::int32_t>(len);
            queryOffset += static_cast<std::int32_t>(len);
            if (ConsumesReference(type)) {
                refAdvance += static_cast<std::int32_t>(len);
            }
            ++idx;
        } else {
            // Partial consumption — remainder becomes first op of result
            const std::uint32_t remainder{len - static_cast<std::uint32_t>(frontRemove)};
            if (ConsumesReference(type)) {
                refAdvance += frontRemove;
            }
            queryOffset += frontRemove;
            frontRemove = 0;

            // Build result starting with the partial op
            const std::vector<CigarOp> remaining{RemainingOps(cigar, idx + 1)};
            std::vector<CigarOp> result;
            result.reserve(1 + std::size(remaining));
            result.emplace_back(type, remainder);
            result.insert(std::end(result), std::begin(remaining), std::end(remaining));

            // --- Back trim ---
            TrimBack(result, backRemove);

            return MakeQueryClipResult(std::move(result), queryOffset, queryStart, queryEnd,
                                       refPos + refAdvance);
        }
    }

    // If we get here, frontRemove was fully consumed by whole ops (or was 0).
    // Copy remaining ops from idx onward.
    std::vector<CigarOp> result{RemainingOps(cigar, idx)};

    // --- Back trim ---
    TrimBack(result, backRemove);

    return MakeQueryClipResult(std::move(result), queryOffset, queryStart, queryEnd,
                               refPos + refAdvance);
}

ClipResult ClipCigarToReference(std::span<const CigarOp> cigar, std::int32_t refStart,
                                std::int32_t refEnd, std::int32_t origRefPos, bool isReverse,
                                bool exciseFlankingInserts)
{
    const std::int64_t refLen{ReferenceLength(cigar)};
    const std::int64_t totalQueryLen{QueryLength(cigar)};

    std::int32_t frontRefRemove{std::max(std::int32_t{0}, refStart - origRefPos)};
    std::int32_t backRefRemove{
        std::max(std::int32_t{0}, static_cast<std::int32_t>(origRefPos + refLen) - refEnd)};

    // Build a mutable working copy of the CIGAR
    std::vector<CigarOp> ops{std::begin(cigar), std::end(cigar)};

    // --- Front trim: consume frontRefRemove reference bases ---
    std::size_t queryRemovedFront{0};
    std::size_t idx{0};

    while ((idx < std::size(ops)) && (frontRefRemove > 0)) {
        const CigarOp op{ops[idx]};
        const CigarOpType type{op.Type()};
        const std::uint32_t len{op.Length()};

        if (!ConsumesReference(type)) {
            // Non-ref ops (S, I, H, P) are popped entirely during the walk
            if (ConsumesQuery(type)) {
                queryRemovedFront += len;
            }
            ++idx;
            continue;
        }

        // Ref-consuming op
        if (len <= static_cast<std::uint32_t>(frontRefRemove)) {
            // Consume entire op
            frontRefRemove -= static_cast<std::int32_t>(len);
            if (ConsumesQuery(type)) {
                queryRemovedFront += len;
            }
            ++idx;
        } else {
            // Partial: shrink this op
            ops[idx] = CigarOp{type, len - static_cast<std::uint32_t>(frontRefRemove)};
            if (ConsumesQuery(type)) {
                queryRemovedFront += static_cast<std::uint32_t>(frontRefRemove);
            }
            frontRefRemove = 0;
        }
    }

    // Remove consumed front ops
    ops.erase(std::begin(ops), std::begin(ops) + static_cast<std::ptrdiff_t>(idx));

    // --- Back trim: consume backRefRemove reference bases from end ---
    std::size_t queryRemovedBack{0};

    while (!std::empty(ops) && (backRefRemove > 0)) {
        const CigarOp op{ops.back()};
        const CigarOpType type{op.Type()};
        const std::uint32_t len{op.Length()};

        if (!ConsumesReference(type)) {
            // Non-ref ops popped entirely
            if (ConsumesQuery(type)) {
                queryRemovedBack += len;
            }
            ops.pop_back();
            continue;
        }

        // Ref-consuming op
        if (len <= static_cast<std::uint32_t>(backRefRemove)) {
            backRefRemove -= static_cast<std::int32_t>(len);
            if (ConsumesQuery(type)) {
                queryRemovedBack += len;
            }
            ops.pop_back();
        } else {
            ops.back() = CigarOp{type, len - static_cast<std::uint32_t>(backRefRemove)};
            if (ConsumesQuery(type)) {
                queryRemovedBack += static_cast<std::uint32_t>(backRefRemove);
            }
            backRefRemove = 0;
        }
    }

    // --- Excise flanking inserts ---
    if (exciseFlankingInserts) {
        if (!std::empty(ops) && (ops.front().Type() == CigarOpType::I)) {
            queryRemovedFront += ops.front().Length();
            ops.erase(std::begin(ops));
        }
        if (!std::empty(ops) && (ops.back().Type() == CigarOpType::I)) {
            queryRemovedBack += ops.back().Length();
            ops.pop_back();
        }
    }

    // --- Reverse strand: swap front/back query removal ---
    if (isReverse) {
        std::swap(queryRemovedFront, queryRemovedBack);
    }

    const std::size_t clipOffset{queryRemovedFront};
    const std::size_t totalQuery{static_cast<std::size_t>(totalQueryLen)};
    const std::size_t totalRemoved{queryRemovedFront + queryRemovedBack};
    const std::size_t clipLength{(totalRemoved <= totalQuery) ? (totalQuery - totalRemoved) : 0U};
    const std::int32_t newPos{std::max(origRefPos, refStart)};

    return ClipResult{
        std::move(ops),
        clipOffset,
        clipLength,
        newPos,
    };
}

}  // namespace Samoa
}  // namespace PacBio
