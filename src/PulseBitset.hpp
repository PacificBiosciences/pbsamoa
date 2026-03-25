#ifndef PBSAMOA_PULSEBITSET_HPP
#define PBSAMOA_PULSEBITSET_HPP

#include <string_view>
#include <vector>

#include <bit>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace PacBio {
namespace Samoa {

/// \brief Bitset built from a PacBio pulse-call string (the `pc` tag).
///
/// Uppercase characters represent basecalled pulses (bit set),
/// lowercase characters represent squashed pulses (bit clear).
/// Provides hardware-accelerated FindNext/FindNthBase via
/// std::countr_zero and std::popcount.
class PulseBitset
{
public:
    static constexpr std::size_t NPOS{std::numeric_limits<std::size_t>::max()};
    static constexpr std::size_t BITS_PER_BLOCK{64};

    /// \brief Construct from a pulse-call string.
    /// \param[in] pulseCalls  The `pc` tag value where uppercase = basecalled.
    explicit PulseBitset(std::string_view pulseCalls) : size_{std::size(pulseCalls)}
    {
        const std::size_t numBlocks{(size_ + (BITS_PER_BLOCK - 1)) / BITS_PER_BLOCK};
        blocks_.resize(numBlocks, 0);

        for (std::size_t i{0}; i < size_; ++i) {
            if (std::isupper(static_cast<unsigned char>(pulseCalls[i]))) {
                const std::size_t blockIdx{i / BITS_PER_BLOCK};
                const std::size_t bitIdx{i % BITS_PER_BLOCK};
                blocks_[blockIdx] |= (std::uint64_t{1} << bitIdx);
                ++count_;
            }
        }
    }

    /// \brief Total number of pulses.
    std::size_t Size() const { return size_; }

    /// \brief Number of basecalled pulses (set bits).
    std::size_t Count() const { return count_; }

    /// \brief Test whether the pulse at pos is basecalled.
    bool IsSet(std::size_t pos) const
    {
        if (pos >= size_) {
            return false;
        }
        const std::size_t blockIdx{pos / BITS_PER_BLOCK};
        const std::size_t bitIdx{pos % BITS_PER_BLOCK};
        return (blocks_[blockIdx] & (std::uint64_t{1} << bitIdx)) != 0;
    }

    /// \brief Find the index of the first basecalled pulse.
    /// \returns Pulse index, or NPOS if none.
    std::size_t FindFirst() const
    {
        const std::size_t blockCount{std::size(blocks_)};
        for (std::size_t blockIndex{0}; blockIndex < blockCount; ++blockIndex) {
            if (blocks_[blockIndex] != 0) {
                return blockIndex * BITS_PER_BLOCK +
                       static_cast<std::size_t>(std::countr_zero(blocks_[blockIndex]));
            }
        }
        return NPOS;
    }

    /// \brief Find the next basecalled pulse after pos.
    /// \returns Pulse index, or NPOS if none.
    std::size_t FindNext(std::size_t pos) const
    {
        const std::size_t nextPos{pos + 1};
        if (nextPos >= size_) {
            return NPOS;
        }
        const std::size_t blockCount{std::size(blocks_)};
        std::size_t blockIdx{nextPos / BITS_PER_BLOCK};
        const std::size_t bitIdx{nextPos % BITS_PER_BLOCK};

        // Mask off bits below start within the first block
        const std::uint64_t masked{blocks_[blockIdx] >> bitIdx};
        if (masked != 0) {
            return blockIdx * BITS_PER_BLOCK + bitIdx +
                   static_cast<std::size_t>(std::countr_zero(masked));
        }

        // Scan remaining blocks
        for (++blockIdx; blockIdx < blockCount; ++blockIdx) {
            if (blocks_[blockIdx] != 0) {
                return blockIdx * BITS_PER_BLOCK +
                       static_cast<std::size_t>(std::countr_zero(blocks_[blockIdx]));
            }
        }
        return NPOS;
    }

    /// \brief Find the pulse index of the Nth basecalled pulse (0-indexed).
    /// \returns Pulse index, or NPOS if fewer than n+1 basecalled pulses.
    std::size_t FindNthBase(std::size_t n) const
    {
        if (n >= count_) {
            return NPOS;
        }

        std::size_t remaining{n};
        const std::size_t blockCount{std::size(blocks_)};
        for (std::size_t blockIndex{0}; blockIndex < blockCount; ++blockIndex) {
            const std::size_t pop{static_cast<std::size_t>(std::popcount(blocks_[blockIndex]))};
            if (remaining < pop) {
                // The target bit is in this block
                std::uint64_t word{blocks_[blockIndex]};
                for (std::size_t i{0}; i < remaining; ++i) {
                    word &= word - 1;  // clear lowest set bit
                }
                return blockIndex * BITS_PER_BLOCK +
                       static_cast<std::size_t>(std::countr_zero(word));
            }
            remaining -= pop;
        }
        return NPOS;
    }

private:
    std::vector<std::uint64_t> blocks_;
    std::size_t size_;
    std::size_t count_{0};
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_PULSEBITSET_HPP
