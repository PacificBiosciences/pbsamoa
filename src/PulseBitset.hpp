#ifndef PBSAMOA_PULSEBITSET_HPP
#define PBSAMOA_PULSEBITSET_HPP

#include <string_view>
#include <vector>

#include <bit>
#include <cctype>
#include <cstddef>
#include <cstdint>

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
    static constexpr std::size_t NPOS{static_cast<std::size_t>(-1)};

    /// \brief Construct from a pulse-call string.
    /// \param[in] pulseCalls  The `pc` tag value where uppercase = basecalled.
    inline explicit PulseBitset(std::string_view pulseCalls)
        : size_{std::size(pulseCalls)}, count_{0}
    {
        const std::size_t numBlocks{(size_ + 63) / 64};
        blocks_.resize(numBlocks, 0);

        for (std::size_t i{0}; i < size_; ++i) {
            if (std::isupper(static_cast<unsigned char>(pulseCalls[i]))) {
                const std::size_t blockIdx{i / 64};
                const std::size_t bitIdx{i % 64};
                blocks_[blockIdx] |= (std::uint64_t{1} << bitIdx);
            }
        }

        for (std::size_t b{0}; b < std::size(blocks_); ++b) {
            count_ += static_cast<std::size_t>(std::popcount(blocks_[b]));
        }
    }

    /// \brief Total number of pulses.
    inline std::size_t Size() const { return size_; }

    /// \brief Number of basecalled pulses (set bits).
    inline std::size_t Count() const { return count_; }

    /// \brief Test whether the pulse at pos is basecalled.
    inline bool IsSet(std::size_t pos) const
    {
        if (pos >= size_) {
            return false;
        }
        const std::size_t blockIdx{pos / 64};
        const std::size_t bitIdx{pos % 64};
        return (blocks_[blockIdx] & (std::uint64_t{1} << bitIdx)) != 0;
    }

    /// \brief Find the index of the first basecalled pulse.
    /// \returns Pulse index, or NPOS if none.
    inline std::size_t FindFirst() const
    {
        for (std::size_t b{0}; b < std::size(blocks_); ++b) {
            if (blocks_[b] != 0) {
                const std::size_t pos{b * 64 +
                                      static_cast<std::size_t>(std::countr_zero(blocks_[b]))};
                return (pos < size_) ? pos : NPOS;
            }
        }
        return NPOS;
    }

    /// \brief Find the next basecalled pulse after pos.
    /// \returns Pulse index, or NPOS if none.
    inline std::size_t FindNext(std::size_t pos) const
    {
        if (pos + 1 >= size_) {
            return NPOS;
        }
        const std::size_t start{pos + 1};
        std::size_t blockIdx{start / 64};
        const std::size_t bitIdx{start % 64};

        // Mask off bits below start within the first block
        const std::uint64_t masked{blocks_[blockIdx] >> bitIdx};
        if (masked != 0) {
            const std::size_t result{blockIdx * 64 + bitIdx +
                                     static_cast<std::size_t>(std::countr_zero(masked))};
            return (result < size_) ? result : NPOS;
        }

        // Scan remaining blocks
        for (++blockIdx; blockIdx < std::size(blocks_); ++blockIdx) {
            if (blocks_[blockIdx] != 0) {
                const std::size_t result{
                    blockIdx * 64 + static_cast<std::size_t>(std::countr_zero(blocks_[blockIdx]))};
                return (result < size_) ? result : NPOS;
            }
        }
        return NPOS;
    }

    /// \brief Find the pulse index of the Nth basecalled pulse (0-indexed).
    /// \returns Pulse index, or NPOS if fewer than n+1 basecalled pulses.
    inline std::size_t FindNthBase(std::size_t n) const
    {
        if (n >= count_) {
            return NPOS;
        }

        std::size_t remaining{n};
        for (std::size_t b{0}; b < std::size(blocks_); ++b) {
            const std::size_t pop{static_cast<std::size_t>(std::popcount(blocks_[b]))};
            if (remaining < pop) {
                // The target bit is in this block
                std::uint64_t word{blocks_[b]};
                for (std::size_t i{0}; i < remaining; ++i) {
                    word &= word - 1;  // clear lowest set bit
                }
                return b * 64 + static_cast<std::size_t>(std::countr_zero(word));
            }
            remaining -= pop;
        }
        return NPOS;  // unreachable if count_ is correct
    }

private:
    std::vector<std::uint64_t> blocks_;
    std::size_t size_;
    std::size_t count_;
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_PULSEBITSET_HPP
