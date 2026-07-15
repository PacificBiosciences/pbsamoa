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
/// Provides hardware-accelerated FindNthBase via
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
        blocks_.resize((size_ + BITS_PER_BLOCK - 1) / BITS_PER_BLOCK, 0);

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
    std::size_t size_{0};
    std::size_t count_{0};
};

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_PULSEBITSET_HPP
