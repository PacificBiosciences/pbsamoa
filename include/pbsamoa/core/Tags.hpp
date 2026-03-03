#ifndef PBSAMOA_CORE_TAGS_HPP
#define PBSAMOA_CORE_TAGS_HPP

#include <algorithm>
#include <initializer_list>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include <cstddef>
#include <cstdint>
#include <cstring>

namespace PacBio {
namespace Samoa {

/// \brief Two-character auxiliary tag key, stored as uint16_t.
class TagKey
{
public:
    constexpr TagKey();
    constexpr TagKey(char first, char second);

    constexpr char First() const;
    constexpr char Second() const;
    constexpr std::uint16_t Value() const;
    std::string ToString() const;

    constexpr bool operator==(const TagKey&) const = default;

private:
    std::uint16_t value_;
};

/// \brief Hex byte string (type 'H' in SAM/BAM).
struct HexString
{
    std::string value;
    bool operator==(const HexString&) const = default;
};

/// \brief Typed array for tag type 'B'. Stores raw element bytes + type code.
class TagArray
{
public:
    explicit TagArray(char elementType);

    char ElementType() const;
    std::uint32_t Count() const;
    std::size_t ElementSize() const;
    std::span<const std::byte> Data() const;
    std::span<std::byte> MutableData();

    void AppendInt8(std::int8_t v);
    void AppendUInt8(std::uint8_t v);
    void AppendInt16(std::int16_t v);
    void AppendUInt16(std::uint16_t v);
    void AppendInt32(std::int32_t v);
    void AppendUInt32(std::uint32_t v);
    void AppendFloat(float v);

    /// \brief Resize backing storage for count elements.
    void Resize(std::uint32_t count);

    bool operator==(const TagArray&) const = default;

private:
    char elementType_;
    std::vector<std::byte> data_;
};

/// \brief Auxiliary tag value: variant over all SAM/BAM tag types.
///
/// Integer types (c/C/s/S/i/I in BAM, i in SAM) stored uniformly as int64_t.
/// Serialization to BAM picks the smallest sufficient type.
using TagValue = std::variant<char, std::int64_t, float, std::string, HexString, TagArray>;

/// \brief Ordered collection of tag key-value pairs with linear-scan lookup.
///
/// Records rarely exceed ~30 tags; linear scan beats hash due to cache effects.
class TagMap
{
public:
    using Entry = std::pair<TagKey, TagValue>;

    const TagValue* Get(TagKey key) const;
    void Set(TagKey key, TagValue value);
    bool Remove(TagKey key);
    bool Contains(TagKey key) const;
    /// \brief Append a key-value pair without checking for duplicates. Caller must guarantee the key does not already exist.
    void Append(TagKey key, TagValue value);

    std::span<const Entry> Entries() const;
    std::size_t Size() const;
    bool Empty() const;

private:
    std::vector<Entry> entries_;
};

// --- Tag parsing and serialization (defined in Tags.cpp) ---

/// \brief Parse auxiliary tags from raw BAM bytes.
TagMap ParseTagsFromBam(std::span<const std::byte> data);

/// \brief Parse a single tag from SAM text (e.g., "NM:i:5").
std::optional<std::pair<TagKey, TagValue>> ParseTagFromSam(std::string_view text);

/// \brief Serialize tags to BAM binary encoding.
std::vector<std::byte> SerializeTagsToBam(const TagMap& tags);

/// \brief Serialize a single tag to SAM text (e.g., "NM:i:5").
std::string SerializeTagToSam(TagKey key, const TagValue& value);

/// \brief Serialize raw BAM auxiliary bytes directly to SAM text.
/// Walks binary tag data and appends tab-prefixed SAM tag fields to out.
/// Avoids creating TagMap/TagValue intermediaries for zero-allocation output.
void SerializeRawTagsToSam(std::span<const std::byte> auxData, std::string& out);

// --- Tag filters for ToOwned() ---

/// \brief Sorted set of TagKeys with O(log n) lookup via binary search.
class SortedTagKeySet
{
public:
    SortedTagKeySet(std::initializer_list<TagKey> keys);
    bool Contains(TagKey key) const;

private:
    std::vector<TagKey> keys_;
};

/// \brief Drop listed tags during ToOwned() conversion.
class DropTags
{
public:
    DropTags(std::initializer_list<TagKey> keys);
    bool ShouldDrop(TagKey key) const;

private:
    SortedTagKeySet keys_;
};

/// \brief Keep only listed tags during ToOwned() conversion.
class KeepTags
{
public:
    KeepTags(std::initializer_list<TagKey> keys);
    bool ShouldKeep(TagKey key) const;

private:
    SortedTagKeySet keys_;
};

// --- Zero-allocation tag serialization for BamRecord::SerializeToBam ---

/// \brief Compute the serialized BAM byte size of all tags without allocating.
///
/// \param[in] tags  tag map to measure
/// \returns total byte count that AppendTagsToBam would write
std::size_t SerializedBamSize(const TagMap& tags);

/// \brief Write serialized BAM tag bytes directly into a destination buffer.
///
/// \param[in] tags  tag map to serialize
/// \param[in] dest  destination buffer, must have at least SerializedBamSize(tags) bytes
void AppendTagsToBam(const TagMap& tags, std::byte* dest);

// --- inline constexpr definitions ---

constexpr TagKey::TagKey() : value_{0} {}

constexpr TagKey::TagKey(char first, char second)
    : value_((static_cast<std::uint8_t>(first) << 8) | static_cast<std::uint8_t>(second))
{
}

constexpr char TagKey::First() const { return value_ >> 8; }

constexpr char TagKey::Second() const { return value_ & 0xFF; }

constexpr std::uint16_t TagKey::Value() const { return value_; }

}  // namespace Samoa
}  // namespace PacBio

#endif  // PBSAMOA_CORE_TAGS_HPP
