# API Reference

All types live in `namespace PacBio::Samoa`. Headers are under
`include/pbsamoa/`, organized into `core/`, `io/`, and `index/`
subdirectories.

---

## BamRawReader

```cpp
#include <pbsamoa/io/BamRawReader.hpp>
```

Reads BAM files. Supports sequential iteration, batch reading, and BAI
region queries. Returns non-owning views into raw BAM bytes.

### Construction

```cpp
// Basic — single-threaded decompression
BamRawReader reader{"input.bam"};

// With parallel BGZF decompression (4 worker threads)
BamRawReader reader{"input.bam", BamRawReaderConfig{.BgzfWorkers = 4}};

// Stop after 1000 records
BamRawReader reader{"input.bam", BamRawReaderConfig{.RecordLimit = 1000}};

// Read chunk 2 of 4 (requires .zmi or .pbi alongside the BAM; ChunkNum is 1-based)
BamRawReader reader{"input.bam", BamRawReaderConfig{.ChunkNum = 2, .TotalChunks = 4}};
```

`BamRawReaderConfig` fields:

|     Field     | Default |                 Description                  |
| ------------- | ------- | -------------------------------------------- |
| `BgzfWorkers` | `0`     | BGZF decompression threads (0 = synchronous) |
| `RecordLimit` | `0`     | Stop after N records (0 = unlimited)         |
| `ChunkNum`    | `0`     | 1-based chunk number (0 = no chunking)       |
| `TotalChunks` | `0`     | Total chunk count (0 = no chunking)          |

### Header access

```cpp
const SamHeader& header = reader.Header();
```

### Range iteration

```cpp
#include <print>

for (const auto& view : reader.Records()) {
    std::println("{}", view.Name());
}
```

Returns `RawRecord` references. Each view is valid until the next
iteration step (views share an internal buffer).

### Batch reading

```cpp
using namespace PacBio::Samoa::Literals;

while (auto batch = reader.ReadBatch(128_MiB)) {
    for (std::size_t i{0}; i < batch->RecordCount(); ++i) {
        const RawRecord view{batch->RecordData(i)};
        // process
    }
}
```

Returns `std::optional<RawRecordBatch>`. `nullopt` at EOF.

### Region query

```cpp
auto index = BaiIndex::FromFile("sorted.bam.bai");
std::int32_t refId = reader.Header().ReferenceId("chr1");

for (const auto& view : reader.Query(index, refId, 1000, 2000)) {
    // records overlapping chr1:1000-2000
}
```

### Random access

```cpp
VirtualOffset pos = reader.Tell();
// ... read some records ...
reader.Seek(pos);  // jump back
```

---

## BamRecordReader

```cpp
#include <pbsamoa/io/BamRecordReader.hpp>
```

Wraps `BamRawReader` and delivers pre-decoded owned `BamRecord` objects.
A background producer thread reads batches, decodes views to records
(optionally in parallel via a thread pool), and pushes them through an
SPSC queue for consumption.

### Construction

```cpp
// Default config — 4 decode workers, 4 MiB batch budget
BamRecordReader reader{"input.bam"};

// Serial decoding, no parallelism
BamRecordReader reader{"input.bam", BamRecordReaderConfig{.DecodeWorkers = 0}};

// Parallel BGZF decompression + parallel ToOwned decode
BamRecordReader reader{"input.bam", BamRecordReaderConfig{
    .ViewConfig = {.BgzfWorkers = 4},
    .DecodeWorkers = 8,
}};

// With tag filtering — drop large PacBio kinetics tags during decode
BamRecordReader reader{"input.bam", BamRecordReaderConfig{
    .TagFilter = DropTags{TagKey{'i','p'}, TagKey{'p','w'}},
}};
```

### Configuration

|      Field       |                   Type                   | Default |                     Description                      |
| ---------------- | ---------------------------------------- | ------- | ---------------------------------------------------- |
| `ViewConfig`     | `BamRawReaderConfig`                     | `{}`    | Config for underlying `BamRawReader`                 |
| `DecodeWorkers`  | `size_t`                                 | `4`     | Threads for parallel `ToOwned()` decode (0 = serial) |
| `BatchBudget`    | `ByteLimit`                              | `4_MiB` | Memory budget per batch read from view reader        |
| `OutputCapacity` | `size_t`                                 | `4096`  | SPSC queue capacity (records)                        |
| `TagFilter`      | `variant<monostate, DropTags, KeepTags>` | `{}`    | Tag filter applied during `ToOwned()`                |

### Reading records

```cpp
while (auto record = reader.ReadRecord()) {
    // record is an owned BamRecord
}
```

### Range iteration

```cpp
#include <print>

for (const auto& record : reader.Records()) {
    std::println("{}", record.Name());
}
```

---

## SamReader

```cpp
#include <pbsamoa/io/SamReader.hpp>
```

Reads text SAM files, producing owned `BamRecord` objects.

```cpp
#include <print>
#include <pbsamoa/io/SamReader.hpp>

SamReader reader{"input.sam"};

for (const auto& record : reader.Records()) {
    std::println("{}", record.Name());
}
```

---

## BamWriter

```cpp
#include <pbsamoa/io/BamWriter.hpp>
```

Writes BGZF-compressed BAM files.

### Construction

```cpp
// Default compression (level 6)
BamWriter writer{"output.bam", header};

// Custom compression level (1 = fastest, 12 = smallest)
BamWriter writer{"output.bam", header, 1};

// With index callback (for building BAI alongside writing)
BamWriter writer{"output.bam", header, 6, [](std::int64_t vOffset,
                                              std::span<const std::byte> raw) {
    // track offsets for index building
}};
```

### Writing records

```cpp
writer.Write(view);    // from RawRecord
writer.Write(record);  // serializes from BamRecord
writer.Write(rawData); // raw span<const byte>
writer.WriteBatch(batch);  // batch of records
writer.Close();        // flush and finalize (also called by destructor)
```

---

## SamWriter

```cpp
#include <pbsamoa/io/SamWriter.hpp>
```

Writes text SAM files.

```cpp
SamWriter writer{"output.sam", header};
writer.Write(view);    // accepts RawRecord
writer.Write(record);  // accepts BamRecord
writer.Close();
```

---

## RawRecord

```cpp
#include <pbsamoa/core/RawRecord.hpp>
```

Owning type that copies raw BAM bytes. Decodes on demand — all decode
logic is inline in the header.

### Fixed fields

|      Method      | Return type |                Description                |
| ---------------- | ----------- | ----------------------------------------- |
| `RefId()`        | `int32_t`   | Reference sequence ID (-1 = unmapped)     |
| `Pos()`          | `int32_t`   | 0-based leftmost position (-1 = unmapped) |
| `MapQ()`         | `uint8_t`   | Mapping quality                           |
| `Flag()`         | `uint16_t`  | Bitwise flag                              |
| `NextRefId()`    | `int32_t`   | Mate reference ID                         |
| `NextPos()`      | `int32_t`   | Mate 0-based position                     |
| `Tlen()`         | `int32_t`   | Template length                           |
| `Bin()`          | `uint16_t`  | BAI bin                                   |
| `NameLength()`   | `uint8_t`   | Length of read name (incl. NUL)           |
| `CigarOpCount()` | `uint16_t`  | Number of CIGAR operations                |
| `SeqLength()`    | `uint32_t`  | Sequence length                           |

### Variable-length fields

|    Method    |      Return type      |               Description                |
| ------------ | --------------------- | ---------------------------------------- |
| `Name()`     | `string_view`         | Read name (NUL-terminated in buffer)     |
| `CigarOps()` | `CigarView`           | Span of CIGAR operations                 |
| `Seq()`      | `SequenceView`        | 4-bit packed sequence (decode on access) |
| `Qual()`     | `span<const uint8_t>` | Base qualities                           |
| `AuxData()`  | `span<const byte>`    | Raw auxiliary tag bytes                  |

### Derived fields

|       Method        | Return type |              Description              |
| ------------------- | ----------- | ------------------------------------- |
| `IsMapped()`        | `bool`      | `(Flag() & 4) == 0`                   |
| `IsReverseStrand()` | `bool`      | `(Flag() & 16) != 0`                  |
| `IsPrimary()`       | `bool`      | Not secondary and not supplementary   |
| `ReferenceLength()` | `int64_t`   | Consumed reference bases (from CIGAR) |
| `QueryLength()`     | `int64_t`   | Consumed query bases (from CIGAR)     |

### Conversion

```cpp
BamRecord record = view.ToOwned();
BamRecord record = view.ToOwned(DropTags{TagKey{'i','p'}, TagKey{'p','w'}});
BamRecord record = view.ToOwned(KeepTags{TagKey{'R','G'}});
TagMap tags = view.ParseTags();
```

---

## BamRecord

```cpp
#include <pbsamoa/core/BamRecord.hpp>
```

Owned, mutable record. Fields stored as C++ types. Mutations are plain
assignments; serialization happens once in the writer.

### Accessors

`Name()`, `Flag()`, `RefId()`, `Pos()`, `MapQ()`, `NextRefId()`,
`NextPos()`, `Tlen()`, `Tags()` have the same names as `RawRecord`. Three
accessors differ:

|  `BamRecord`  | `RawRecord`  |
| ------------- | ------------ |
| `Cigar()`     | `CigarOps()` |
| `Sequence()`  | `Seq()`      |
| `Qualities()` | `Qual()`     |

### Fluent mutators

Every accessor has a corresponding mutator returning `BamRecord&`:

```cpp
record.Name("read1")
      .Flag(99)
      .RefId(0)
      .Pos(1000)
      .MapQ(60)
      .Cigar({CigarOp{CigarOpType::M, 100}})
      .Sequence("ACGT...")
      .Qualities({30, 30, 30, 30});
```

### Serialization

```cpp
std::vector<std::byte> raw = record.SerializeToBam();
```

Normally you don't call this directly — `BamWriter::Write(record)` handles
it.

---

## RawRecordBatch

```cpp
#include <pbsamoa/core/RawRecord.hpp>
```

Owns a decompressed buffer and provides views into it.

|     Method      |    Return type     |            Description            |
| --------------- | ------------------ | --------------------------------- |
| `RecordData(i)` | `span<const byte>` | Raw data for record `i`           |
| `RecordCount()` | `size_t`           | Number of records                 |
| `BufferSize()`  | `size_t`           | Decompressed buffer size in bytes |

Views are invalidated when the batch is destroyed.

---

## SamHeader

```cpp
#include <pbsamoa/core/SamHeader.hpp>
```

Parsed SAM/BAM header.

### Construction

```cpp
auto header = SamHeader::FromText(samText);
auto header = SamHeader::FromBamHeaderBlock(rawBytes);
```

### Serialization

```cpp
std::string text = header.ToText();
std::vector<std::byte> bytes = header.ToBamHeaderBlock();
```

### Header fields (`@HD`)

|     Method     |                   Description                    |
| -------------- | ------------------------------------------------ |
| `Version()`    | Format version (e.g. "1.6")                      |
| `SortOrder()`  | "coordinate", "queryname", "unsorted", "unknown" |
| `GroupOrder()` | Grouping of alignments                           |
| `SubSort()`    | Sub-sorting order within groups                  |

### Reference sequences (`@SQ`)

```cpp
header.AddReferenceSequence(ReferenceSequence{"chr1", 248956422});

std::int32_t id = header.ReferenceId("chr1");       // name → ID
std::string_view name = header.ReferenceName(id);    // ID → name
std::int32_t n = header.NumReferences();
```

### Read groups, programs, comments

```cpp
header.AddReadGroup(ReadGroup{"rg1"});
header.AddProgramRecord(ProgramRecord{"pbsamoa"});
header.AddComment("processed by pbsamoa");
```

---

## BaiIndex

```cpp
#include <pbsamoa/index/BaiIndex.hpp>
```

BAM index: read, write, build, and query.

### Load and query

```cpp
auto index = BaiIndex::FromFile("sorted.bam.bai");
auto chunks = index.Query(refId, 1000, 2000);
```

### Build from BAM

```cpp
auto index = BaiIndex::Build("sorted.bam");
index.ToFile("sorted.bam.bai");
```

### Statistics

|      Method       |                Description                 |
| ----------------- | ------------------------------------------ |
| `NumReferences()` | Number of reference sequences in the index |
| `MappedCount()`   | Total mapped reads (from metadata bin)     |
| `UnmappedCount()` | Total unmapped reads (from metadata bin)   |

---

## CigarOp

```cpp
#include <pbsamoa/core/CigarOp.hpp>
```

32-bit CIGAR operation matching BAM binary layout.

```cpp
CigarOp op{CigarOpType::M, 100};  // 100M
```

### Operations

| Enum value | SAM char |             Description             | Consumes query | Consumes ref |
| ---------- | -------- | ----------------------------------- | :------------: | :----------: |
| `M`        | `M`      | Alignment match (match or mismatch) |      yes       |     yes      |
| `I`        | `I`      | Insertion to the reference          |      yes       |              |
| `D`        | `D`      | Deletion from the reference         |                |     yes      |
| `N`        | `N`      | Skipped region from the reference   |                |     yes      |
| `S`        | `S`      | Soft clipping                       |      yes       |              |
| `H`        | `H`      | Hard clipping                       |                |              |
| `P`        | `P`      | Padding                             |                |              |
| `EQ`       | `=`      | Sequence match                      |      yes       |     yes      |
| `X`        | `X`      | Sequence mismatch                   |      yes       |     yes      |

### Functions

```cpp
std::int64_t refLen = ReferenceLength(cigarOps);
std::int64_t qLen = QueryLength(cigarOps);
std::string str = CigarToString(cigarOps);         // "100M50I25M"
std::vector<CigarOp> ops = ParseCigar("100M50I25M");
std::uint16_t bin = Reg2Bin(beg, end);
```

---

## Tags

```cpp
#include <pbsamoa/core/Tags.hpp>
```

### TagKey

Two-character key stored as `uint16_t`:

```cpp
constexpr TagKey nm{'N', 'M'};
```

### TagValue

Variant: `char | int64_t | float | string | HexString | TagArray`.

Integer types (`c/C/s/S/i/I` in BAM) are stored uniformly as `int64_t`.
Serialization picks the smallest sufficient type.

### TagMap

```cpp
TagMap tags;
tags.Set(TagKey{'N','M'}, std::int64_t{5});

if (const auto* val = tags.Get(TagKey{'N','M'})) {
    // val is a TagValue*
}

for (const auto& [key, value] : tags.Entries()) {
    // iterate all tags
}
```

### Parsing and serialization

```cpp
TagMap tags = ParseTagsFromBam(auxBytes);
auto kv = ParseTagFromSam("NM:i:5");
std::vector<std::byte> raw = SerializeTagsToBam(tags);
std::string text = SerializeTagToSam(TagKey{'N','M'}, TagValue{std::int64_t{5}});
```

### Tag filters

```cpp
DropTags drop{TagKey{'i','p'}, TagKey{'p','w'}};
KeepTags keep{TagKey{'R','G'}, TagKey{'N','M'}};

auto record = view.ToOwned(drop);  // copy everything except ip, pw
auto record = view.ToOwned(keep);  // copy only RG, NM
```

---

## VirtualOffset

```cpp
#include <pbsamoa/core/Bgzf.hpp>
```

Type-safe 64-bit BGZF virtual offset.

```cpp
VirtualOffset vo{blockOffset, withinBlockOffset};
std::uint64_t block = vo.BlockOffset();
std::uint16_t within = vo.WithinBlockOffset();
```

Supports `<=>` comparison. No arithmetic — this is by design.

---

## BgzfReader

```cpp
#include <pbsamoa/core/Bgzf.hpp>
```

Block-by-block BGZF decompression.

```cpp
BgzfReader reader{"file.bam"};
std::vector<std::byte> buf(65536);

while (auto n = reader.ReadBlock(buf)) {
    // *n bytes of decompressed data in buf
}
```

---

## BgzfWriter

```cpp
#include <pbsamoa/core/Bgzf.hpp>
```

BGZF compression. Accumulates into 64 KiB blocks, compresses with
libdeflate.

```cpp
BgzfWriter writer{"output.gz", 6};  // compression level 1-12
writer.Write(data);
writer.Close();  // appends EOF marker
```

---

## ByteLimit

```cpp
#include <pbsamoa/io/BamRawReader.hpp>
```

Memory budget for batch reading. Default 256 MiB.

```cpp
using namespace PacBio::Samoa::Literals;

ByteLimit limit = 128_MiB;
ByteLimit small = 512_KiB;
ByteLimit precise{1024 * 1024 * 64};  // 64 MiB
```

---

## Error Policies

```cpp
#include <pbsamoa/io/BamRawReader.hpp>
```

|    Policy     |                        Behavior                         |
| ------------- | ------------------------------------------------------- |
| `ThrowPolicy` | Throws `std::runtime_error` on corrupt records          |
| `SkipPolicy`  | Logs and skips. `SkippedCount()` returns total skipped. |

---

## ZmwIndex

```cpp
#include <pbsamoa/index/ZmwIndex.hpp>
```

PacBio ZMW index — maps ZMW hole numbers to virtual offsets.

```cpp
auto index = ZmwIndex::Open("movie.bam");  // auto-detects .zmi or .pbi
std::vector<std::int64_t> offsets = index.Find(12345);  // raw virtual offsets by hole number
```

## ZmwGroup

```cpp
#include <pbsamoa/io/BamZmwReader.hpp>
```

A single ZMW's worth of owned records.

```cpp
struct ZmwGroup {
    ZmwIdentity zmw;
    std::vector<BamRecord> records;
};
```

---

## BamZmwReader

```cpp
#include <pbsamoa/io/BamZmwReader.hpp>
```

Groups consecutive records by ZMW identity. Wraps a `BamRecordReader` and
does synchronous grouping — no background threads. Assumes the source
yields records grouped by ZMW (standard PacBio BAM ordering).

### Construction

```cpp
// Takes ownership of BamRecordReader
BamZmwReader reader{BamRecordReader{"movie.bam"}};

// With custom config
BamZmwReader reader{BamRecordReader{"movie.bam",
    BamRecordReaderConfig{.DecodeWorkers = 0}}};
```

### Reading

```cpp
std::vector<BamRecord> records;
while (reader.GetNext(records)) {
    ZmwIdentity zmw = reader.CurrentZmw();
    // records contains all alignments for one ZMW
}
```

## ZmiBamWriter

```cpp
#include <pbsamoa/io/ZmiBamWriter.hpp>
```

Writes BAM and `.zmi` index simultaneously.

```cpp
ZmiBamWriter writer{"output.bam", header};
writer.Write(record);  // updates both BAM and ZMI
writer.Close();
```
