# API Reference

All types live in `namespace PacBio::Samoa`. Headers are under
`include/pbsamoa/`, organized into `core/`, `io/`, `index/`, and `cram/`
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

|     Field     | Default |                  Description                  |
| ------------- | ------- | --------------------------------------------- |
| `BgzfWorkers` | `0`     | BGZF decompression threads (0 = synchronous)  |
| `RecordLimit` | `0`     | Stop after N records (0 = unlimited)          |
| `ChunkNum`    | `0`     | 1-based chunk number (0 = no chunking)        |
| `TotalChunks` | `0`     | Total chunk count (0 = no chunking)           |
| `Whitelist`   | `{}`    | Optional `ZmwWhitelist` for whitelist queries |

### Header access

```cpp
const SamHeader& header = reader.Header();
```

### Single-record reading

```cpp
while (auto view = reader.ReadRecord()) {
    // view is an owning RawRecord, safe to hold across calls
}
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

### ZMW whitelist query

```cpp
#include <pbsamoa/index/ZmwWhitelist.hpp>

ZmwWhitelist wl{std::vector<std::int32_t>{12345, 23456}};
for (const auto& view : reader.Whitelist(wl)) {
    // records for ZMWs 12345 and 23456 via seek-per-record
}
```

The `WhitelistRange` owns a dedicated sync reader internally (pipeline mode
would restart per seek).

### Random access

```cpp
VirtualOffset pos = reader.Tell();
// ... read some records ...
reader.Seek(pos);  // jump back
```

### Metrics

```cpp
BgzfMetrics m = reader.GetMetrics();
// m.BytesRead, m.BlocksRead, m.BytesDecompressed, m.RecordsConsumed, ...
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
    .RawReaderConfig = {.BgzfWorkers = 4},
    .DecodeWorkers = 8,
}};

// With tag filtering — drop large PacBio kinetics tags during decode
BamRecordReader reader{"input.bam", BamRecordReaderConfig{
    .TagFilter = DropTags{TagKey{'i','p'}, TagKey{'p','w'}},
}};
```

### Configuration

|       Field       |                   Type                   | Default |                     Description                      |
| ----------------- | ---------------------------------------- | ------- | ---------------------------------------------------- |
| `RawReaderConfig` | `BamRawReaderConfig`                     | `{}`    | Config for underlying `BamRawReader`                 |
| `DecodeWorkers`   | `size_t`                                 | `4`     | Threads for parallel `ToOwned()` decode (0 = serial) |
| `BatchBudget`     | `ByteLimit`                              | `4_MiB` | Memory budget per batch read from view reader        |
| `OutputCapacity`  | `size_t`                                 | `4096`  | SPSC queue capacity (records)                        |
| `TagFilter`       | `variant<monostate, DropTags, KeepTags>` | `{}`    | Tag filter applied during `ToOwned()`                |

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

### Metrics

```cpp
ReaderMetrics m = reader.GetMetrics();
// m.Bgzf.BytesRead, m.Decode.RecordsDecoded, m.TotalRecordsRead, ...
```

---

## CramReader

```cpp
#include <pbsamoa/io/CramReader.hpp>
```

Reads CRAM v3.x files and yields owned `BamRecord` objects, with optional
`RawRecord` accessors for BAM-layout output.

### Construction

```cpp
// Default config
CramReader reader{"input.cram"};

// With explicit reference FASTA and parallel decompression
CramReader reader{"input.cram", CramReaderConfig{
    .ReferencePath = "ref.fa",
    .DecompressionWorkers = 4,
}};
```

`CramReaderConfig` fields:

| Field | Default | Description |
| ----- | ------- | ----------- |
| `ReferencePath` | `{}` | Optional FASTA path for reference-based slices |
| `DecompressionWorkers` | `0` | CRAM block decompression workers for slice data blocks (`0` = synchronous) |

### Header access

```cpp
const SamHeader& header = reader.Header();
```

### Single-record reading

```cpp
while (auto record = reader.ReadRecord()) {
    // record is an owned BamRecord
}
```

### Single-record reading as RawRecord

```cpp
while (auto raw = reader.ReadRawRecord()) {
    // raw is an owning RawRecord (BAM binary layout)
}
```

### Range iteration

```cpp
for (const auto& record : reader.Records()) {
    // process record
}
```

### Range iteration as RawRecord

```cpp
for (const auto& raw : reader.RawRecords()) {
    // process raw BAM-layout records
}
```

### Region query with CRAI

```cpp
#include <pbsamoa/index/CraiIndex.hpp>

const CraiIndex index = CraiIndex::FromFile("input.cram.crai");
const std::int32_t refId = reader.Header().ReferenceId("chr1");
const auto records = reader.Query(index, refId, 1000, 2000);  // [beg,end), 0-based
const auto rawRecords = reader.QueryRaw(index, refId, 1000, 2000);
```

Use `refId = -1` for unmapped query rows.

---

## CramWriter

```cpp
#include <pbsamoa/io/CramWriter.hpp>
```

Writes CRAM v3.0/v3.1 files from `BamRecord` and `RawRecord` input.

### Construction

```cpp
// Default config
CramWriter writer{"output.cram", header};

// Tuned config with CRAI sidecar
CramWriter writer{"output.cram", header, CramWriterConfig{
    .BlockCompressionMethod = CramBlockMethod::RANS4X8,
    .RecordsPerSlice = 10000,
    .SlicesPerContainer = 8,
    .CompressionLevel = 1,  // gzip levels: [0, 12]
    .UseTempFile = true,
    .WriteCrai = true,
    .CompressionWorkers = 4,
}};
```

`CramWriterConfig` fields:

| Field | Type | Default | Description |
| ----- | ---- | ------- | ----------- |
| `BlockCompressionMethod` | `CramBlockMethod` | `GZIP` | Default block method |
| `DataSeriesCompressionMethods` | `unordered_map<CramDataSeries, CramBlockMethod>` | `{}` | Per-data-series method overrides |
| `RecordsPerSlice` | `int32_t` | `10000` | Max records per slice |
| `SlicesPerContainer` | `int32_t` | `1` | Max slices per container (`<= 0` means unlimited until `Close()`) |
| `MajorVersion` | `uint8_t` | `3` | CRAM major version (`3` only) |
| `MinorVersion` | `uint8_t` | `0` | CRAM minor version (`0` or `1`) |
| `WriteCrai` | `bool` | `false` | Write `<output>.crai` sidecar |
| `CraiPath` | `optional<path>` | `{}` | CRAI output path override |
| `CompressionWorkers` | `size_t` | `0` | CRAM block compression threads (`0` = synchronous) |
| `CompressionLevel` | `optional<int>` | `{}` | Optional gzip compression level (`[0,12]`) |
| `UseTempFile` | `bool` | `false` | Write to a temporary path and rename atomically on `Close()` |

Notes:
- `MinorVersion` is automatically raised to `1` when selected codecs require CRAM 3.1.
- `CompressionLevel` is currently applied to gzip-compressed blocks; non-gzip codecs ignore it.

### Containerization and CRAI behavior

- Records are buffered into slices of up to `RecordsPerSlice`.
- Slices are buffered into containers of up to `SlicesPerContainer`.
- With `SlicesPerContainer <= 0`, containers are flushed only at `Close()`.
- With `WriteCrai = true`, CRAI rows are emitted per serialized slice offset.
- Mixed-reference slices (`RefSeqId == -2`) emit one CRAI row per referenced sequence,
  plus one unmapped row if unmapped records are present.

### Writing records

```cpp
std::vector<BamRecord> records = /* ... */;
for (const auto& record : records) {
    writer.Write(record);
}
writer.Write(rawRecord);
writer.WriteBatch(rawBatch);
writer.Close();
```

### Low-level CRAM container/slice APIs

```cpp
#include <pbsamoa/cram/CramStructs.hpp>
```

Use these when working directly with serialized CRAM container payloads:

- `CramSlice` and `CramContainer` aggregate parsed wire-format pieces.
- `SerializeSlice(slice)` and `SerializedSliceSize(slice)` for slice-level output/size accounting.
- `SerializeContainer(container)` for complete container bytes (header + payload).
- `ParseContainer(header, payload)` to parse a full payload into slices.

`ParseContainer` parses and decompresses compression/slice headers, while
core/external slice data blocks remain in their on-disk compressed form.

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
// Default writer config
BamWriter defaultWriter{"output.bam", header};

// Custom BGZF settings
BamWriterConfig cfg{
    .BgzfConfig = {
        .CompressionLevel = 1,      // 1 = fastest, 12 = smallest
        .BgzfWorkers = 4,
        .InputQueueCapacity = 256,
        .BlocksPerBatch = 32,
    },
    .UseTempFile = true,  // write via temp + atomic rename on Close()
};
BamWriter tunedWriter{"output.bam", header, cfg};

// With index callback (for building BAI alongside writing)
BamWriter indexedWriter{"output.bam", header, cfg, [](std::int64_t vOffset,
                                                      std::span<const std::byte> raw) {
    // track offsets for index building
}};
```

Index callback notes:
- Callback runs asynchronously on the BGZF IO writer thread.
- Callback can be invoked after `Write()` returns (often during `Close()`).
- Call `Close()` explicitly to surface write/close errors (destructor close suppresses errors).

### Writing records

```cpp
writer.Write(view);    // from RawRecord
writer.Write(record);  // serializes from BamRecord
writer.Write(rawData); // raw span<const byte>
writer.WriteBatch(batch);  // batch of records
writer.Close();        // flush and finalize (also called by destructor)
```

### Metrics

```cpp
WriterMetrics m = writer.GetMetrics();
// m.Bgzf.BytesCompressed, m.Bgzf.BlocksWritten, m.TotalRecordsWritten, ...
```

---

## SamWriter

```cpp
#include <pbsamoa/io/SamWriter.hpp>
```

Writes text SAM files.

```cpp
SamWriter writer{"output.sam", header, SamWriterConfig{
    .UseTempFile = true,
}};
writer.Write(view);         // accepts RawRecord
writer.Write(record);       // accepts BamRecord
writer.WriteBatch(batch);   // accepts RawRecordBatch
writer.Close();
```

---

## RawRecord

```cpp
#include <pbsamoa/core/RawRecord.hpp>
```

Owning type that copies raw BAM bytes. CIGAR ops are eagerly copied into
an aligned buffer on construction. Other fields decode on demand — all
decode logic is inline in the header.

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

### Derived fields

|       Method        | Return type |             Description             |
| ------------------- | ----------- | ----------------------------------- |
| `IsMapped()`        | `bool`      | `(Flag() & 4) == 0`                 |
| `IsReverseStrand()` | `bool`      | `(Flag() & 16) != 0`                |
| `IsPrimary()`       | `bool`      | Not secondary and not supplementary |
| `ReferenceEnd()`    | `int32_t`   | `Pos() + ReferenceLength(Cigar())`  |
| `MutableTags()`     | `TagMap&`   | Mutable access to the tag map       |

### Clipping

```cpp
#include <pbsamoa/core/BamRecord.hpp>    // ClipType
#include <pbsamoa/core/TagClipping.hpp>  // TagClipper

// Clip in-place to reference coordinates
record.Clip(ClipType::CLIP_TO_REFERENCE, 1000, 2000);

// Clip with PacBio tag clipping strategies
auto clipper = TagClipper::PacBioDefault();
record.Clip(ClipType::CLIP_TO_REFERENCE, 1000, 2000, clipper);

// Non-mutating variant
BamRecord clipped = record.Clipped(ClipType::CLIP_TO_QUERY, 10, 500);
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

## CraiIndex

```cpp
#include <pbsamoa/index/CraiIndex.hpp>
```

CRAM index: read and query gzip-compressed CRAI entries.

```cpp
const CraiIndex index = CraiIndex::FromFile("input.cram.crai");
const auto& all = index.Entries();
const auto chr1 = index.EntriesForReference(0);
const auto unmapped = index.EntriesForReference(-1);
```

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

// ReadBlock returns std::optional<size_t>:
//   - optional{N>0}: N decompressed bytes
//   - optional{0}: EOF reached
//   - nullopt: error
std::optional<std::size_t> n;
while ((n = reader.ReadBlock(buf)) && *n > 0) {
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
BgzfWriterConfig cfg{
    .CompressionLevel = 6,
    .BgzfWorkers = 4,
    .InputQueueCapacity = 256,
    .BlocksPerBatch = 32,
    .UseTempFile = true,
};
BgzfWriter writer{"output.gz", cfg};
writer.Write(data);
writer.Close();  // call explicitly to observe errors
```

Close semantics:
- `Close()` transitions the writer to a terminal state; later `Write()` calls throw.
- If `Close()` throws, treat the writer as closed/failed and do not reuse it.
- Destructor still attempts close, but suppresses exceptions.

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
| `SkipPolicy`  | Silently skips. `SkippedCount()` returns total skipped. |

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

### Loading

|   Method    |               Description               |
| ----------- | --------------------------------------- |
| `Open(bam)` | Auto-detect `.zmi` then `.pbi` from BAM |
| `FromZmi()` | Load from `.zmi` file directly          |
| `FromPbi()` | Load from legacy `.pbi` file            |

### Queries

|             Method              |      Return type      |                Description                |
| ------------------------------- | --------------------- | ----------------------------------------- |
| `Find(int32_t zmw)`             | `vector<int64_t>`     | Offsets for hole number (all read groups) |
| `Find(ZmwIdentity)`             | `vector<int64_t>`     | Offsets for exact (rgId, zmw)             |
| `Find(span<const ZmwIdentity>)` | `vector<int64_t>`     | Batch query                               |
| `FirstOffset(int32_t zmw)`      | `int64_t`             | First offset for hole number              |
| `FirstOffset(ZmwIdentity)`      | `int64_t`             | First offset for exact (rgId, zmw)        |
| `UniqueZmws()`                  | `vector<ZmwIdentity>` | All unique (rgId, zmw) in file order      |
| `NumRecords()`                  | `uint64_t`            | Total records in the index                |
| `NumZmws()`                     | `uint64_t`            | Unique ZMW count                          |

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
uses a background producer thread to prefetch complete ZMW groups into a
bounded queue. Assumes the source yields records grouped by ZMW
(standard PacBio BAM ordering).

### Configuration

```cpp
struct BamZmwReaderConfig
{
    BamRecordReaderConfig Reader{};
    std::size_t PrefetchCapacityZmws{4};
};
```

Use `Reader` for the underlying `BamRecordReader` settings and
`PrefetchCapacityZmws` for the ZMW-group queue depth.

`PrefetchCapacityZmws` is measured in complete prefetched ZMW groups, not
records or bytes.

### Construction

```cpp
// Takes ownership of BamRecordReader
BamZmwReader reader{
    BamRecordReader{"movie.bam"},
    BamZmwReaderConfig{.PrefetchCapacityZmws = 4},
};

// With custom config
BamZmwReader reader{
    BamRecordReader{"movie.bam", BamRecordReaderConfig{.DecodeWorkers = 0}},
    BamZmwReaderConfig{.PrefetchCapacityZmws = 2},
};

// Convenience constructor from paths
BamZmwReader reader{
    std::vector<std::filesystem::path>{"movie1.bam", "movie2.bam"},
    BamZmwReaderConfig{
        .Reader = BamRecordReaderConfig{.DecodeWorkers = 0},
        .PrefetchCapacityZmws = 2,
    },
};
```

### Reading

```cpp
const SamHeader& header = reader.Header();

std::vector<BamRecord> records;
while (reader.GetNext(records)) {
    ZmwIdentity zmw = reader.CurrentZmw();
    // records contains all alignments for one ZMW
}
```

### Metrics

```cpp
ZmwReaderMetrics m = reader.GetMetrics();
double fullness = static_cast<double>(m.QueueDepth) / m.ConfiguredCapacity;
```

`ZmwReaderMetrics` fields:

- Capacity/occupancy: `ConfiguredCapacity`, `QueueDepth`, `PeakQueueDepth`
- Throughput: `GroupsProduced`, `GroupsConsumed`
- Stalls: `ProducerStalls`, `ConsumerStalls`

## ZmiBamWriter

```cpp
#include <pbsamoa/io/ZmiBamWriter.hpp>
```

Writes BAM and `.zmi` index simultaneously.

```cpp
ZmiBamWriter writer{"output.bam", header};
writer.Write(record);       // updates both BAM and ZMI
writer.WriteBatch(batch);   // batch write to both BAM and ZMI
writer.Close();
```

---

## ZmwWhitelist

```cpp
#include <pbsamoa/index/ZmwWhitelist.hpp>
```

A set of ZMW hole numbers or identities to select from a BAM file.

```cpp
// From bare hole numbers (matches any read group)
ZmwWhitelist wl{std::vector<std::int32_t>{12345, 23456}};

// From (rgId, zmw) pairs (exact identity matching)
ZmwWhitelist wl{std::vector<ZmwIdentity>{{0, 12345}, {0, 23456}}};

// Resolve against an index (sorted, deduplicated virtual offsets)
auto offsets = wl.Resolve(index);
```

---

## ZmiWriter

```cpp
#include <pbsamoa/io/ZmiWriter.hpp>
```

Streaming writer for `.zmi` files. BGZF-compressed.

```cpp
ZmiWriter writer{"output.bam.zmi", ZmiWriterConfig{
    .UseTempFile = true,
}};
writer.AddRecord(rgId, zmw, virtualOffset);
writer.Close();
```

---

## SequenceView

```cpp
#include <pbsamoa/core/Sequence.hpp>
```

Non-owning view over 4-bit packed BAM sequence bytes. Returned by
`RawRecord::Seq()`.

|     Method      | Return type |               Description                |
| --------------- | ----------- | ---------------------------------------- |
| `operator[](i)` | `char`      | Decoded base at position `i`             |
| `Size()`        | `uint32_t`  | Sequence length                          |
| `ToString()`    | `string`    | Decode full sequence to string           |
| `WriteTo(out)`  | `void`      | Append decoded sequence to string buffer |

### Utility functions

```cpp
std::vector<std::byte> packed = PackSequence("ACGT");
std::string seq = UnpackSequence(packed, 4);
std::string rc = ReverseComplement("ACGT");
ReverseComplementInPlace(seq);
```

---

## TagClipping

```cpp
#include <pbsamoa/core/TagClipping.hpp>
```

Extensible tag clipping via strategy pattern. `TagClipper` maps tag keys
to `TagClipStrategy` subclasses.

### Built-in strategies

|            Strategy            |                                     Tags                                     |          Behavior          |
| ------------------------------ | ---------------------------------------------------------------------------- | -------------------------- |
| `SubstringClipStrategy`        | `dq`, `iq`, `mq`, `sq`, `dt`, `st`, `ip`, `pw`, `fi`, `fp`                   | Simple substring           |
| `ReverseSubstringClipStrategy` | `ri`, `rp`                                                                   | Mirrored offset substring  |
| `PulseClipStrategy`            | `pc`, `pt`, `pq`, `pv`, `pg`, `pa`, `pm`, `ps`, `pi`, `pd`, `px`, `pe`, `sf` | Pulse-space via `pc` tag   |
| `BasemodClipStrategy`          | `MM`, `ML`                                                                   | SAMv1.6 base modifications |
| `PileupClipStrategy`           | `sa`                                                                         | RLE pileup coverage        |

### Usage

```cpp
auto clipper = TagClipper::PacBioDefault();  // all standard strategies
clipper.ClipTags(tags, clipOffset, clipLength, seqLength, sequence);
```

---

## Metrics

```cpp
#include <pbsamoa/core/Metrics.hpp>
```

Runtime performance snapshots for all reader/writer pipelines. Captured
via relaxed atomics — diff two snapshots for per-second rates.

### Types

|        Type        |             Source              |
| ------------------ | ------------------------------- |
| `BgzfMetrics`      | `BamRawReader::GetMetrics()`    |
| `DecodeMetrics`    | (internal to `ReaderMetrics`)   |
| `ReaderMetrics`    | `BamRecordReader::GetMetrics()` |
| `ZmwReaderMetrics` | `BamZmwReader::GetMetrics()`    |
| `BgzfWriteMetrics` | (internal to `WriterMetrics`)   |
| `WriterMetrics`    | `BamWriter::GetMetrics()`       |

### BgzfMetrics fields

Throughput: `BytesRead`, `BlocksRead`, `BytesDecompressed`.
Pool (via `PoolMetrics Pool`): `QueueDepth`, `PeakQueueDepth`, `ActiveWorkers`,
`PeakActiveWorkers`, `ResultQueueDepth`, `PeakResultQueueDepth`.
Queue: `RecordsProduced`, `RecordsConsumed`.
Stalls: `IoStalls`, `ConsumerStalls`, `ReaderStalls`.
Timing (ns): `IoReadNs`, `DecompressNs`, `RecordParseNs`.

### ZmwReaderMetrics fields

Capacity/occupancy: `ConfiguredCapacity`, `QueueDepth`, `PeakQueueDepth`.
Throughput: `GroupsProduced`, `GroupsConsumed`.
Stalls: `ProducerStalls`, `ConsumerStalls`.
