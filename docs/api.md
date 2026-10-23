# API Reference

All types live in `namespace PacBio::Samoa`. Headers are under
`include/pbsamoa/`, organized into `core/`, `io/`, `index/`, and `cram/`
subdirectories. `PbSamoaLibraryInfo.hpp` is at the top level.

The snippets below are fragments. They assume `using namespace PacBio::Samoa;`.
Names such as `header`, `record`, `view`, and `batch` stand for objects of the types
that the section describes.

---

## BamRawReader

```cpp
#include <pbsamoa/io/BamRawReader.hpp>
```

Reads BAM files. Supports sequential iteration, batch reading, and BAI
region queries. `ReadRecord()` and `Records()` give owning `RawRecord`s.
`ReadBatch()` with `View(i)` gives non-owning views and does not copy the records.

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

// Scatter chunk 2 of 4: a chaotic-but-deterministic sample across the whole file
BamRawReader reader{"input.bam", BamRawReaderConfig{.ChunkNum = 2,
                                                    .TotalChunks = 4,
                                                    .ChunkingMode = ChunkMode::SCATTER,
                                                    .ChunkTileZmws = 1,
                                                    .ChunkSeed = 0}};

// Several BAM files as one stream (the headers must be compatible)
BamRawReader reader{std::vector<std::filesystem::path>{"a.bam", "b.bam"}};
```

The constructor also accepts a `BamCollection` or a `std::vector<BamFile>` (see
[BamCollection and BamFile](#bamcollection-and-bamfile)).

`BamRawReaderConfig` fields:

|      Field      |   Default    |                    Description                    |
| --------------- | ------------ | ------------------------------------------------- |
| `BgzfWorkers`   | `0`          | BGZF decompression threads (0 = synchronous)      |
| `RecordLimit`   | `0`          | Stop after N records (0 = unlimited)              |
| `ChunkNum`      | `0`          | 1-based chunk number (0 = no chunking)            |
| `TotalChunks`   | `0`          | Total chunk count (0 = no chunking)               |
| `ChunkingMode`  | `CONTIGUOUS` | `CONTIGUOUS` ZMW ranges, or `SCATTER` (see below) |
| `ChunkTileZmws` | `100`        | Scatter only: read up to M consecutive ZMWs/seek  |
| `ChunkSeed`     | `42`         | Scatter only: deterministic shuffle seed          |
| `Whitelist`     | `{}`         | `std::optional<ZmwWhitelist>` (see note below)    |

`Whitelist` does not filter `ReadRecord()`, `Records()`, or `ReadBatch()`. Use
`reader.Whitelist(wl)` to select ZMWs. The constructor throws `std::invalid_argument`
if `Whitelist` is set together with chunking or with several input files.
Chunking also needs one input file. Set both `ChunkNum` and `TotalChunks` or neither:
if either is set, the constructor throws `std::invalid_argument` unless
`1 <= ChunkNum <= TotalChunks`.

With `CONTIGUOUS` chunking, the chunk's record count replaces `RecordLimit`. In
`SCATTER` mode, `RecordLimit` caps the records of the chunk.

Both chunk modes are ZMW-keyed (a ZMW's records stay together) and partition the
file exactly: all `TotalChunks` chunks together visit every record once.
`CONTIGUOUS` gives each chunk one contiguous ZMW range. `SCATTER` assigns tiles
of up to `ChunkTileZmws` consecutive ZMWs to chunks via a seeded balanced
shuffle, so each chunk samples ZMWs across the whole file; the same
`(file, ChunkNum, TotalChunks, ChunkTileZmws, ChunkSeed)` reproduces the same
chunk on every platform. Scatter reads synchronously (it ignores `BgzfWorkers`);
both `Records()` and `ReadBatch` honor the scatter plan, so `BamRecordReader`
works in scatter mode too.

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

Returns references to iterator-owned `RawRecord`s. Each reference is valid until the next
iteration step.

### Batch reading

```cpp
using namespace PacBio::Samoa::Literals;

while (auto batch = reader.ReadBatch(128_MiB)) {
    for (std::size_t i{0}; i < batch->RecordCount(); ++i) {
        const RawRecordView view{batch->View(i)};
        // process
    }
}
```

Returns `std::optional<RawRecordBatch>`. `nullopt` at EOF. The default limit is
`ByteLimit{}` (256 MiB).

With `BgzfWorkers == 0`, the batch copies the records and grows up to the limit.
With `BgzfWorkers > 0`, a batch holds the records of at most one group of
decompressed BGZF blocks (about 2 MiB). If the group fits in the limit, the batch
takes over the group buffer and does not copy the records. A larger limit does not
join groups, so a large limit gives batches of about 2 MiB.

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
would restart per seek). `Whitelist()` needs a single input file. It throws
`std::logic_error` for a multi-file reader, and `std::runtime_error` when the BAM has
no `.zmi` or `.pbi`. Records come in file order, not whitelist order.

### ZMW count

```cpp
const std::int32_t numZmws = reader.NumZmws();  // -1 if not available
```

With chunking, `NumZmws()` gives the ZMW count of the chunk (for an empty
`CONTIGUOUS` chunk, the count of the whole file). Otherwise it opens the `.zmi` or
`.pbi` of the BAM. It gives `-1` for multi-file input or when no index loads.

### Random access

```cpp
VirtualOffset pos = reader.Tell();
// ... read some records ...
reader.Seek(pos);  // jump back
```

`Seek`, `Tell`, and `Query(index, ...)` throw `std::logic_error` for multi-file input.

### Metrics

```cpp
BgzfMetrics m = reader.GetMetrics();
// m.BytesRead, m.BlocksRead, m.BytesDecompressed, m.RecordsConsumed, ...
```

With `BgzfWorkers == 0`, only `RecordsConsumed` has a value.

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

As `BamRawReader`, the constructor also accepts a `BamCollection`, a
`std::vector<std::filesystem::path>`, or a `std::vector<BamFile>`. `Header()` and
`NumZmws()` behave as in `BamRawReader`.

### Configuration

|       Field       |                   Type                   | Default |                     Description                      |
| ----------------- | ---------------------------------------- | ------- | ---------------------------------------------------- |
| `RawReaderConfig` | `BamRawReaderConfig`                     | `{}`    | Config for underlying `BamRawReader`                 |
| `DecodeWorkers`   | `size_t`                                 | `4`     | Threads for parallel `ToOwned()` decode (0 = serial) |
| `BatchBudget`     | `ByteLimit`                              | `4_MiB` | Memory budget per batch read from view reader        |
| `OutputCapacity`  | `size_t`                                 | `1024`  | SPSC queue capacity (records)                        |
| `TagFilter`       | `variant<monostate, DropTags, KeepTags>` | `{}`    | Tag filter applied during `ToOwned()`                |

### Reading records

```cpp
while (auto record = reader.ReadRecord()) {
    // record is an owned BamRecord
}
```

When the producer thread fails, `ReadRecord()` rethrows its exception with the original
type.

### Range iteration

```cpp
#include <print>

for (const auto& record : reader.Records()) {
    std::println("{}", record.Name());
}
```

### Region query

```cpp
for (const BamRecord& record : reader.Query("chr1", 1000, 2000)) {
    // decoded records that overlap chr1:[1000, 2000)
}
const GenomicInterval interval{"chr1", 1000, 2000};
for (const BamRecord& record : reader.Query(interval)) {
}
```

`Query` opens the sibling `.bai` of each input BAM (`<input>.bai`) and its own readers.
It does not change the sequential iteration. With several inputs, it merges the
records of each file in coordinate order. It throws `std::runtime_error` when the
reference name is not in the header or a `.bai` file is missing. The tag filter of the
reader applies; the rest of `RawReaderConfig` (BGZF workers, chunking, `RecordLimit`)
does not, since the query readers use the default config.

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

|         Field          | Default |                                                                                     Description                                                                                      |
| ---------------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `ReferencePath`        | `{}`    | Optional FASTA path for reference-based slices                                                                                                                                       |
| `Reference`            | `{}`    | Optional `CramReference` for reference-based slices, read in windows                                                                                                                 |
| `DecompressionWorkers` | `0`     | Workers that decode up to 2 containers each ahead of `ReadRecord`/`ReadRawRecord`, in file order; for `Query`/`ReadSlice`, they decompress the blocks of a slice (`0` = synchronous) |

With `ReferencePath`, the reader opens the file as a [`FastaFile`](#fastafile): it maps the file and
reads only the bases it needs. With `Reference` (a
`shared_ptr<const CramReference>`, one sequence per @SQ line, in header order), the reader keeps
the bases of one window only: the span of the current slice, or 64 KiB or more for a record of a
multi-reference slice. Several readers can share one `Reference`. The reader makes the bases upper
case. Set `ReferencePath` or `Reference`, not both: the constructor throws
`std::invalid_argument`, and also when `Reference` does not have one sequence per @SQ line.
With neither set, the reader uses the FASTA named by the environment variable
`PBSAMOA_CRAM_REFERENCE`, if set. FASTA sequences are matched to @SQ lines by name.

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

A query reads on the stream of the reader and then puts the stream back, so the sequential
iteration (`ReadRecord`, `Records()`) goes on after it unchanged. A query does not open the file
again, parse the header again or load the reference again.

To decode each slice once for many queries, take the slices of a region from the index and
decode each one whole:

```cpp
for (const CraiEntry& row : index.Query(refId, 1000, 2000)) {  // one row per slice, file order
    const std::vector<RawRecord> records = reader.ReadSliceRaw(row);  // or ReadSlice
}
```

`ReadSlice` reads on the stream of the reader as `Query` does. A multi-reference slice also has
the records of other references.

---

## CramWriter

```cpp
#include <pbsamoa/io/CramWriter.hpp>
```

Writes CRAM v3.0/v3.1 files from `BamRecord`, `RawRecordView`, and `RawRecord` input.

### Construction

```cpp
// Default config
CramWriter writer{"output.cram", header};

// Tuned config with CRAI sidecar
CramWriter writer{"output.cram", header, CramWriterConfig{
    .BlockCompressionMethod = CramBlockMethod::RANS4X8,
    .RecordsPerSlice = 10000,
    .SlicesPerContainer = 8,
    .WriteCrai = true,
    .CompressionWorkers = 4,
    .CompressionLevel = 1,  // gzip levels: [0, 12]
    .UseTempFile = true,
}};
```

`CramWriterConfig` fields:

|             Field              |                       Type                       | Default |                            Description                            |
| ------------------------------ | ------------------------------------------------ | ------- | ----------------------------------------------------------------- |
| `BlockCompressionMethod`       | `CramBlockMethod`                                | `GZIP`  | Default block method                                              |
| `DataSeriesCompressionMethods` | `unordered_map<CramDataSeries, CramBlockMethod>` | `{}`    | Per-data-series method overrides                                  |
| `RecordsPerSlice`              | `int32_t`                                        | `10000` | Max records per slice                                             |
| `BasesPerSlice`                | `int64_t`                                        | `20000000` | Max bases per slice: sequence bases plus BAM aux bytes (`0` = no limit) |
| `SlicesPerContainer`           | `int32_t`                                        | `1`     | Max slices per container (`<= 0` means unlimited until `Close()`) |
| `MajorVersion`                 | `uint8_t`                                        | `3`     | CRAM major version (must be `3`)                                  |
| `MinorVersion`                 | `uint8_t`                                        | `0`     | CRAM minor version (`0` or `1`)                                   |
| `WriteCrai`                    | `bool`                                           | `false` | Write `<output>.crai` sidecar                                     |
| `CraiPath`                     | `optional<path>`                                 | `{}`    | CRAI output path override                                         |
| `CompressionWorkers`           | `size_t`                                         | `0`     | CRAM encode threads: with 2 or more, up to this many containers encode at once, in file order (`0` = synchronous) |
| `CompressionLevel`             | `optional<int>`                                  | `{}`    | Optional gzip compression level (`[0,12]`)                        |
| `UseTempFile`                  | `bool`                                           | `false` | Write to a temporary path and rename atomically on `Close()`      |
| `ReferenceSequences`           | `vector<string_view>`                            | `{}`    | The sequence of each @SQ line: reference-based CRAM (see below)   |
| `Reference`                    | `shared_ptr<const CramReference>`                | `{}`    | The reference in windows, in place of `ReferenceSequences`        |

Notes:
- `MinorVersion` is automatically raised to `1` when selected codecs require CRAM 3.1.
- `CompressionLevel` is currently applied to gzip-compressed blocks; non-gzip codecs ignore it.
- A slice ends with the record that reaches `RecordsPerSlice` or `BasesPerSlice`, so it holds at
  least one record. The bases limit counts as htslib's `bases_per_slice` does; htslib's default is 5000000.
- The constructor throws `std::runtime_error` when `MajorVersion != 3`, `RecordsPerSlice <= 0`,
  `BasesPerSlice < 0`, `CompressionLevel` is outside `[0, 12]`, or a `ReferenceSequences` entry does not have the
  length of its @SQ line.

### Reference-based CRAM

With `ReferenceSequences` (one view per @SQ line, in header order, of the @SQ length; any case;
valid until `Close()`), a mapped record stores only its differences to the reference: a
substitution (`X`, default substitution matrix) where both bases are in ACGTN, else a run of bases
(`b`), for example an IUPAC code or a base past the end of the reference. The preservation map
says `RR`, each single-reference slice carries the MD5 of its reference bases, and @SQ lines
without `M5` get one (`CramReferenceMd5(sequence)`: the MD5 of the bases in upper case). Readers
need the reference (`CramReaderConfig::ReferencePath`, `samtools view -T`). As with every CRAM,
readers give `M` for `=` and `X`, and htslib adds `MD` unless `decode_md=0`.

`Reference` gives the same file without the text of the whole reference in memory: a
`CramReference` gives the number of sequences, their lengths, and the bases of a window
(`Bases(refId, start, length, buffer)`, called from several threads; the result can be a view of
`buffer`). The writer asks for the bases under each mapped record and of each single-reference
slice. `CramReferenceMd5(reference, refId)` gives the M5 from windows of 64 KiB. Set
`ReferenceSequences` or `Reference`, not both.

### Encoding containers on several threads

```cpp
// any thread: one container (slices of RecordsPerSlice records or BasesPerSlice bases)
CramEncodedContainer encoded = writer.EncodeContainer(std::span<const BamRecord>{records});
// one thread, in output order: record counters and CRAI offsets continue
writer.WriteContainer(std::move(encoded));
```

`EncodeContainer` uses the calling thread only (thread-local compressors, not the compression
workers), so several threads can encode while another writes. `WriteContainer` first writes the
records given one by one before it.

To also serialize on the worker threads, use `Serialize` and `WriteSerialized`:

```cpp
// any thread: the container follows recordCounter records in the file
CramSerializedContainer bytes = writer.Serialize(writer.EncodeContainer(records), recordCounter);
// one thread, in output order
writer.WriteSerialized(std::move(bytes));
```

`CramSerializedContainer` holds `Bytes`, `RecordCounter`, `NumRecords`, and the CRAI rows
(`Index`). `WriteSerialized` throws `std::runtime_error` if `RecordCounter` is not the
number of records that the writer wrote before it.

### Containerization and CRAI behavior

- Records are buffered into slices of up to `RecordsPerSlice` records and `BasesPerSlice` bases.
- Slices are buffered into containers of up to `SlicesPerContainer`.
- With `SlicesPerContainer <= 0`, containers are flushed only at `Close()`.
- With `WriteCrai = true`, CRAI rows are emitted per serialized slice offset.
- Mixed-reference slices (`RefSeqId == -2`) emit one CRAI row per referenced sequence,
  plus one unmapped row if unmapped records are present.

### Writing records

```cpp
std::vector<BamRecord> records{/* ... */};
for (const auto& record : records) {
    writer.Write(record);
}
BamRecord temporary{};
writer.Write(std::move(temporary));  // BamRecord&& overload
writer.Write(rawRecord);             // RawRecord or RawRecordView
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
- `ParseFileDefinition`/`SerializeFileDefinition`, `ParseContainerHeader`/`SerializeContainerHeader`,
  `ParseBlock`/`SerializeBlock`, `ParseSliceHeader`/`SerializeSliceHeader`, and
  `ParseCompressionHeader`/`SerializeCompressionHeader` read and write each part.
- `ReadItf8`/`WriteItf8` and `ReadLtf8`/`WriteLtf8` encode the CRAM variable-length integers.
- `IsCramEofMarker(bytes)` checks for `CRAM_EOF_MARKER`. `CRAM_MAGIC` is the file magic.
- Enums: `CramBlockMethod` (`RAW`, `GZIP`, `BZIP2`, `LZMA`, `RANS4X8`, `RANS4X16`,
  `ADAPTIVE_ARITH`, `FQZCOMP`, `NAME_TOKENISER`), `CramBlockContentType`, `CramCodecId`, and
  `CramDataSeries` (the two-letter data-series keys, for example `CramDataSeries::QS`).
- `CramFileDefinition`, `CramContainerHeader`, `CramBlock`, `CramSliceHeader`,
  `CramPreservationMap`, `CramEncodingDescriptor`, and `CramCompressionHeader` are the parts.

`ParseContainer` parses and decompresses compression/slice headers, while
core/external slice data blocks remain in their on-disk compressed form.

```cpp
std::size_t used{0};
const CramContainerHeader containerHeader = ParseContainerHeader(bytes, used);
const CramContainer container =
    ParseContainer(containerHeader, bytes.subspan(used, containerHeader.Length));
```

`CramEncodedContainer` (the result of `CramWriter::EncodeContainer`) holds a `CramContainer`
and the CRAI rows of each slice (`SliceIndex`).

### Block compression

```cpp
#include <pbsamoa/cram/CramCompression.hpp>
```

- `CompressCramBlock(block, method, level)` and `DecompressCramBlock(block)` change a
  `CramBlock` in place.
- `CompressionRegistry::Instance()` holds the codecs of the block methods 0 to 8 (`BZIP2`
  and `LZMA` only when the library is built with bzip2 and liblzma). `Register`
  adds a custom method ID above 8. It throws `std::runtime_error` for a built-in ID.
  Do not register a method while other threads use the registry.
- `CramGzipCompress`, `CramGzipDecompress`, and `CramFqzcompCompress` call one codec directly.

This header includes `<libdeflate.h>`, because some overloads take a reusable
`libdeflate_compressor*` or `libdeflate_decompressor*`.

`<pbsamoa/cram/CramCodec.hpp>` holds the CRAM data-series codecs (`CramCodec`, `CreateCodec`,
`NullCodec`, `ExternalCodec`, `HuffmanCodec`, `BetaCodec`, `GammaCodec`, `SubexpCodec`,
`ByteArrayLenCodec`, `ByteArrayStopCodec`), their `CramCodecDecodeKind` (`SCALAR` or
`BYTE_ARRAY`), and the streams that they use (`CramBitReader`, `CramBitWriter`,
`CramExternalBlockStore`). The reader and the
writer use these types. Most programs do not need them.

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

`Header()` gives the parsed header. `ReadRecord()` gives the next record, or
`std::nullopt` at the end of the file.

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
- Call `Close()` explicitly to surface write/close errors. The destructors of all writers
  (`BamWriter`, `SamWriter`, `CramWriter`, `BgzfWriter`) close and suppress errors.
- The callback supplies each record's begin virtual offset; after `Close()`,
  `writer.EndVirtualOffset()` gives the offset just past the last record (the end offset
  of the final record). Feed both to `BaiStreamBuilder` to build a BAI on the fly — see
  `BaiIndex` below.

### Writing records

```cpp
writer.Write(view);    // from RawRecordView or RawRecord
writer.Write(record);  // serializes from BamRecord
writer.Write(rawData); // raw span<const byte>
writer.WriteBatch(batch);  // batch of records
writer.Close();        // flush and finalize (also called by destructor)
```

`rawData` holds the record bytes without the 4-byte `block_size`; the writer prepends it.
`Write()` throws `std::invalid_argument` for a record of more than `INT32_MAX` bytes.

### Metrics

```cpp
WriterMetrics m = writer.GetMetrics();
// m.Bgzf.BytesCompressed, m.Bgzf.BlocksWritten, m.TotalRecordsWritten, ...
```

### Compression on the caller threads

`BgzfBlockEncoder` and `BgzfBlockWriter` (`<pbsamoa/core/Bgzf.hpp>`) move the
compression to the threads that make the records. Each thread fills its own
encoder; the writer only appends the blocks, in the order of the `Write()` calls.

```cpp
BgzfBlockWriter out{"output.bam", BgzfBlockWriterConfig{.UseTempFile = true}};

BgzfBlockEncoder headerEncoder{};
headerEncoder.Write(header.ToBamHeaderBlock());
headerEncoder.Finish();
out.Write(headerEncoder.TakeBlocks());

// On each worker thread, for one batch:
BgzfBlockEncoder encoder{/*compressionLevel=*/2};
for (const BamRecord& record : batch) {
    EncodeBamRecord(record, encoder);
}
encoder.Finish();
out.Write(encoder.TakeBlocks());  // thread-safe

out.Close();  // writes the EOF marker
```

`EncodeBamRecord` (`<pbsamoa/io/BamWriter.hpp>`) appends the `block_size` field and the
record bytes. It throws `std::invalid_argument` if the record is too large for BAM.
`BgzfBlockEncoder` accepts levels 0 to 12 and throws `std::invalid_argument` for other
values. Level 0 selects a fast lossless mode, not stored blocks.

`Write()` blocks while the queue holds `MaxPendingBytes` (default 64 MiB) that
the IO thread did not write yet. Thus slow storage cannot make the queue grow
without limit. The writer does not track virtual offsets and has no index
callback. Build a BAI from the finished file.

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
writer.Write(view);         // accepts RawRecordView or RawRecord
writer.Write(record);       // accepts BamRecord
writer.WriteBatch(batch);   // accepts RawRecordBatch
writer.Close();
```

---

## RawRecordView and RawRecord

```cpp
#include <pbsamoa/core/RawRecord.hpp>
```

`RawRecordView` borrows raw BAM bytes. The source buffer must outlive the view and
the spans that its accessors return. `RawRecord` provides the same accessor API
but owns a byte copy, so it can cross threads and outlive the source. `RawRecord`'s
explicit `View()` accessor exposes a non-owning view and keeps the lifetime relationship
visible. You cannot call `View()` on a temporary `RawRecord`.
Both types copy CIGAR ops into aligned storage because BAM does not
guarantee their alignment. Other fields decode on demand.

Construct either type from the record bytes after `block_size` (the span starts at
`refID`). The constructor throws `std::invalid_argument` for a malformed record. Examples
are fewer than 32 bytes, `l_read_name` of 0, or a name without NUL. Fields past the end of
the span, or a mapped CIGAR whose query length is not `l_seq`, also cause it.
`RawData()` gives these bytes again. If the CIGAR is in a `CG` tag (more than 65535
operations), the record expands it, and `ParseTags()` leaves out `CG`.

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
| `CigarOps()` | `CigarView`           | Span of CIGAR operations (see below)     |
| `Seq()`      | `SequenceView`        | 4-bit packed sequence (decode on access) |
| `Qual()`     | `span<const uint8_t>` | Base qualities                           |
| `AuxData()`  | `span<const byte>`    | Raw auxiliary tag bytes                  |

The span from `CigarOps()` points into the view's own storage, so it becomes invalid when
the view or `RawRecord` is moved or destroyed. Call `CigarOps()` where you use it. The
other spans of a `RawRecord` survive a move.

### Derived fields

|       Method        | Return type |              Description              |
| ------------------- | ----------- | ------------------------------------- |
| `IsMapped()`        | `bool`      | `(Flag() & 4) == 0`                   |
| `IsReverseStrand()` | `bool`      | `(Flag() & 16) != 0`                  |
| `IsPrimary()`       | `bool`      | Not secondary and not supplementary   |
| `IsSecondary()`     | `bool`      | `(Flag() & 0x100) != 0`               |
| `IsSupplementary()` | `bool`      | `(Flag() & 0x800) != 0`               |
| `ReferenceLength()` | `int64_t`   | Consumed reference bases (from CIGAR) |
| `QueryLength()`     | `int64_t`   | Consumed query bases (from CIGAR)     |

### Conversion

```cpp
BamRecord record = view.ToOwned();
BamRecord noKinetics = view.ToOwned(DropTags{TagKey{'i','p'}, TagKey{'p','w'}});
BamRecord rgOnly = view.ToOwned(KeepTags{TagKey{'R','G'}});
TagMap tags = view.ParseTags();
TagMap some = view.ParseTagsFiltered(filter, keep);  // keep(filter, key) selects the tags
RawRecord raw = RawRecord::FromBamRecord(record);  // serialized in place, no second copy
```

`RawRecord::SetNext(nextRefId, nextPos)` and `SetTlen(tlen)` set the mate fields of an owned
record in place (the CRAM reader resolves attached mates on raw records).

---

## BamRecord

```cpp
#include <pbsamoa/core/BamRecord.hpp>
```

Owned, mutable record. Fields stored as C++ types. Mutations are plain
assignments; serialization happens once in the writer.

### Accessors

`Name()`, `Flag()`, `RefId()`, `Pos()`, `MapQ()`, `NextRefId()`,
`NextPos()`, `Tlen()` have the same names as `RawRecord`. Four
accessors differ:

|  `BamRecord`  | `RawRecord`   |
| ------------- | ------------- |
| `Cigar()`     | `CigarOps()`  |
| `Sequence()`  | `Seq()`       |
| `Qualities()` | `Qual()`      |
| `Tags()`      | `ParseTags()` |

`Sequence()` gives a `string_view` of the bases. `Qualities()` gives a
`span<const uint8_t>` and is empty when the record has no qualities. `Tags()` gives a
`const TagMap&`. A default-constructed record has `RefId() == -1` and `Pos() == -1`.
`Name(name)` throws `std::invalid_argument` for a name of more than 254 characters.

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
| `IsSecondary()`     | `bool`      | `(Flag() & 0x100) != 0`             |
| `IsSupplementary()` | `bool`      | `(Flag() & 0x800) != 0`             |
| `ReferenceEnd()`    | `int32_t`   | `Pos() + ReferenceLength(Cigar())`  |
| `AlignedStart()`    | `int32_t`   | First aligned query position (-1)   |
| `AlignedEnd()`      | `int32_t`   | One past last aligned position (-1) |
| `MutableTags()`     | `TagMap&`   | Mutable access to the tag map       |

`AlignedStart()` and `AlignedEnd()` use polymerase coordinates: `QueryStart()` or
`QueryEnd()` corrected by the soft clips. They give -1 if a hard clip hides part of the query.

### Mapping

```cpp
// Make an unmapped record a mapped record on the reverse strand
record.Map(/*refId=*/0, /*pos=*/1000, /*reverse=*/true, {CigarOp{CigarOpType::M, 100}},
           /*mapQ=*/60);
```

`Map` clears the unmapped flag and sets the position, CIGAR, MAPQ, and strand. On the
reverse strand, it reverse-complements SEQ and reverses QUAL. It does not change per-base
tags such as `ip`, `pw`, and `ML`. The record must have forward-strand SEQ and QUAL.

### PacBio accessors

```cpp
std::string movie = record.MovieName();   // "movie" of "movie/123/0_1000"
std::int32_t zmw = record.HoleNumber();   // 123
std::int32_t qs = record.QueryStart();    // qs tag, else from the name, else 0
std::int32_t qe = record.QueryEnd();      // qe tag, else from the name, else sequence length
std::string rg = record.ReadGroupId();    // RG tag
std::optional<PacBio::Data::Frames> pw = record.PulseWidth();  // nullopt if absent
```

`HoleNumber()`, `ReadGroupId()`, `SignalToNoise()`, and `ReadAccuracy()` throw
`std::runtime_error` when the name or the tag is missing or malformed. The optional
accessors give `nullopt` for a missing tag but throw `std::runtime_error` for a tag of the
wrong type: `QueryStart()`, `QueryEnd()`, `WallStart()`, `WallEnd()`, and
`LocalContextFlags()` for a non-integer tag, `PulseWidth()` and `IPD()` for a tag that is not
an array of a supported element type.

| Method                  | Return type                    | Tag  |
| ----------------------- | ------------------------------ | ---- |
| `FullName()`            | `string`                       |      |
| `LocalContextFlags()`   | `optional<LocalContextFlags>`  | `cx` |
| `LocalContextFlags(f)`  | `void` (sets the tag)          | `cx` |
| `SignalToNoise()`       | `Data::SNR`                    | `sn` |
| `ReadAccuracy()`        | `Data::Accuracy`               | `rq` |
| `PulseWidth()`, `IPD()` | `optional<Data::Frames>`       | `pw`, `ip` |
| `WallStart()`, `WallEnd()` | `optional<int32_t>`         | `ws`, `we` |

The `Data::` types come from pbcopper (`<pbcopper/data/...>`). Thus programs that include
`BamRecord.hpp`, `SamHeader.hpp` (`Data::Strand`), or headers that include them, such as
`Bgzf.hpp`, also need the pbcopper headers.

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

// Also remove insertions at the new ends of the alignment
record.Clip(ClipType::CLIP_TO_REFERENCE, 1000, 2000, /*exciseFlankingInserts=*/true);
```

`CLIP_TO_QUERY` uses polymerase (ZMW) coordinates. `CLIP_TO_REFERENCE` uses reference
coordinates. Each overload takes an optional `bool exciseFlankingInserts` (default `false`)
as its last argument. Only `CLIP_TO_REFERENCE` uses it. A reference clip of an unmapped
record does nothing. Without a clipper, `Clip` uses `TagClipper::PacBioDefault()`.
See also [CigarClipping](#cigarclipping) and [LosslessClipping](#losslessclipping).

### Serialization

```cpp
std::vector<std::byte> raw = record.SerializeToBam();  // without the block_size prefix
```

Normally you don't call this directly — `BamWriter::Write(record)` handles
it. `SerializeToBam()` and `RawRecord::FromBamRecord` throw `std::invalid_argument` when the
qualities are not empty and differ in length from the sequence, or when a CIGAR of more than
65535 operations spans 2^28 reference bases or more.

---

## RawRecordBatch

```cpp
#include <pbsamoa/core/RawRecord.hpp>
```

Owns a decompressed buffer and provides views into it. Construct it from a buffer and a
`std::vector<RawRecordBatch::RecordExtent>` (`offset` and `size` of each record).

|     Method      |    Return type     |            Description            |
| --------------- | ------------------ | --------------------------------- |
| `View(i)`       | `RawRecordView`    | Non-owning decoder for record `i` |
| `RecordData(i)` | `span<const byte>` | Raw data for record `i`           |
| `RecordCount()` | `size_t`           | Number of records                 |
| `BufferSize()`  | `size_t`           | Decompressed buffer size in bytes |

Views are invalidated when the batch is destroyed. The constructor does not check the
extents, and `View(i)` and `RecordData(i)` do not check `i`.

---

## SamHeader

```cpp
#include <pbsamoa/core/SamHeader.hpp>
```

Parsed SAM/BAM header.

### Construction

```cpp
std::expected<SamHeader, std::string> parsed = SamHeader::FromText(samText);
if (!parsed) {
    // parsed.error() describes the problem
}
SamHeader header = *parsed;

std::expected<SamHeader, std::string> fromBam = SamHeader::FromBamHeaderBlock(rawBytes);
```

Both factories return `std::expected<SamHeader, std::string>`. On bad input, `error()`
holds the message. In a BAM header block, the binary reference list is authoritative.

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

Each getter gives a `string_view`. `SetVersion`, `SetSortOrder`, `SetGroupOrder`, and
`SetSubSort` change the fields.

### Reference sequences (`@SQ`)

```cpp
header.AddReferenceSequence(ReferenceSequence{"chr1", 248956422});

std::int32_t id = header.ReferenceId("chr1");       // name → ID
std::string_view name = header.ReferenceName(id);    // ID → name
std::int32_t n = header.NumReferences();

for (const ReferenceSequence& sq : header.ReferenceSequences()) {
    const std::string* md5 = sq.GetTag("M5");  // nullptr if absent
}
```

`ReferenceSequence` has `Name()` and `Length()`. `ReferenceId` gives -1 for an unknown name. `ReferenceName(-1)` gives `"*"`.
`ReferenceName` throws `std::out_of_range` for other invalid IDs.

### Read groups, programs, comments

```cpp
header.AddReadGroup(ReadGroup{"rg1"});
header.AddProgramRecord(ProgramRecord{"pbsamoa"});
header.AddComment("processed by pbsamoa");

// Add a chained @PG; the ID gets a ".N" suffix if it is already in use
AppendProgramRecord(header, "pbsamoa", std::optional<std::string>{"pbsamoa sort in.bam"});
```

`ReadGroups()`, `ProgramRecords()`, and `Comments()` give spans. The non-const
`ReferenceSequences()`, `ReadGroups()`, and `ProgramRecords()` give the vectors for change.
`ReferenceSequence`, `ReadGroup`, and `ProgramRecord` keep optional fields as text tags:
`GetTag(key)` (nullptr if absent), `SetTag(key, value)`, and `CustomTags()`.
`ProgramRecord` also has `RemoveTag(key)`. `ReadGroup::Id()` and `ProgramRecord::Id()` give
the ID. `AppendProgramRecord` sets `PN=pbsamoa`, `VN` to the library release, `CL` when a
command line is given, and `PP` to the ID of the last @PG.

### Read group IDs

```cpp
const std::string id = MakeReadGroupId("m84001_230601_123456", "CCS");  // 8 hex digits
const std::string barcoded = MakeReadGroupId("m84001_230601_123456", "CCS", "bc1001--bc1001");
const std::string_view base = ReadGroupBaseId(barcoded);  // == id

ReadGroup rg{id};
rg.SetDsField("READTYPE", "CCS");
const std::string* readType = rg.ReadType();
```

`MakeReadGroupId` gives the first 8 hex digits of the MD5 of `"{movie}//{readType}"`. An
optional `Data::Strand` argument appends `//fwd` or `//rev`. `ReadGroup` also reads the
PacBio `DS` fields: `DsField(key)`, `SetDsField`, `ParsedDsFields`, `SetDsFields`,
`ReadType()`, `BindingKit()`, `SequencingKit()`, `BasecallerVersion()`, and `MovieName()`
(the `PU` tag).

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

auto faster = BaiIndex::Build("sorted.bam", /*numWorkers=*/4);  // same bytes
```

`Build` throws `std::runtime_error` if the BAM is not coordinate-sorted. `FromFile` and
`ToFile` throw `std::runtime_error` on I/O or format errors.

`Query(refId, beg, end)` gives the `Chunk`s (`Begin`, `End` virtual offsets) that overlap
`[beg, end)`, sorted and merged. It gives an empty vector for an unknown `refId`.
`Reference(refId)` gives the `ReferenceIndex` (`bins` and `linearIndex`) of one reference.
It throws `std::out_of_range` for an invalid `refId`.

### Build on the fly (streaming)

`BaiStreamBuilder` is the shared indexing core that `Build()` drives from a file scan.
Drive it instead from a `BamWriter` index callback to index while writing — no second
pass. Both paths produce byte-identical output. Records must arrive in coordinate order;
each record's end offset is the begin offset of the next, so `Finalize()` takes the end
offset of the last record (`BamWriter::EndVirtualOffset()`).

```cpp
#include <pbsamoa/io/BamWriter.hpp>

BaiStreamBuilder builder{header.NumReferences()};
BamWriter writer{"sorted.bam", header, cfg,
                 [&](std::int64_t beginVo, std::span<const std::byte> recordBytes) {
                     builder.Observe(VirtualOffset{static_cast<std::uint64_t>(beginVo)},
                                     recordBytes);  // runs on the IO thread
                 }};
// ... writer.Write(...) per record ...
writer.Close();                                       // joins the IO thread
builder.Finalize(writer.EndVirtualOffset()).ToFile("sorted.bam.bai");
```

`pbsamoa merge --bai` (and `MergeConfig::BaiOutput`) use exactly this path for
coordinate-sorted heap merges. Concat, queryname/tag merges, and disjoint-chain
passthrough are rejected; run `bai-build` on those outputs instead. `Observe()`
throws if the records are not coordinate-sorted.

### Statistics

|      Method       |                 Description                 |
| ----------------- | ------------------------------------------- |
| `NumReferences()` | Number of reference sequences in the index  |
| `MappedCount()`   | Total mapped reads (from metadata bin)      |
| `UnmappedCount()` | Unmapped reads: metadata bins + `n_no_coor` |

---

## CraiIndex

```cpp
#include <pbsamoa/index/CraiIndex.hpp>
```

CRAM index: read and query gzip-compressed CRAI entries.

```cpp
const CraiIndex index = CraiIndex::FromFile("input.cram.crai");
const std::vector<CraiEntry> chr1 = index.EntriesForReference(0);
const std::vector<CraiEntry> unmapped = index.EntriesForReference(-1);
const std::vector<CraiEntry> slices = index.Query(0, 1000, 2000);  // [beg, end), one row per slice
```

A `CraiEntry` holds `SequenceId`, `AlignmentStart`, `AlignmentSpan`, `ContainerOffset`,
`SliceOffset`, and `SliceSize`. `FromFile` throws `std::runtime_error` on parse or I/O
errors. There is no accessor for all rows.

---

## CigarOp

```cpp
#include <pbsamoa/core/CigarOp.hpp>
```

32-bit CIGAR operation matching BAM binary layout.

```cpp
CigarOp op{CigarOpType::M, 100};  // 100M
CigarOpType type = op.Type();     // CigarOpType::M
std::uint32_t length = op.Length();
std::uint32_t raw = op.RawValue();  // length << 4 | op, as in BAM
CigarOp fromRaw{raw};               // explicit, from the BAM value
```

A default `CigarOp` is `0M`. `CigarOp` supports `==`.

`CigarView` is `std::span<const CigarOp>`.

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
std::int32_t mismatches = CountMismatches(cigarOps);  // sum of X lengths
std::string str = CigarToString(cigarOps);         // "100M50I25M"
std::expected<std::vector<CigarOp>, std::string> ops = ParseCigar("100M50I25M");
if (!ops) {
    // ops.error() describes the problem
}
std::uint16_t bin = Reg2Bin(beg, end);
std::vector<std::uint16_t> bins = Reg2Bins(beg, end);  // all bins that overlap
CigarStats stats = ComputeCigarStats(cigarOps);       // Identity, Mismatches, ...
```

`ParseCigar("*")` gives an empty vector. `CigarToString` of an empty CIGAR gives `"*"`.
`WriteCigarTo(cigar, out)` appends the text to a string. `CigarOpToChar`, `CharToCigarOp`
(gives `std::optional`), `ConsumesQuery`, and `ConsumesReference` work on one operation type.
`CigarStats` has `NumAlignedBases`, `Span`, `Mismatches`, `Insertions`, `Deletions`,
`InsertionEvents`, `DeletionEvents`, `Identity`, and `IdentityGapComp`, as in pbcopper.

---

## Tags

```cpp
#include <pbsamoa/core/Tags.hpp>
```

### TagKey

Two-character key stored as `uint16_t`:

```cpp
constexpr TagKey nm{'N', 'M'};
std::string text = nm.ToString();  // "NM"
```

`First()`, `Second()`, and `Value()` give the parts. `RG_TAG` is `TagKey{'R', 'G'}`.

### TagValue

Variant: `char | int64_t | float | string | HexString | TagArray`.

Integer types (`c/C/s/S/i/I` in BAM) are stored uniformly as `int64_t`.
Serialization picks the smallest sufficient type.

`HexString` holds the text of an `H` tag in `value`. `TagArray` holds a `B` array:

```cpp
TagArray qualities{'C'};  // element type c, C, s, S, i, I, or f
qualities.AppendUInt8(30);
qualities.AppendUInt8(31);
std::uint32_t count = qualities.Count();
std::span<const std::byte> bytes = qualities.Data();  // little-endian elements
```

`ElementType()` and `ElementSize()` give the type code and its size. There is an `Append`
function for each type (`AppendInt8`, `AppendUInt8`, `AppendInt16`, `AppendUInt16`,
`AppendInt32`, `AppendUInt32`, `AppendFloat`). `Resize(count)` and `MutableData()` fill the
elements in place. The constructor does not check the type code: for an invalid one,
`ElementSize()` and `Count()` give 0.

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

`TagMap` keeps the entries sorted by key. It also has `Remove(key)` (gives `true` if the
key was there), `Contains(key)`, `Size()`, and `Empty()`. `Append(key, value)` does not
check for a duplicate key. Use it only for a key that is not in the map.

### Parsing and serialization

```cpp
TagMap tags = ParseTagsFromBam(auxBytes);
auto kv = ParseTagFromSam("NM:i:5");  // std::optional<std::pair<TagKey, TagValue>>
std::vector<std::byte> raw = SerializeTagsToBam(tags);
std::string text = SerializeTagToSam(TagKey{'N','M'}, TagValue{std::int64_t{5}});

std::string line;
SerializeRawTagsToSam(auxBytes, line);  // appends "\tNM:i:5"... without a TagMap
std::size_t entrySize = BamAuxEntrySize(auxBytes, 0);  // bytes of the first entry
```

`ParseTagsFromBam`, `SerializeRawTagsToSam`, and `BamAuxEntrySize` throw
`std::runtime_error` on malformed or unsupported tag data. `SerializedBamSize(tags)` and
`AppendTagsToBam(tags, dest)` write tags into a buffer that you allocate.

### Tag filters

```cpp
#include <pbsamoa/core/RawRecord.hpp>  // RawRecordView::ToOwned

DropTags drop{TagKey{'i','p'}, TagKey{'p','w'}};
KeepTags keep{TagKey{'R','G'}, TagKey{'N','M'}};

auto withoutKinetics = view.ToOwned(drop);  // copy everything except ip, pw
auto rgAndNm = view.ToOwned(keep);          // copy only RG, NM
```

`DropTags::ShouldDrop(key)` and `KeepTags::ShouldKeep(key)` test one key. Both hold a
`SortedTagKeySet` (binary search, `Contains(key)`).

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

VirtualOffset fromRaw{std::uint64_t{0x10000}};  // explicit, from the packed 64-bit value
std::uint64_t packed = fromRaw.Value();
```

A default `VirtualOffset` is 0. Supports `<=>` comparison. No arithmetic — this is by design.

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

The one-argument constructor gives block-only mode. The two-argument constructor
`BgzfReader{path, numWorkers}` gives BAM mode: it parses the BAM header at once and adds
`Header()`, `ReadRecord()`, and `ReadBatch(maxBytes, maxRecords)`. With `numWorkers > 0`,
it decompresses on a thread pool. These BAM-mode functions throw `std::runtime_error`
in block-only mode. Both modes have `Seek()`, `Tell()`, `HasEofMarker()`, and `GetMetrics()`.
The constructors throw `std::runtime_error` if the file does not open. The two-argument
constructor also throws it for an invalid or empty BAM header. `Tell()` gives the start of
the last block read in block-only mode, and the exact offset of the next record in BAM mode.

```cpp
BgzfReader bam{"file.bam", /*numWorkers=*/0};
while (std::optional<RawRecord> record = bam.ReadRecord()) {
}
```

Free functions work on single blocks: `ParseBgzfBlockHeader` (gives
`std::optional<BgzfBlockInfo>`), `DecompressBgzfBlock`, `CompressBgzfBlock(input, level)`,
and `IsBgzfEofMarker`. `BGZF_EOF_MARKER` holds the 28-byte EOF block.
`BGZF_MAX_UNCOMPRESSED_BLOCK` (65280) is the largest payload of one block.

---

## BgzfWriter

```cpp
#include <pbsamoa/core/Bgzf.hpp>
```

BGZF compression. Accumulates into blocks of up to `BGZF_MAX_UNCOMPRESSED_BLOCK` (65280)
bytes, compresses with libdeflate.

`BgzfWriterConfig` fields:

|        Field         | Default |                       Description                        |
| -------------------- | ------- | -------------------------------------------------------- |
| `CompressionLevel`   | `6`     | libdeflate level, `[1, 12]`                              |
| `BgzfWorkers`        | `4`     | Compression threads; must be > 0 when `Pool` is null     |
| `InputQueueCapacity` | `256`   | Input queue slots; must be > 0 when `Pool` is null       |
| `BlocksPerBatch`     | `32`    | Blocks in one compression task                           |
| `UseTempFile`        | `false` | Write to a temporary path and rename on `Close()`        |
| `Pool`               | `{}`    | `shared_ptr<BgzfWriterPool>`; shared mode when not null  |

The constructor throws `std::invalid_argument` for a value out of range.
`Write` also accepts a `std::vector<std::byte>&&` and does not copy it then. Both overloads
take an optional `PendingCallback` (within-block offset, record bytes, `rgId`, `zmw`) for the
index callback; `BamWriter` uses it.
`SetCallback(fn)` sets an index callback, and `EndVirtualOffset()` gives the end offset
after `Close()`. `GetMetrics()` gives `BgzfWriteMetrics`.

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

To write many destinations without creating threads per file, inject one pool:

```cpp
auto pool = std::make_shared<BgzfWriterPool>(8);  // 8 compressors + 1 ordered output thread
BgzfWriterConfig sharedConfig{.Pool = pool};
BgzfWriter first{"first.bgzf", sharedConfig};
BgzfWriter second{"second.bgzf", sharedConfig};
```

Shared mode packs on the caller. Calls to a particular writer must remain serial;
different writers may submit concurrently. `BgzfWorkers` and `InputQueueCapacity`
only configure the privately owned pipeline when `Pool` is null. The pool bounds
both work and result queues with its `queueMultiplier` constructor argument
(default 3, multiplied by worker count). Each writer retains its own partially
filled block batch. Writers retain shared ownership of the pool and close only
their own output; the last owner shuts down the pool. Pool metrics describe the
whole shared pool, and input-queue stall counters are zero in shared mode.

Results use the existing global FIFO consumer: a slow batch or output file can
delay other files. Callbacks execute on that consumer thread and must not submit
to or close writers on the same pool, which could wait for the consumer itself.

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
std::size_t bytes = precise.Value();
```

The `std::size_t` constructor is `explicit`. The literals are in
`PacBio::Samoa::Literals`.

---

## Error Policies

```cpp
#include <pbsamoa/io/BamRawReader.hpp>
```

|    Policy     |                                  Behavior                                   |
| ------------- | --------------------------------------------------------------------------- |
| `ThrowPolicy` | `OnCorruptRecord(message)` throws `std::runtime_error` with the message     |
| `SkipPolicy`  | `OnCorruptRecord(message)` counts the call. `SkippedCount()` gives the count |

No reader in this release takes a policy. The readers throw an exception on corrupt
data: usually `std::runtime_error`, and `std::invalid_argument` from the `RawRecordView`
field checks.

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

|             Method              |      Return type      |                Description                 |
| ------------------------------- | --------------------- | ------------------------------------------ |
| `Find(int32_t zmw)`             | `vector<int64_t>`     | Offsets for hole number (all read groups)  |
| `Find(ZmwIdentity)`             | `vector<int64_t>`     | Offsets for exact (rgId, zmw)              |
| `Find(span<const ZmwIdentity>)` | `vector<int64_t>`     | Batch query                                |
| `FirstOffset(int32_t zmw)`      | `int64_t`             | First offset for hole number               |
| `FirstOffset(ZmwIdentity)`      | `int64_t`             | First offset for exact (rgId, zmw)         |
| `UniqueZmws()`                  | `vector<ZmwIdentity>` | All unique (rgId, zmw) in file order       |
| `ContiguousChunk(i, n)`         | `ZmwChunk`            | Seek offset and sizes for chunk `i` of `n` |
| `NumRecords()`                  | `uint64_t`            | Total records in the index                 |
| `NumZmws()`                     | `uint64_t`            | Unique ZMW count                           |
| `Identity(row)`                 | `ZmwIdentity`         | (rgId, zmw) of record `row` in file order  |
| `QueryLengths()`                | `span<const int32_t>` | `qEnd - qStart` per record; `.pbi` only    |
| `HasBarcodeData()`              | `bool`                | The index has barcode fields (ZMI 1.1)     |
| `FindBarcodes(BarcodePair)`     | `span<const int64_t>` | Offsets of one barcode pair, file order    |
| `BarcodePairs()`                | `vector<BarcodePair>` | Distinct barcode pairs, sorted             |

Offsets are raw virtual offsets (`std::int64_t`). `ZmwIdentity` holds `rgId` and `zmw`.
`ZmwChunk` holds `FirstOffset`, `NumRecords`, and `NumZmws`. `ContiguousChunk` uses a 1-based
chunk number and throws `std::invalid_argument` for an invalid chunk.
`FirstOffset` throws `std::runtime_error` if the ZMW is not in the index. `ZmwIndex` is
move-only.

### Barcode lookup

ZMI 1.1 records the ordered forward/reverse barcode IDs from `bc:B:S` or
`bc:B:s`. `ZmiBamWriter`, `zmi-build`, and `zmi-index` capture them automatically.
`BarcodePair` has `Forward` and `Reverse`, both `-1` by default.
An absent tag is `BarcodePair{-1, -1}`; signed -1 components remain unassigned.
Unsigned IDs, including 65535, are preserved. Malformed barcode tags are rejected.

```cpp
#include <bit>
#include <pbsamoa/io/BamRawReader.hpp>  // reader, VirtualOffset

BamRawReader reader{"movie.bam"};
const ZmwIndex index = ZmwIndex::Open("movie.bam");
for (const BarcodePair pair : index.BarcodePairs()) {
    for (const std::int64_t offset : index.FindBarcodes(pair)) {
        reader.Seek(VirtualOffset{std::bit_cast<std::uint64_t>(offset)});
        // Read this record; offsets retain source file order.
    }
}
```

`HasBarcodeData()` distinguishes an index without barcode fields from a new
index whose records are all unassigned. `FindBarcodes()` returns a borrowed
`span<const int64_t>`; `BarcodePairs()` returns distinct pairs in lexicographic
order. Pairs are directional: `{a, b}` and `{b, a}` remain distinct. Queries
are safe concurrently on an unchanged index; moving or destroying it invalidates
borrowed results. Loading groups offsets once, using one additional int64 offset
per record plus per-pair map/vector storage, so each lookup avoids scanning the
whole index.

The BGZF-compressed format retains the 64-byte header and `ZMI\1` magic. Version
at byte 4 is now `0x010100`; entry size at byte 8 is 24. Each little-endian entry
retains `rgId:int32`, `zmw:int32`, and `virtualOffset:int64` at offsets 0, 4, and 8,
and appends forward/reverse `int32` IDs at offsets 16 and 20. IDs range from -1
through 65535. Older readers that honor entry size can skip the extra fields.
Legacy 16-byte ZMI entries remain readable for ZMW queries. Barcode queries throw
when fields are unavailable, including through `FromPbi()`, which still reads
only PBI's ZMW columns; rebuild the ZMI to enable barcode lookup.

---

## BamZmwReader

```cpp
#include <pbsamoa/io/BamZmwReader.hpp>
```

Groups consecutive records by ZMW identity. Wraps a `BamRecordReader` and
uses a background producer thread to prefetch complete ZMW groups into a
bounded queue. Assumes the source yields records grouped by ZMW
(standard PacBio BAM ordering). The `rgId` of the identity comes from the base read-group
ID in the `RG` tag. It is `0` when the record has no such tag. The `ZmwGroup` struct in
the header is internal; no public function takes or returns it.

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

The constructor also accepts a `BamCollection` or a `std::vector<BamFile>`.
`NumZmws()` gives the ZMW count of the chunk, or of the file from its ZMW index. It gives
-1 when no count is available.

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
- Underlying reader: `Reader` (`ReaderMetrics`)

---

## BamSort

```cpp
#include <pbsamoa/io/BamSort.hpp>
```

External (RAM-bounded) sort of a BAM file. The record content matches `samtools sort`.
The BGZF bytes can be different.

```cpp
SortStats stats = SortBam("in.bam", "sorted.bam", SortConfig{
    .Order = SortOrder::COORDINATE,
    .MaxMemory = ByteLimit{std::size_t{2} * 1024 * 1024 * 1024},
    .NumThreads = 8,
});
// stats.NumRecords, stats.NumRuns (0: sorted in memory)

SortConfig byTag{.Order = SortOrder::TAG, .Tag = {'H', 'P'}};
```

`SortOrder` values:
- `COORDINATE`: by reference, position, and strand; unmapped last (`SO:coordinate`).
- `QUERY_NAME`: natural read-name order (`SO:queryname`).
- `TAG`: by a 2-char tag, then by coordinate (`SO:unknown`).

`SortConfig` fields:

|         Field         | Default      |                               Description                               |
| --------------------- | ------------ | ----------------------------------------------------------------------- |
| `Order`               | `COORDINATE` | Sort order                                                              |
| `Tag`                 | `{}`         | `std::array<char, 2>`; required for `TAG`                               |
| `MaxMemory`           | `768 MiB`    | Peak record-run buffer (one run at a time)                              |
| `NumThreads`          | `0`          | BGZF threads; `0` = `min(hardware_concurrency, 8)`                      |
| `CompressionLevel`    | `6`          | Output level `[1, 12]`; temporary runs use level 1                      |
| `TempDir`             | `{}`         | Directory of the run files; default: the directory of the output, or the current directory for stream output |
| `Minimise`            | `false`      | Cluster unmapped reads by minimiser (`samtools sort -M -R`); coordinate only |
| `CommandLine`         | `{}`         | `CL` value of the added `@PG` line                                      |
| `DirectToOutput`      | `false`      | Write to the output directly, not through a temporary file and rename   |
| `MergeFanIn`          | `64`         | Run readers per merge pass; at least 2                                  |
| `OutputIndexCallback` | `{}`         | `IndexCallback` for each record of the final output                     |

Set `DirectToOutput` for `/dev/stdout`, `/dev/fd/*`, or a FIFO. The sort also detects an
existing output that is not a regular file. `OutputIndexCallback` runs on the IO thread of
the output writer, only for the final file. All calls finish before `SortBam` returns.
`SortBam` throws `std::runtime_error` on I/O errors, a missing input, or `Minimise` with an
order other than `COORDINATE`. It throws `std::invalid_argument` if `MergeFanIn < 2` or the
record source is empty.

### Selected records

`SortBam(header, nextRecord, output, config)` accepts a
`std::function<std::optional<RawRecord>()>` source and the destination header.
Return nullopt at end of input. This permits ZMI barcode extraction directly
into the existing external sort without creating a per-barcode input BAM.
The function releases its source after ingestion, before merging runs. The
caller must keep borrowed offsets and other externally owned source data alive.
Source exceptions propagate and leave an existing regular output unchanged.

`SortConfig::MergeFanIn` defaults to 64 and must be at least 2. Concurrent callers
must divide their available descriptor budget across jobs, allowing additional
headroom for each output and any source handles they retain. `MaxMemory` bounds
each sort's record run buffer; it does not include external indexes or codec state.

---

## BamMerge

```cpp
#include <pbsamoa/io/BamMerge.hpp>
```

Merges several BAM files into one. All inputs must have the same `@SQ` list. The merged
header joins their `@RG`, `@PG`, and `@CO` lines. The output is written atomically.

```cpp
const std::vector<std::filesystem::path> inputs{"a.bam", "b.bam"};
MergeStats stats = MergeBam(inputs, "merged.bam", MergeConfig{
    .NumThreads = 8,
    .BaiOutput = std::filesystem::path{"merged.bam.bai"},
});
// stats.NumRecords (-1 for concat), stats.NumInputs, stats.Passthrough, stats.Runtime
```

Modes:
- Sorted merge: a k-way merge of inputs sorted in the same order. Equal keys keep the input
  order. If coordinate-sorted inputs do not overlap, the merge orders them by their first
  record and copies the BGZF blocks (`MergeStats::Passthrough`).
- Concat (`Concat = true`, or all inputs unsorted): copies the BGZF blocks and marks the
  output `SO:unsorted`.

`MergeConfig` fields:

|         Field         | Default   |                               Description                                |
| --------------------- | --------- | ------------------------------------------------------------------------ |
| `Order`               | `nullopt` | `optional<SortOrder>`; `nullopt` detects the order from `@HD SO`         |
| `Tag`                 | `{}`      | `std::array<char, 2>`; required for `TAG`                                |
| `NumThreads`          | `0`       | Size of both worker pools; `0` = `min(hardware_concurrency, 8)`          |
| `DecodeThreads`       | `0`       | Input decompression threads; `0` = `NumThreads`                          |
| `CompressThreads`     | `0`       | Output compression threads; `0` = `NumThreads`                           |
| `CompressionLevel`    | `6`       | Output level `[1, 12]`                                                   |
| `ReadAheadMemory`     | `768 MiB` | Read-ahead budget for all inputs (sorted merge only)                     |
| `BatchBytes`          | `0`       | Record bytes per input batch; `0` = 256 KiB                              |
| `WriterQueueCapacity` | `0`       | Output writer queue slots; `0` = 256                                     |
| `CommandLine`         | `{}`      | `CL` value of the added `@PG` line                                       |
| `Concat`              | `false`   | Concatenate the inputs and ignore their sort order                       |
| `BaiOutput`           | `{}`      | `optional<path>`; build a BAI during a coordinate heap merge             |

`ReadAheadMemory` is not a limit for the process. Peak resident memory is about 1.4 times
the budget, plus one batch for each input and the writer queue. Under a hard memory limit,
set the limit to at least 2 times the budget.

`MergeBam` throws `std::runtime_error` on I/O errors, a missing input, mixed sort orders,
or different reference lists. `BaiOutput` together with concat, a passthrough merge, or a
queryname or tag order is an error. `MergeRuntimeStats` (`stats.Runtime`) holds wall times,
wait times, and byte counts for a bottleneck report.

### MergeReadAhead

```cpp
#include <pbsamoa/io/MergeReadAhead.hpp>
```

The parallel read-ahead of the sorted merge. Each input has a producer thread. A shared
thread pool decompresses the blocks. One budget limits the decompressed bytes in flight.

```cpp
using namespace PacBio::Samoa::Literals;

MergeReadAhead readAhead{inputs, /*decodeWorkers=*/4, 256_MiB, /*batchBytes=*/0};
for (std::size_t i{0}; i < readAhead.NumSources(); ++i) {
    while (std::optional<std::span<const std::byte>> bytes = readAhead.Next(i)) {
        // record bytes after block_size; valid until the next Next(i)
    }
}
MergeReadAheadStats stats = readAhead.Stats();  // after all sources reach their end
```

`decodeWorkers` is at least 1. `batchBytes` is the record bytes of one batch per input
(`0` = 256 KiB); each input holds one batch outside the budget.
`Next(i)` gives the records of input `i` in file order. It throws the exception of the
producer on a read error. The constructor throws `std::runtime_error` if an input does
not open. `MergeReadAhead` cannot be copied or moved.

## ZmiBamWriter

```cpp
#include <pbsamoa/io/ZmiBamWriter.hpp>
```

Writes BAM and `.zmi` index simultaneously. The index goes to `<output>.zmi`.

```cpp
ZmiBamWriter writer{"output.bam", header};
writer.Write(record);       // updates both BAM and ZMI
writer.WriteBatch(batch);   // batch write to both BAM and ZMI
writer.Close();

ZmiBamWriter tuned{"tuned.bam", header, BamWriterConfig{.UseTempFile = true}};
```

`Write` also accepts a `RawRecordView` or a `RawRecord`. For its ZMI entry the writer takes
the base read-group ID from `RG`, the hole number from the read name (`movie/<zmw>/...`, 0
if the name does not match), and the barcodes from `bc`.

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
writer.AddRecord(rgId, zmw, virtualOffset);                     // unassigned barcodes
writer.AddRecord(rgId, zmw, virtualOffset, BarcodePair{12, 12});  // forward, reverse
writer.Close();
```

Add the entries in BAM write order. `AddRecord` throws `std::invalid_argument` for a
barcode ID outside `[-1, 65535]`.

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

std::vector<std::byte> dest(2);
PackSequenceInto("ACGT", dest.data());  // dest needs (size + 1) / 2 bytes
std::string out;
WriteSequenceTo(packed, 4, out);        // appends "ACGT"
```

The tables `BASE_TO_CODE` (ASCII to 4-bit code; unknown characters give 15, `N`),
`CODE_TO_BASE`, and `PACKED_TO_BASES` (one packed byte to two bases) are public.

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
| `SubstringClipStrategy`        | `dq`, `iq`, `mq`, `sq`, `dt`, `st`, `ip`, `pw`, `fi`, `fp`, `sm`, `sx`       | Simple substring           |
| `ReverseSubstringClipStrategy` | `ri`, `rp`                                                                   | Mirrored offset substring  |
| `PulseClipStrategy`            | `pc`, `pt`, `pq`, `pv`, `pg`, `pa`, `pm`, `ps`, `pi`, `pd`, `px`, `pe`, `sf` | Pulse-space via `pc` tag   |
| `BasemodClipStrategy`          | `MM`, `ML`                                                                   | SAMv1.6 base modifications |
| `PileupClipStrategy`           | `sa`                                                                         | RLE pileup coverage        |

### Usage

```cpp
auto clipper = TagClipper::PacBioDefault();  // all standard strategies
clipper.ClipTags(tags, clipOffset, clipLength, seqLength, sequence, /*isReverse=*/false);

// Add a strategy for more tags; the clipper keeps a reference to the strategy
static const SubstringClipStrategy substring;
clipper.Register({TagKey{'x', 'q'}}, substring);
```

`clipOffset` and `clipLength` give the kept query window. `sequence` is the sequence
before the clip. Without it, the clipper removes `MM` and `ML`. Set `isReverse` for a
reverse-strand record. `ClipTags` does not change tags that have no strategy, except that
it sets an existing `MN` to `clipLength`. A
strategy derives from `TagClipStrategy` and implements `Clip(value, clipOffset,
clipLength, seqLength, ctx)`. It returns `false` to remove a malformed tag. `ClipContext`
gives a strategy access to the other tags.
The strategy must stay valid while the clipper uses it.

---

## Metrics

```cpp
#include <pbsamoa/core/Metrics.hpp>
```

Runtime performance snapshots for all reader/writer pipelines. Captured
via relaxed atomics — diff two snapshots for per-second rates.

### Types

|        Type        |                      Source                      |
| ------------------ | ------------------------------------------------ |
| `BgzfMetrics`      | `BamRawReader::GetMetrics()`, `BgzfReader`       |
| `DecodeMetrics`    | `ReaderMetrics::Decode`                          |
| `ReaderMetrics`    | `BamRecordReader::GetMetrics()`                  |
| `ZmwReaderMetrics` | `BamZmwReader::GetMetrics()`                     |
| `BgzfWriteMetrics` | `WriterMetrics::Bgzf`, `BgzfWriter::GetMetrics()` |
| `WriterMetrics`    | `BamWriter::GetMetrics()`                        |
| `PoolMetrics`      | The `Pool` member of the other types             |

### ReaderMetrics fields

`Bgzf` (`BgzfMetrics`), `Decode` (`DecodeMetrics`), `TotalRecordsRead`, `ParallelBgzf`,
`ParallelDecode`.

### DecodeMetrics fields

Pool: `Pool`. Throughput: `BatchesDecoded`, `RecordsDecoded`.
Queue: `RecordsProduced`, `RecordsConsumed`. Stalls: `ProducerStalls`, `ConsumerStalls`.
Timing (ns): `DecodeNs`, `BatchReadNs`.

### WriterMetrics fields

`Bgzf` (`BgzfWriteMetrics`) and `TotalRecordsWritten`. `BgzfWriteMetrics` has
stalls: `CallerStalls`, `PackerStalls`, `WriterStalls`, `CallerStallNs`; throughput:
`BytesCompressed`, `BlocksWritten`; timing (ns): `CompressNs`, `IoWriteNs`, `CallbackNs`;
pool: `Pool`.

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
Underlying reader: `Reader` (`ReaderMetrics`).

---

## BamCollection and BamFile

```cpp
#include <pbsamoa/io/BamCollection.hpp>
#include <pbsamoa/io/BamFile.hpp>
```

`BamFile` wraps one BAM path. The constructor reads the header.

```cpp
const BamFile file{"movie.bam"};
const SamHeader& header = file.Header();
const bool complete = file.HasEOF();             // ends with the BGZF EOF marker
const std::filesystem::path& path = file.Filename();
const std::filesystem::path bai = file.StandardIndexFilename();  // "movie.bam.bai"
file.EnsureStandardIndexExists();                // builds the .bai if it is missing or old
const std::uint32_t chr1Length = file.ReferenceLength("chr1");  // 0 if unknown
```

`CreateStandardIndex()` always builds the `.bai`. `StandardIndexExists()` and
`StandardIndexIsNewer()` test the sidecar. `ReferenceName(id)` gives an empty string for an
unknown ID.

`BamCollection` is an ordered set of BAM files with one merged header. The readers accept
it in place of one path.

```cpp
#include <pbsamoa/io/BamRecordReader.hpp>

const BamCollection collection{std::vector<std::filesystem::path>{"a.bam", "b.bam"}};
const SamHeader& merged = collection.Header();
for (const BamFile& member : collection.Files()) {  // collection.Size() files
}
BamRecordReader reader{collection};
```

The constructor also takes one path or one `BamFile`. The files must have the same `@SQ`
names and lengths, in the same order. The merged header
joins the `@RG`, `@PG`, and `@CO` lines. With more than one file, its sort order is
`unknown`. The constructor throws `std::runtime_error` for conflicting headers and
`std::invalid_argument` for an empty list.

---

## GenomicInterval

```cpp
#include <pbsamoa/core/GenomicInterval.hpp>
```

A half-open interval `[start, stop)` on a named reference. `BamRecordReader::Query`
uses it.

```cpp
GenomicInterval interval{"chr1", 1000, 2000};
interval.Start(1500).Stop(2500);
const bool empty = interval.Empty();  // start == stop
```

The constructor and the setters throw `std::invalid_argument` if `start < 0` or
`stop < start`.

---

## CramReference

```cpp
#include <pbsamoa/cram/CramReference.hpp>
```

An abstract reference for `CramWriterConfig::Reference` and `CramReaderConfig::Reference`.
Implement it to give the bases in windows, for example from a packed file. Several threads
call `Bases` at the same time.

```cpp
#include <pbsamoa/io/CramWriter.hpp>  // CramWriterConfig

class InMemoryReference final : public CramReference
{
public:
    explicit InMemoryReference(std::vector<std::string> sequences)
        : sequences_{std::move(sequences)}
    {}

    std::int32_t NumSequences() const override { return std::ssize(sequences_); }

    std::int64_t Length(const std::int32_t refId) const override
    {
        return std::ssize(sequences_[refId]);
    }

    std::string_view Bases(const std::int32_t refId, const std::int64_t start,
                           const std::int64_t length, std::string& /*buffer*/) const override
    {
        return std::string_view{sequences_[refId]}.substr(start, length);
    }

private:
    std::vector<std::string> sequences_;
};

auto reference = std::make_shared<const InMemoryReference>(std::vector<std::string>{"ACGT"});
CramWriterConfig config{.Reference = reference};
```

Give one sequence per `@SQ` line, in header order. `Bases` can return a view of `buffer`,
which it can change. The case of the bases does not matter.

## FastaFile

```cpp
#include <pbsamoa/io/FastaReference.hpp>
```

A FASTA file as a `CramReference`, without the text of the whole genome in memory. The
constructor maps the file and takes the offset of each sequence from `<fasta>.fai`, or, without
one, from one scan of the file (a few seconds for GRCh38). A sequence whose lines are uneven
(different lengths, blank lines, spaces), which an index cannot describe, is read into memory. A
gzip-compressed or empty file has no sequences. The constructor throws `std::runtime_error` when
the file cannot be opened or mapped, or when the `.fai` does not fit the file.

```cpp
const FastaFile fasta{"ref.fa"};
fasta.NumSequences();  // sequences in the file
CramWriterConfig config{.Reference = fasta.Reference(header, MissingFastaSequence::THROW)};
```

`Reference(header, missing)` gives the sequences of the `@SQ` lines, in header order, matched by
name (the first word of the FASTA header line; the first of two equal names wins). An `@SQ`
sequence that is not in the file has length 0 (`MissingFastaSequence::EMPTY`, the default) or
throws `std::runtime_error` (`THROW`). The reference keeps the mapping alive after the
`FastaFile` is gone. A window on one line of the file is a view of the mapping; others are copied
to `buffer`.

---

## CigarClipping

```cpp
#include <pbsamoa/core/CigarClipping.hpp>
```

The CIGAR part of `BamRecord::Clip`. Both functions return a `ClipResult`: the new
`cigar`, the query offset `clipOffset`, the kept query length `clipLength`, and the new
reference position `newPos`.

```cpp
ClipResult byQuery = ClipCigarToQuery(cigar, /*queryStart=*/100, /*queryEnd=*/900,
                                      /*origQueryStart=*/0, /*origQueryEnd=*/1000,
                                      /*refPos=*/5000, /*isReverse=*/false);
ClipResult byRef = ClipCigarToReference(cigar, /*refStart=*/5100, /*refEnd=*/5800,
                                        /*origRefPos=*/5000, /*isReverse=*/false,
                                        /*exciseFlankingInserts=*/false);
```

---

## LosslessClipping

```cpp
#include <pbsamoa/core/LosslessClipping.hpp>
```

Clips a record to a query range and keeps the removed flanks, so that you can undo the clip.

```cpp
ClipToQueryLossless(record, /*clipLeft=*/24, /*clipRight=*/1976);
// ...
const bool restored = RestoreFromLossless(record);  // false if the record has no ls tag
```

`ClipToQueryLossless` stores the removed sequence, qualities, kinetics, base modifications,
and pileup data in an `ls` tag. The `ls` format is the format of lima, so `lima-undo` can
also restore the read. A record that already has `ls` keeps the old value under a `nested`
key. The clip does nothing if the range is empty, out of bounds, or the full sequence.
It throws `std::runtime_error` if the record has raw-subread pulse tags (`pc`, `pt`, ...)
or legacy QV tags (`dq`, `iq`, ...). Use it for unmapped reads: the undo does not restore
CIGAR and POS.

`RestoreFromLossless` rebuilds the record from the full `ls` chain. It then removes the
demultiplexing tags (`bc`, `bq`, `bx`, `bl`, `bt`, `ql`, `qt`, `ls`) and `qs`/`qe`. It also
removes the barcode bits in `cx` and the barcode suffix of `RG`.

---

## Basemods

```cpp
#include <pbsamoa/core/Basemods.hpp>
```

Parses and writes `MM` base-modification strings (SAMtags).

```cpp
std::vector<BasemodRecord> mods = ParseBasemodString("C+m,1,3,0;A+a,2;");
// mods[0].Prefix == "C+m", mods[0].Skips == {1, 3, 0}
std::string mm = WriteBasemodString(mods);
std::size_t stride = ModCodeCount("C+mh");  // 2 ML values per site
BasemodClipWindow window = ClipBasemodRecord(mods[0], "ACCCCCCCGT", /*clipOffset=*/2,
                                             /*clipLength=*/4, /*reverse=*/false);
```

`ClipBasemodRecord` splits the sites of one record into three groups. `FrontRemoved` counts
the sites before the window. `Retained` counts the sites in it, and `RetainedSkips` holds
their new skips. The other sites are after the window. `PrefixLost` is the count of canonical bases between the last front-removed site and the
clip start; undo adds it back to the first kept skip. For a reverse-strand record, set `reverse`, because `MM` uses the original
read orientation.

---

## Endian

```cpp
#include <pbsamoa/core/Endian.hpp>
```

`ReadLE<T>(p)`, `ReadU16LE`, `ReadU32LE`, `ReadI32LE`, `ReadU64LE`, and `ReadI64LE` read
little-endian integers from unaligned bytes. The header requires a little-endian platform
and a 64-bit `std::size_t` (`static_assert`).

---

## Library information

```cpp
#include <pbsamoa/PbSamoaLibraryInfo.hpp>
```

```cpp
const LibraryInfo info = GetLibraryInfo();  // Name, Release, GitSha1 (string_view)
const std::string version = LibraryFormattedVersion();  // "<release> (commit <sha1>)"
```
