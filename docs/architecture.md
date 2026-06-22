# Architecture

pbsamoa is organized in four layers. Each layer depends only on the layers
below it.

```
┌────────────────────────────────────────────────────────┐
│  API                                                   │
│  BamRawReader  BamRecordReader  SamReader  CramReader  │
│  BamWriter  SamWriter  CramWriter  BaiIndex  CraiIndex │
│  BamZmwReader  ZmiBamWriter                            │
├────────────────────────────────────────────────────────┤
│  Record                                                │
│  RawRecord  BamRecord  RawRecordBatch                  │
│  ZmwGroup                                              │
│  CigarOp  SequenceView  TagMap  SamHeader              │
├────────────────────────────────────────────────────────┤
│  Compression                                           │
│  BgzfReader  BgzfWriter  VirtualOffset  CramCodec      │
│  CramCompression                                       │
└────────────────────────────────────────────────────────┘
```

## Class Dependencies

Arrows show composition ("owns") or required dependencies ("uses"). Dotted
arrows show decode-on-demand relationships. The only inheritance hierarchy
is `TagClipStrategy` and its subclasses (see `TagClipping.hpp`).

```
┌─ API ──────────────────────────────────────────────────────────────────┐
│                                                                        │
│  BamRawReader ──┬──► BgzfReader (sync or parallel via numWorkers)      │
│                 ├──► SamHeader                                         │
│                 └──► RawRecordBatch                                    │
│                       └──► RecordData(i) → span<const byte>            │
│                                                                        │
│  BamRecordReader ┬──► BamRawReader                                     │
│                  ├──► SPSCQueue<BamRecord> (pre-decode pipeline)       │
│                  ├──► std::jthread (background producer)               │
│                  └──► ThreadPool (optional, for parallel ToOwned)      │
│                                                                        │
│  BamZmwReader ───► BamRecordReader (synchronous ZMW grouping)          │
│                                                                        │
│  BamWriter ──┬──► BgzfWriter                                           │
│              └──► SamHeader                                            │
│                                                                        │
│  SamReader ──┬──► SamHeader                                            │
│              └──► BamRecord (produces per record)                      │
│                                                                        │
│  SamWriter ─────► SamHeader                                            │
│                                                                        │
│  CramReader ──┬──► SamHeader                                           │
│               ├──► CramCodec + CramCompression                         │
│               └──► CraiIndex (for region query)                        │
│                                                                        │
│  CramWriter ──┬──► SamHeader                                           │
│               ├──► CramCodec + CramCompression                         │
│               └──► optional CRAI writer                                │
│                                                                        │
│  BaiIndex ──────► Chunk ──► VirtualOffset                              │
│  CraiIndex ─────► CraiEntry                                            │
│                                                                        │
│  ZmiBamWriter ┬──► BamWriter                                           │
│               └──► ZmiWriter ──► BgzfWriter                            │
│                                                                        │
├─ Record ───────────────────────────────────────────────────────────────┤
│                                                                        │
│  RawRecord ─────···► SequenceView, TagMap (decode on demand)           │
│             ────────► CigarOp (eagerly copied on construction)         │
│                                                                        │
│  BamRecord ──┬──► std::vector<CigarOp>                                 │
│              └──► TagMap ──► TagKey + TagValue                         │
│                                                                        │
│  RawRecordBatch ────► buffer + RecordExtent vector                     │
│                                                                        │
│  ZmwGroup ───┬──► ZmwIdentity                                          │
│              └──► std::vector<BamRecord>                               │
│                                                                        │
│  SamHeader ──┬──► ReferenceSequence (vector)                           │
│              ├──► ReadGroup (vector)                                   │
│              └──► ProgramRecord (vector)                               │
│                                                                        │
├─ Compression ──────────────────────────────────────────────────────────┤
│                                                                        │
│  BgzfReader ─────► VirtualOffset (current position)                    │
│                                                                        │
│  BgzfWriter ─────► parallel compression pipeline + callback dispatch   │
│                                                                        │
└────────────────────────────────────────────────────────────────────────┘
```

Key relationships:

- **BamRawReader** owns a `BgzfReader` (synchronous when `BgzfWorkers == 0`,
  parallel pipeline when `BgzfWorkers > 0`) for decompression, and a
  `SamHeader` parsed from the BAM header block. It produces `RawRecordBatch`
  objects that own their backing buffers and provide indexed access to raw
  record spans via `RecordData(i)`. Supports chunked reading via
  `ChunkNum`/`TotalChunks` in config and ZMW whitelist filtering via
  `Whitelist()`.

- **BamRecordReader** wraps a `BamRawReader` and pre-decodes views into
  owned `BamRecord` objects in a background `std::jthread`. Records are
  pushed through an SPSC queue for consumption. Optional `ThreadPool`
  parallelizes the `ToOwned()` decode step across batch views.

- **BamZmwReader** wraps a `BamRecordReader` and groups consecutive records
  by ZMW hole number, parsed from the read name (`movie/zmw/…`). Grouping is
  synchronous — no background threads. `CurrentZmw()` returns the current
  group's `ZmwIdentity` with `rgId` always set to 0.

- **BamWriter** owns a `BgzfWriter` for compression and a `SamHeader` for
  the output file header. Accepts `BamRecord`, `RawRecord`, raw
  `span<const byte>`, and `RawRecordBatch` (via `WriteBatch()`).

- **CramReader** owns CRAM container/slice decode state and produces owned
  `BamRecord` objects. Container payload parsing is centralized through
  `ParseContainer(...)`, then slices are decompressed/decoded in order.
  It supports optional reference-based decoding and region queries via
  `CraiIndex`, plus `RawRecord` accessors that serialize decoded records into
  BAM layout on demand.

- **CramWriter** uses a two-level write buffer:
  `pendingRecords` (records for the current slice) and `pendingSlices`
  (completed slices for the current container). `RecordsPerSlice` controls
  slice boundaries; `SlicesPerContainer` controls container flush boundaries.
  Slice encode/compress can run in parallel, while container writes are kept
  ordered with an async write future. It can emit a samtools-compatible CRAI
  sidecar with per-slice offsets.

- **ZmiBamWriter** composes a `BamWriter` and a `ZmiWriter` so every record
  written is simultaneously indexed.

- **SortBam / MergeBam** (free functions in `io/BamSort.hpp`, `io/BamMerge.hpp`)
  combine files rather than streaming a single one. `SortBam` is a RAM-bounded
  external merge sort (spill sorted runs, then the k-way `MergeRuns` core, which
  caps simultaneously-open run files and falls back to a multi-pass merge when the
  run count would exceed the process open-file limit). `MergeBam` detects the mode
  from the inputs' `@HD SO`: a sorted merge (coordinate/queryname; deterministic via
  source-index tie-break) or, for unsorted inputs / `--concat`, byte-level
  concatenation through BGZF block passthrough (compressed blocks copied
  input→output via `CompressBgzfBlock` only for the header and the one block
  straddling each input's header/record boundary). A coordinate merge whose inputs
  form a strictly disjoint chain is detected by `ProbeDisjointChain` (reads each
  input's first record for its minimum key, then verifies — in min-sorted order —
  that each input's records are non-decreasing and strictly precede the next input's
  minimum, aborting early on the first overlapping or out-of-order record). When the
  chain holds, the merge reuses that same passthrough copy in min-sorted order,
  eliminating recompression; otherwise it falls back to the heap merge below. The
  strict `max[i] < min[i+1]` test guarantees no cross-file equal keys, so the
  passthrough output is record-identical to the heap merge.

- **MergeReadAhead** (`io/MergeReadAhead.hpp`) backs the sorted merge. Each input
  gets a producer thread that reads raw BGZF blocks and decompresses them on a
  shared `ThreadPool`, framing records into per-input batches handed to the merge
  heap as zero-copy views (the decompressed buffer is moved into the batch, not
  copied per record); the heap writes those spans straight to the writer, which
  compresses in parallel. A single `ReadAheadMemory` budget caps the summed
  in-flight decompressed bytes across all inputs, with condition-variable
  backpressure (no busy-spin) and a per-source head-batch guarantee so a sub-batch
  budget cannot stall the merge. Per-source file order is preserved, so the merged
  output is independent of thread count, decode timing, and budget. Avoiding the
  per-record allocation is what lets input decompression actually overlap the
  merge — on a 575 MB four-way coordinate merge it runs ~1.7x faster than
  `samtools merge -@8`.

- **RawRecord** owns a copy of raw BAM bytes. `CigarOp` values are eagerly
  copied into an aligned buffer on construction (BAM does not guarantee
  alignment). `SequenceView` and `TagMap` are decoded on demand. All decode
  logic is inline in the header for maximum performance.

## Compression

BAM files use BGZF (Blocked GZip Format) — a series of independently
compressed gzip blocks, each at most 64 KiB decompressed. This enables
random access via virtual offsets.

CRAM uses container/slice structures and block-level codecs selected from
CRAM v3.x methods (`raw`, `gzip`, `bzip2`, `lzma`, `rANS`, `arith`,
`fqzcomp`, `name tokeniser`) via `CramCodec` + `CramCompression`.

### BgzfReader

Reads BGZF blocks sequentially, decompressing each with libdeflate.
Supports `Seek()` to any `VirtualOffset`.

### BgzfWriter

Accumulates data into 64 KiB blocks and compresses each with libdeflate in a
parallel writer pipeline. Compression level 1–12 (default 6). Appends the
standard 28-byte EOF marker on close.

When a callback is configured (via `BamWriter`), callback dispatch happens on
the BGZF IO writer thread after block file offsets are known, not on the caller
thread. Callbacks can be deferred until `Close()`.

### VirtualOffset

A 64-bit value encoding position within a BGZF file: upper 48 bits are the
compressed block offset, lower 16 bits are the offset within the
decompressed block. Supports comparison via `<=>` but not arithmetic — this
is intentional per the BAM spec.

## Record

Two representations for alignment records, chosen by use case:

### RawRecord — fast, read-only

An owning type that copies raw BAM bytes into a `vector<byte>`. CIGAR ops
are eagerly copied into an aligned buffer on construction (BAM does not
guarantee 4-byte alignment). Other accessors decode fields on demand
directly from the binary layout. All decode logic is inline in the header
for maximum inlining.

Use raw records when you're reading records and don't need to modify them —
this is the common case.

### BamRecord — owned, mutable

Stores fields as structured C++ types (strings, vectors). Mutations are
plain field assignments with a fluent interface:

```cpp
record.Name("read1").MapQ(60).Pos(1000);
```

Serialization to BAM binary happens once, in the writer, via
`SerializeToBam()`.

### Conversion between views and records

`RawRecord::ToOwned()` produces a `BamRecord`. Optional tag filters
(`DropTags`, `KeepTags`) let you skip large tags during conversion. This
matters for PacBio BAM files where kinetics arrays can be megabytes per
record.

Both `BamWriter` and `SamWriter` accept either type (plus `WriteBatch()`
for `RawRecordBatch`). `BamWriter` also accepts raw `span<const byte>`.
Writing a view to `BamWriter` is zero-copy; writing a record serializes it.

### RawRecordBatch

A batch owns the decompressed buffer and provides indexed access to raw
record spans via `RecordData(i)`. Batches are sized by memory budget
(`ByteLimit`), not record count. The default budget is 256 MiB.

```cpp
using namespace PacBio::Samoa::Literals;
while (auto batch = reader.ReadBatch(128_MiB)) {
    for (std::size_t i{0}; i < batch->RecordCount(); ++i) {
        const RawRecord view{batch->RecordData(i)};
    }
}
```

Views within a batch are invalidated when the batch is destroyed.

### Supporting types

**CigarOp** — 32-bit value matching BAM binary layout (`op_len << 4 | op`).
Constexpr operations for reference/query length, string conversion, and bin
calculation.

**SequenceView** — Decodes 4-bit packed sequence on demand. Supports
indexing and iteration without unpacking the whole sequence.

**TagMap** — Ordered key-value collection for auxiliary tags. Linear-scan
lookup (records rarely exceed ~30 tags, so this beats a hash map on cache
effects). `TagKey` stores two characters as a `uint16_t`.

**SamHeader** — Parsed SAM header with `@HD`, `@SQ`, `@RG`, `@PG`, `@CO`
sections. Provides `ReferenceId(name)` and `ReferenceName(id)` lookups.
Round-trips between SAM text and BAM binary.

## API

### Readers

**BamRawReader** provides zero-copy view-based iteration:

1. **Range** — `for (const auto& view : reader.Records())` — wraps batch
   reading internally, simplest to use.

2. **Batch** — `reader.ReadBatch(limit)` — returns `RawRecordBatch` for
   explicit parallel processing control.

Region queries via `reader.Query(index, refId, beg, end)` return a filtered
range using a BAI index.

**BamRecordReader** wraps `BamRawReader` and delivers pre-decoded owned
`BamRecord` objects. A background producer thread reads batches, decodes
views to records (optionally in parallel), and pushes them through an SPSC
queue. Supports range-based iteration via `reader.Records()`.

**SamReader** reads text SAM files, producing owned `BamRecord` objects (not
views, since there's no persistent binary buffer to view into).

### Writers

**BamWriter** writes BGZF-compressed BAM. Accepts `BamRecord` (serializes),
`RawRecord`, raw `span<const byte>`, and `RawRecordBatch` (via
`WriteBatch()`). An optional `IndexCallback` is invoked after each record
write for downstream index building.

**SamWriter** writes text SAM. Accepts `BamRecord`, `RawRecord`, and
`RawRecordBatch` (via `WriteBatch()`).

### Indexing

**BaiIndex** reads, writes, builds, and queries BAM index files (`.bai`).
`Build()` scans a coordinate-sorted BAM to create an index.
`Query(refId, beg, end)` returns merged, sorted chunks of virtual offsets
for records overlapping a genomic region.

### PacBio extensions

**ZmwIndex** — In-memory index mapping ZMW hole numbers to virtual offsets.
Reads `.zmi` (native format) or `.pbi` (legacy PacBio index). Supports
point queries (`Find(zmw)`, `Find(ZmwIdentity)`), batch queries
(`Find(span<const ZmwIdentity>)`), and chunking via `UniqueZmws()`.

**ZmwWhitelist** — A set of ZMW hole numbers or identities to select from
a BAM file. Resolves against a `ZmwIndex` to produce sorted virtual offsets.
Used via `BamRawReader::Whitelist()` for seek-per-record iteration.

**ZmiWriter** — Streaming writer for `.zmi` files. Appends `(rgId, zmw,
virtualOffset)` entries in BAM write order. BGZF-compressed.

**BamZmwReader** — Groups consecutive records by ZMW identity, producing
groups of `BamRecord` objects (one group per ZMW). Wraps a `BamRecordReader`
and prefetches complete groups into a bounded queue on a background thread.

**ZmiBamWriter** — Writes a BAM and its `.zmi` index simultaneously.

## Reader Pipeline

The three-layer reader stack:

```
 BamZmwReader::GetNext()        (consume prefetched ZMW groups)
      │
      ▼
 BamZmwReader jthread           (group records, fill ZMW queue)
      │
      ▼
 BamRecordReader::ReadRecord()  (pops from SPSC queue)
      │
      └── BamRecordReader jthread (producer) ───────────────────┐
                │                                               │
                ▼                                               │
          BamRawReader::ReadBatch()                             │
                │                                               │
                ▼                                               │
          BgzfReader (sync or parallel via BgzfWorkers)         │
                │                                               │
                ▼                                               │
          Parallel::Dispatch (optional ThreadPool)              │
                │  ┌── ToOwned() per view ───────────┐          │
                │  │  decode views to BamRecord      │          │
                │  └─────────────────────────────────┘          │
                │                                               │
                └──────── push BamRecord → SPSC queue ──────────┘
```

**BamZmwReader** runs a producer thread that calls `BamRecordReader::ReadRecord()`
in a loop, accumulates records with the same ZMW hole number parsed from the
read name (`movie/zmw/start_end`), and publishes each complete group to a
bounded queue. `GetNext()` consumes those prefetched groups, so downstream work
can overlap with grouping of upcoming ZMWs.

---

### Chunking

Chunking is built into `BamRawReader` via `BamRawReaderConfig::ChunkNum`
and `TotalChunks`. When set, the reader automatically partitions the file by
ZMW identities using the `.zmi` or `.pbi` index alongside the BAM:

```cpp
// Read chunk 2 of 4 — reader resolves boundaries automatically
BamRawReader reader{bamPath, BamRawReaderConfig{.ChunkNum = 2, .TotalChunks = 4}};
while (auto view = reader.ReadRecord()) { /* … */ }
```

`ChunkNum` is 1-based (1 through `TotalChunks`). Internally the reader calls
`ZmwIndex::UniqueZmws()` to get sorted ZMW identities, partitions them
proportionally, then seeks to `ZmwIndex::FirstOffset()` and sets an internal
`RecordLimit`.

Setting `ChunkingMode = ChunkMode::SCATTER` keeps the same exact-partition
guarantee but spreads each chunk's ZMWs across the whole file (a
chaotic-but-deterministic sample) instead of one contiguous range. The unique
ZMWs are grouped into tiles of up to `ChunkTileZmws` consecutive ZMWs, the tiles
are shuffled by a seeded balanced Fisher–Yates permutation (`ChunkSeed`), and
this chunk takes its slice of the shuffled order. Each tile is read as a seek to
its first record plus its exact record count (counts, not virtual-offset
boundaries, because sync `BgzfReader::Tell()` is block-granular). Scatter reads
synchronously; both `ReadRecord` and `ReadBatch` honor the scatter plan (each
batch stays within one tile), so the `BamRecordReader` pre-decode path works in
scatter mode as well.

`RecordLimit` is also available directly to stop after an exact record count
at a specific file position:

```cpp
BamRawReader reader{bamPath, BamRawReaderConfig{.RecordLimit = 1000}};
reader.Seek(startOffset);
while (auto view = reader.ReadRecord()) { /* exactly 1000 records */ }
```

Both limits are enforced in `ReadRecord()` and `ReadBatch()`.

## Error Handling

Two error policy types are defined in `<pbsamoa/io/BamRawReader.hpp>`:

- **ThrowPolicy** — `OnCorruptRecord()` is `[[noreturn]]`; throws
  `std::runtime_error` on corrupt records.
- **SkipPolicy** — `OnCorruptRecord()` silently skips; `SkippedCount()`
  returns the total number of records skipped.

## Data Flow

A typical read pipeline:

```
file on disk
  → BgzfReader (decompress BGZF blocks, optionally parallel)
    → BamRawReader (parse record boundaries from decompressed bytes)
      → RawRecordBatch (owns buffer, provides RecordData(i) spans)
        → RawRecord (eager CIGAR copy, other fields decode-on-demand)
          → BamRecord (optional: ToOwned() for mutation)
```

A typical write pipeline:

```
BamRecord or RawRecord or raw span<const byte>
  → BamWriter (serialize record if needed, accumulate into BGZF blocks)
    → BgzfWriter (compress blocks with libdeflate)
      → file on disk
```

## Clipping

`BamRecord` supports coordinate-based clipping via `Clip()` (in-place) and
`Clipped()` (copy). Two clip types: `CLIP_TO_QUERY` (polymerase/ZMW
coordinates) and `CLIP_TO_REFERENCE` (genomic coordinates).

Clipping adjusts CIGAR, sequence, qualities, and optionally auxiliary tags.
Tag clipping is driven by a `TagClipper` registry that maps `TagKey` values
to `TagClipStrategy` subclasses:

- **SubstringClipStrategy** — simple substring for per-base tags (`ip`, `pw`, etc.)
- **ReverseSubstringClipStrategy** — mirrored offset for reverse-orientation tags (`ri`, `rp`)
- **PulseClipStrategy** — pulse-space mapping via the `pc` tag
- **BasemodClipStrategy** — MM/ML base-modification tags per SAMv1.6
- **PileupClipStrategy** — RLE-encoded `sa` pileup coverage tag

`TagClipper::PacBioDefault()` returns a pre-configured clipper with all
PacBio-standard strategies registered.

## Metrics

All readers and writers expose `GetMetrics()` for runtime performance
monitoring.

- **BgzfMetrics** — IO throughput, decompression pool stats, SPSC queue
  depth, stall counters, cumulative stage timing (nanoseconds)
- **DecodeMetrics** — decode pool stats, throughput, stall counters, timing
- **ReaderMetrics** — composite of `BgzfMetrics` + `DecodeMetrics` +
  `TotalRecordsRead`, returned by `BamRecordReader::GetMetrics()`
- **ZmwReaderMetrics** — `BamZmwReader` prefetch queue capacity, occupancy,
  throughput, and stall counters
- **BgzfWriteMetrics** — compression pool stats, stall counters, throughput,
  timing
- **WriterMetrics** — composite of `BgzfWriteMetrics` +
  `TotalRecordsWritten`, returned by `BamWriter::GetMetrics()`

`BamRawReader::GetMetrics()` returns `BgzfMetrics` directly. `BamZmwReader`
returns `ZmwReaderMetrics` from its group-prefetch layer. All snapshots are
captured via cheap state reads; diff two snapshots to compute per-second
rates where appropriate.
