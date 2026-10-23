# Architecture

pbsamoa has five layers. The public headers are in four directories:
`core/` (Compression and Record), `cram/`, `index/`, and `io/` (the API).

```
┌────────────────────────────────────────────────────────┐
│  BAM/SAM API (io/)                                     │
│  BamRawReader  BamRecordReader  SamReader  BamFile     │
│  BamWriter  SamWriter  BamCollection  SortBam MergeBam │
│  MergeReadAhead  CramReader  CramWriter                │
│  BamZmwReader  ZmiWriter  ZmiBamWriter                 │
├────────────────────────────────────────────────────────┤
│  Indexes (index/)                                      │
│  BaiIndex  BaiStreamBuilder  CraiIndex                 │
│  ZmwIndex  ZmwWhitelist                                │
├────────────────────────────────────────────────────────┤
│  CRAM v3 (cram/)                                       │
│  CramStructs  CramCodec  CramCompression  CramReference│
├────────────────────────────────────────────────────────┤
│  Record (core/)                                        │
│  RawRecordView  RawRecord  BamRecord  RawRecordBatch   │
│  CigarOp  SequenceView  TagMap  SamHeader  TagClipper  │
│  GenomicInterval  Basemods  CigarClipping              │
│  LosslessClipping                                      │
├────────────────────────────────────────────────────────┤
│  Compression (core/)                                   │
│  BgzfReader  BgzfWriter  BgzfWriterPool  VirtualOffset │
│  BgzfBlockEncoder  BgzfBlockWriter                     │
└────────────────────────────────────────────────────────┘
```

The `core/` headers include only other `core/` headers. The `cram/` headers
include only `cram/CramStructs.hpp`. Of the `index/` headers, only `BaiIndex.hpp`
includes a `core/` header. The `io/` headers use all the layers.
One exception exists: `src/BamRecord.cpp` (core) includes `io/BamWriter.hpp`,
because it defines `EncodeBamRecord()`, which `io/BamWriter.hpp` declares.

### Dependencies

- **pbcopper** gives `Parallel::ThreadPool` and `Parallel::Dispatch`. Most
  worker pools use them. `BaiIndex::Build` uses `Parallel::FireAndForgetIndexed`,
  and `CramWriter` writes containers with `std::async`. The CLI tools also use
  pbcopper CLIv2.
- **libdeflate** does the DEFLATE compression, the decompression, and the CRC32 of BGZF
  and of CRAM gzip blocks.
- **third-party/rigtorp** gives `rigtorp::SPSCQueue`, the lock-free queue of
  the BGZF pipelines and of `BamRecordReader`.
- **third-party/htscodecs** gives the CRAM rANS 4x8, rANS 4x16, adaptive
  arithmetic, fqzcomp, and name tokeniser codecs. On x86-64, meson adds the
  AVX2 rANS and PACK kernels when the compiler accepts `-mavx2` and `cpuid.h`
  has `__cpuid_count` and `__get_cpuid_max`.
- **ISA-L** (optional, meson option `isal`, x86-64 with NASM only) gives
  `BgzfBlockEncoder` its level 0. With `auto` (the default), meson builds it
  on x86-64 when it finds NASM and the ISA-L source archive; an offline build
  skips it. With `enabled`, the build stops when meson cannot build ISA-L.
  The `pbsamoa` dependency exports the variable `isal`.
- **bzip2** and **liblzma** (optional) enable the CRAM block methods `BZIP2`
  and `LZMA` when meson finds them. Without them, the registry does not have
  these methods.
- **threads** (always) and **libm** (when the system has it) are linked too.
- The meson option `tests` (default `true`) builds the tests. It pulls in
  googletest and pbtestutils through their wraps.

## Class Dependencies

Arrows show composition ("owns") or required dependencies ("uses"). Dotted
arrows show decode-on-demand relationships. Inheritance is rare. The public
hierarchies are `TagClipStrategy` (see `TagClipping.hpp`), `CramCodec` (the
CRAM encodings, see `CramCodec.hpp`), and the `CramReference` interface.

```
┌─ API ──────────────────────────────────────────────────────────────────┐
│                                                                        │
│  BamRawReader ──┬──► BgzfReader (sync or parallel via numWorkers)      │
│                 ├──► SamHeader                                         │
│                 ├──► BamRawReader (per file, for a BamCollection)      │
│                 └──► RawRecordBatch                                    │
│                       └──► RecordData(i) → span<const byte>            │
│                                                                        │
│  BamRecordReader ┬──► BamRawReader                                     │
│                  ├──► SPSCQueue<expected<BamRecord, exception_ptr>>    │
│                  ├──► std::jthread (background producer)               │
│                  └──► ThreadPool (optional, for parallel ToOwned)      │
│                                                                        │
│  BamZmwReader ┬──► BamRecordReader                                     │
│               └──► std::jthread + deque<ZmwGroup> (ZMW prefetch)       │
│                                                                        │
│  BamWriter ─────► BgzfWriter (writes the SamHeader once, keeps none)   │
│                                                                        │
│  SamReader ──┬──► SamHeader                                            │
│              └──► BamRecord (produces per record)                      │
│                                                                        │
│  SamWriter ─────► SamHeader                                            │
│                                                                        │
│  CramReader ──┬──► SamHeader                                           │
│               ├──► CramCodec + CramCompression                         │
│               ├──► CramReference (optional, reference-based decode)    │
│               ├──► ThreadPool (optional, container decode read-ahead)  │
│               └──► CraiIndex (for region query)                        │
│                                                                        │
│  CramWriter ──┬──► SamHeader                                           │
│               ├──► CramCodec + CramCompression                         │
│               ├──► CramReference (optional, reference-based encode)    │
│               ├──► ThreadPool (optional, slice compression)            │
│               └──► optional CRAI writer                                │
│                                                                        │
│  BaiIndex ──────► ReferenceIndex ──┬──► bins: Chunk (VirtualOffset)    │
│                                    └──► linearIndex: VirtualOffset     │
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
│             └────► BgzfWriterPool (optional, shared by many writers)   │
│                                                                        │
└────────────────────────────────────────────────────────────────────────┘
```

Not shown: `BamFile`, `BamCollection`, `SortBam`, `MergeBam`, `MergeReadAhead`,
`BgzfBlockEncoder`, `BgzfBlockWriter`, `BaiStreamBuilder`, `ZmwIndex`,
`ZmwWhitelist`, `TagClipper` with its `TagClipStrategy` classes, and
`CompressionRegistry`. The sections below describe them.

Key relationships:

- **BamRawReader** owns a `BgzfReader` (synchronous when `BgzfWorkers == 0`,
  parallel pipeline when `BgzfWorkers > 0`, synchronous always in scatter
  chunking) for decompression, and a
  `SamHeader` parsed from the BAM header block. It produces `RawRecordBatch`
  objects that own their backing buffers and provide indexed access to raw
  record spans via `RecordData(i)`. Supports chunked reading via
  `ChunkNum`/`TotalChunks` in config (see [Chunking](#chunking)) and ZMW
  whitelist filtering via `Whitelist()`. The config field
  `BamRawReaderConfig::Whitelist` does not filter records. For one file, the
  constructor only checks that it is not set with chunking. Setting only one of
  `ChunkNum` and `TotalChunks` throws `std::invalid_argument`. For a
  collection, chunking or a set `Whitelist` throws `std::invalid_argument`, and
  `Whitelist()`, `Query(index, …)`, `Seek()`, and `Tell()` throw
  `std::logic_error`.

- **BamFile** holds a BAM path, its header, and the path of its `.bai`
  sidecar. It can make the `.bai` (`CreateStandardIndex()`).
  **BamCollection** holds an ordered list of `BamFile` inputs and one merged,
  compatible header. `BamRawReader`, `BamRecordReader`, and `BamZmwReader`
  accept a `BamCollection`, a list of paths, or a list of `BamFile` objects.
  The reader then reads the files one after the other, in list order.

- **BamRecordReader** wraps a `BamRawReader` and pre-decodes views into
  owned `BamRecord` objects in a background `std::jthread`. Records are
  pushed through an SPSC queue for consumption. Optional `ThreadPool`
  parallelizes the `ToOwned()` decode step across batch views.
  The defaults are `DecodeWorkers{4}`, `BatchBudget` 4 MiB, and
  `OutputCapacity{1024}` queue slots. `TagFilter` (`DropTags` or `KeepTags`)
  applies to each decode. The queue items are
  `std::expected<BamRecord, std::exception_ptr>`, so a producer error reaches
  the caller in record order. `Query(GenomicInterval)` uses the sibling `.bai`
  of each input. For a collection, it merges the inputs on
  `(refId, pos, name, input order)`.

- **BamZmwReader** wraps a `BamRecordReader` and groups consecutive records
  by `(rgId, zmw)`. The hole number comes from the read name (`movie/zmw/…`).
  The `rgId` comes from the base PacBio read-group ID of the `RG` tag, or is 0
  without one. A background `std::jthread` prefetches complete groups into a
  bounded `std::deque`. `PrefetchCapacityZmws` sets its capacity (default 4).
  `CurrentZmw()` returns the identity of the last returned group.

- **BamWriter** owns a `BgzfWriter` for compression. The constructor writes
  the header block once and does not keep the `SamHeader`. Accepts `BamRecord`, `RawRecordView`, `RawRecord`,
  raw `span<const byte>`, and `RawRecordBatch` (via `WriteBatch()`).

- **CramReader** owns CRAM container/slice decode state and produces owned
  `BamRecord` objects. Container payload parsing is centralized through
  `ParseContainer(...)`, then slices are decompressed/decoded in order.
  It supports optional reference-based decoding and region queries via
  `CraiIndex` (`Query`, `QueryRaw`, `ReadSlice`, `ReadSliceRaw`).
  `ReadRawRecord()`, `RawRecords()`, and the `*Raw` queries give `RawRecord`
  objects. The decode makes each record a `RawRecord` immediately after its
  decode, so a slice is not held two times. The reader then resolves attached
  mates in place with `RawRecord::SetNext` and `RawRecord::SetTlen`.
  See [CRAM](#cram) for the threads and the reference.

- **CramWriter** uses a two-level write buffer:
  `pendingRecords` (records for the current slice) and `pendingSlices`
  (completed slices for the current container). `RecordsPerSlice` controls
  slice boundaries; `SlicesPerContainer` controls container flush boundaries.
  Slice encode/compress can run in parallel, while container writes are kept
  ordered with an async write future. It can emit a samtools-compatible CRAI
  sidecar with per-slice offsets. With reference sequences it writes
  reference-based CRAM. `EncodeContainer` encodes a container on any thread,
  and `WriteContainer` writes it in order. `Serialize` also runs on any
  thread, and `WriteSerialized` writes its bytes in order.

- **ZmiBamWriter** composes a `BamWriter` and a `ZmiWriter` so every record
  written is simultaneously indexed.

- **SortBam / MergeBam** (free functions in `io/BamSort.hpp`, `io/BamMerge.hpp`)
  combine files rather than streaming a single one. `SortBam` is a RAM-bounded
  external merge sort (spill sorted runs, then the k-way `MergeRuns` core).
  `SortConfig::MaxMemory` (default 768 MiB) sets the run size.
  `SortConfig::MergeFanIn` (default 64) sets the maximum number of open runs.
  With more runs, `MergeRuns` merges them in more passes. The intermediate
  files use compression level 1. `MergeBam` detects the mode
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
  budget cannot stall the merge. The defaults are `ReadAheadMemory` 768 MiB,
  `BatchBytes` 256 KiB for each handoff batch, and `WriterQueueCapacity` 256
  output blocks. A producer reads 32 BGZF blocks at a time. `DecodeThreads` and
  `CompressThreads` set the two pools. They default to `NumThreads`, and
  0 means min(hardware threads, 8). Per-source file order is preserved, so the merged
  output is independent of thread count, decode timing, and budget. Avoiding the
  per-record allocation is what lets input decompression actually overlap the
  merge — on a 16-input sorted benchmark with 8,535,828 records,
  `pbsamoa merge --threads 20 --memory 5G` ran in 2m16s versus 8m28s for
  `samtools 1.23 merge -@20` (~3.7x faster), and `--threads 32 --memory 5G`
  ran in 1m03s (~8.0x faster than that samtools run).

- **RawRecordView** borrows raw BAM bytes. **RawRecord** owns a copy. `CigarOp` values are eagerly
  copied into an aligned buffer on construction (BAM does not guarantee
  alignment). `SequenceView` and `TagMap` are decoded on demand. All decode
  logic is inline in the header for maximum performance.

## Compression

BAM files use BGZF (Blocked GZip Format) — a series of independently
compressed gzip blocks, each at most 64 KiB decompressed. This enables
random access via virtual offsets.

CRAM uses container/slice structures and block-level codecs selected from
CRAM v3.x methods (`raw`, `gzip`, `bzip2`, `lzma`, `rANS`, `arith`,
`fqzcomp`, `name tokeniser`) via `CramCompression`. See [CRAM](#cram).

### BgzfReader

Reads BGZF blocks sequentially, decompressing each with libdeflate.
Supports `Seek()` to any `VirtualOffset`.

With `numWorkers > 0`, an IO thread reads groups of 32 blocks, a thread pool
decompresses each group into one buffer, and a consumer thread finds the record
boundaries in the group. The group buffer becomes a `RawRecordBatch`: `ReadBatch()`
gives it to the caller and does not copy the records. A record that spans groups
goes behind the data of the group where it ends. `ReadRecord()` copies one record
out of the current batch. The IO thread and the consumer thread are
`std::jthread` objects that the reader owns. The pool is a pbcopper
`ThreadPool` with `numWorkers` threads. A `rigtorp::SPSCQueue` of 8 batches
(about 2 MiB each) connects the consumer thread to the caller.

With `numWorkers == 0`, the reader decompresses on the calling thread and
starts no threads. `Tell()` then gives the exact virtual offset of the next
record. The one-argument constructor opens a block-only reader for index files.

### BgzfWriter

Accumulates data into blocks of at most 65280 bytes (`BGZF_MAX_UNCOMPRESSED_BLOCK`)
and compresses each with libdeflate in a parallel writer pipeline. Compression
level 1–12 (default 6). Appends the standard 28-byte EOF marker on close.

`BgzfWriterConfig` holds `CompressionLevel{6}`, `BgzfWorkers{4}`,
`InputQueueCapacity{256}`, `BlocksPerBatch{32}`, `UseTempFile`, and `Pool`.
Each writer owns a packer thread, an IO writer thread, and a pbcopper
`ThreadPool` of `BgzfWorkers` threads. A `rigtorp::SPSCQueue` with
`InputQueueCapacity` slots connects the caller to the packer thread. The
packer sends `BlocksPerBatch` blocks to the pool as one task.

With `Pool` set to a shared `BgzfWriterPool`, the writer packs blocks on the
calling thread and starts no threads of its own. The pool has one set of
compression workers and one ordered output thread for all its writers. This
mode is for many open writers, for example demultiplexed output.

When a callback is configured (via `BamWriter`), callback dispatch happens on
the BGZF IO writer thread after block file offsets are known, not on the caller
thread. Callbacks can be deferred until `Close()`.

### BgzfBlockEncoder and BgzfBlockWriter

`BgzfBlockEncoder` compresses bytes into complete BGZF blocks on the calling
thread. `BgzfBlockWriter` takes these blocks from many threads and writes them
in the order of the `Write()` calls on one IO thread. Use the pair when the
producer threads must do the compression, for example when an aligner holds an
output lock: the lock then covers only a queue push. `Write()` blocks while the
queue holds `MaxPendingBytes`, so slow storage stops the producers. The writer does not track
virtual offsets and has no index callback. `EncodeBamRecord()` adds a BAM
record to an encoder. `MaxPendingBytes` is 64 MiB by default.

`BgzfBlockEncoder` takes levels 0–12 (default 6). Levels 1–12 use libdeflate.
Level 0 uses ISA-L igzip level 1 when the build has ISA-L. Without ISA-L,
level 0 uses libdeflate level 1. `BgzfWriter` does not accept level 0.
The ISA-L wrap applies `subprojects/packagefiles/isa-l-write-first-byte-hash.patch`:
without it, the compressed bytes of level 0 depend on the address of the ISA-L stream
struct and change from run to run (the data stays the same).

### VirtualOffset

A 64-bit value encoding position within a BGZF file: upper 48 bits are the
compressed block offset, lower 16 bits are the offset within the
decompressed block. Supports comparison via `<=>` but not arithmetic — this
is intentional per the BAM spec.

## Record

Three representations for alignment records, chosen by use case:

### RawRecordView and RawRecord — fast, read-only

`RawRecordView` borrows raw BAM bytes. `RawRecord` owns a copy for
single-record APIs and for code that crosses threads. CIGAR ops are eagerly copied into
aligned storage on construction (BAM does not guarantee 4-byte alignment).
Other accessors decode fields on demand directly from the binary layout.
All decode logic is inline in the header for maximum inlining.
`RawRecord::View()` makes an owner's borrowed lifetime explicit. You can call it
only on an lvalue `RawRecord`, so a temporary `RawRecord` cannot silently produce a dangling view.

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

`RawRecordView::ToOwned()` and `RawRecord::ToOwned()` produce a `BamRecord`. Optional tag filters
(`DropTags`, `KeepTags`) let you skip large tags during conversion. This
matters for PacBio BAM files where kinetics arrays can be megabytes per
record.

`BamWriter` and `SamWriter` accept `BamRecord`, `RawRecordView`, and `RawRecord`
(plus `WriteBatch()` for `RawRecordBatch`). `BamWriter` also accepts raw
`span<const byte>`. `BamWriter` does not serialize a raw view that you pass to it.
It serializes a `BamRecord` that you pass to it.

### RawRecordBatch

A batch owns the decompressed buffer and provides indexed, non-owning decoders
via `View(i)` plus raw spans via `RecordData(i)`. Batches are sized by memory budget
(`ByteLimit`), not record count. The default budget is 256 MiB.

```cpp
using namespace PacBio::Samoa::Literals;
while (auto batch = reader.ReadBatch(128_MiB)) {
    for (std::size_t i{0}; i < batch->RecordCount(); ++i) {
        const RawRecordView view{batch->View(i)};
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

**TagMap** — Key-value collection for auxiliary tags, a vector sorted by
key. Lookups are a binary search (`std::ranges::lower_bound`).
`TagKey` stores two characters as a `uint16_t`.

**SamHeader** — Parsed SAM header with `@HD`, `@SQ`, `@RG`, `@PG`, `@CO`
sections. Provides `ReferenceId(name)` and `ReferenceName(id)` lookups.
Round-trips between SAM text and BAM binary.

## API

### Readers

**BamRawReader** provides decode-on-demand iteration:

1. **Range** — `for (const auto& view : reader.Records())` — wraps batch
   reading internally, simplest to use. Yields owning `RawRecord` references.

2. **Batch** — `reader.ReadBatch(limit)` — returns `RawRecordBatch` for
   explicit parallel processing control. `batch.View(i)` decodes a record
   and does not copy its bytes.

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
`RawRecordView`, `RawRecord`, raw `span<const byte>`, and `RawRecordBatch` (via
`WriteBatch()`). An optional `IndexCallback` gets each record's offsets for
downstream index building. It runs on the BGZF IO writer thread and can run
after `Write()` returns, as late as `Close()`.

**SamWriter** writes text SAM. Accepts `BamRecord`, `RawRecordView`,
`RawRecord`, and `RawRecordBatch` (via `WriteBatch()`).

### Indexing

**BaiIndex** reads, writes, builds, and queries BAM index files (`.bai`).
`Build()` scans a coordinate-sorted BAM to create an index.
`Query(refId, beg, end)` returns merged, sorted chunks of virtual offsets
for records overlapping a genomic region. `Build(bamPath, numWorkers)`
decompresses blocks in parallel when `numWorkers > 1`.

**BaiStreamBuilder** builds a `BaiIndex` while records go by.
`Observe(beginOffset, recordBytes)` adds one record. `Finalize(endOffset)`
gives the index. `Build()` and the merge option `BaiOutput` give the same
bytes for the same records.

**CraiIndex** reads a gzip-compressed `.crai` file. `Query(refId, beg, end)`
gives one row for each slice that overlaps the region, in file order.
`CramReader::ReadSlice` decodes the slice of one row.

### PacBio extensions

**ZmwIndex** — In-memory index mapping ZMW hole numbers to virtual offsets.
Reads `.zmi` (native format) or `.pbi` (legacy PacBio index). Supports
point queries (`Find(zmw)`, `Find(ZmwIdentity)`), batch queries
(`Find(span<const ZmwIdentity>)`), and chunking via `UniqueZmws()` and
`ContiguousChunk()`. `Open(bamPath)` tries the `.zmi` first, then the `.pbi`.
A ZMI 1.1 index also holds barcodes: `FindBarcodes()` and `BarcodePairs()`
use them. `QueryLengths()` gives the per-record query lengths of a `.pbi`;
a `.zmi` has none.

**ZmwWhitelist** — A set of ZMW hole numbers or identities to select from
a BAM file. Resolves against a `ZmwIndex` to produce sorted virtual offsets.
Used via `BamRawReader::Whitelist()` for seek-per-record iteration.

**ZmiWriter** — Streaming writer for `.zmi` files. Appends `(rgId, zmw,
virtualOffset, forwardBarcode, reverseBarcode)` entries in BAM write order.
BGZF-compressed through a `BgzfWriter`. The format has a 64-byte header with
the magic `ZMI\1` and version 1.1.0, then one 24-byte entry for each record.
The first 16 bytes of an entry are the same as in ZMI 1.0.

**BamZmwReader** — Groups consecutive records by ZMW identity, producing
groups of `BamRecord` objects (one group per ZMW). Wraps a `BamRecordReader`
and prefetches complete groups into a bounded queue on a background thread.

**ZmiBamWriter** — Writes a BAM and its `.zmi` index simultaneously. The
`.zmi` path is `bamPath + ".zmi"`.

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
in a loop, accumulates records with the same `(rgId, zmw)` identity (the
read-group ID from the `RG` tag, the hole number from the read name
`movie/zmw/start_end`), and publishes each complete group to a
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
while (auto record = reader.ReadRecord()) { /* … */ }
```

`ChunkNum` is 1-based (1 through `TotalChunks`). Set both or neither: setting
only one throws `std::invalid_argument`. Internally the reader calls
`ZmwIndex::ContiguousChunk()`. This function scans the columnar index directly,
partitions the file-order ZMW runs proportionally, and returns one seek offset
and an exact record limit. This limit replaces any `RecordLimit` in the config.
Ordered BAM indexes require no hash tables or per-record lookup
vectors during chunk selection.

Setting `ChunkingMode = ChunkMode::SCATTER` keeps the same exact-partition
guarantee but spreads each chunk's ZMWs across the whole file (a
chaotic-but-deterministic sample) instead of one contiguous range. The unique
ZMWs are grouped into tiles of up to `ChunkTileZmws` consecutive ZMWs, the tiles
are shuffled by a seeded balanced Fisher–Yates permutation (`ChunkSeed`), and
this chunk takes its slice of the shuffled order. Each tile is read as a seek to
its first record plus its exact record count. The reader sorts the tiles of
its chunk into file order and joins adjacent tiles into one read. The defaults
are `ChunkTileZmws{100}` and `ChunkSeed{42}`. Scatter reads
synchronously; both `ReadRecord` and `ReadBatch` honor the scatter plan (each
batch stays within one run of joined tiles), so the `BamRecordReader` pre-decode path works in
scatter mode as well.

`RecordLimit` is also available directly to stop after an exact record count
at a specific file position:

```cpp
BamRawReader reader{bamPath, BamRawReaderConfig{.RecordLimit = 1000}};
reader.Seek(startOffset);
while (auto record = reader.ReadRecord()) { /* exactly 1000 records */ }
```

Both limits are enforced in `ReadRecord()` and `ReadBatch()`, in scatter mode too: a
`RecordLimit` from the config caps the records of the scatter chunk.

## CRAM

pbsamoa has its own CRAM v3 implementation. It does not use htslib.
`CramReader` reads CRAM 3.x. `CramWriter` writes CRAM 3.0 or 3.1
(`MinorVersion`). It writes 3.1 whenever a block uses a 3.1 method (`RANS4X16`,
`ADAPTIVE_ARITH`, `FQZCOMP`, or `NAME_TOKENISER`), whatever `MinorVersion`
says. A `MajorVersion` other than 3 throws.

### Codecs

`CramCodec` and its subclasses are the data-series encodings: `NullCodec`,
`ExternalCodec`, `HuffmanCodec`, `BetaCodec`, `GammaCodec`, `SubexpCodec`,
`ByteArrayLenCodec`, and `ByteArrayStopCodec`.

`CramCompression` holds `CompressionRegistry`, the block compression methods.
It registers the CRAM methods 0–8 (`CramBlockMethod`):

| ID | Method | Implementation |
|----|--------|----------------|
| 0 | `RAW` | built in |
| 1 | `GZIP` | libdeflate |
| 2 | `BZIP2` | bzip2, only when meson finds it |
| 3 | `LZMA` | liblzma, only when meson finds it |
| 4 | `RANS4X8` | htscodecs |
| 5 | `RANS4X16` | htscodecs |
| 6 | `ADAPTIVE_ARITH` | htscodecs |
| 7 | `FQZCOMP` | htscodecs |
| 8 | `NAME_TOKENISER` | htscodecs |

`CompressionRegistry::Register()` adds custom methods with IDs above 8. It
cannot replace the methods 0–8.

`CramWriterConfig::BlockCompressionMethod` (default `GZIP`) sets the method of
the data blocks. `DataSeriesCompressionMethods` sets a different method for
one data series, for example `QS`. With `NAME_TOKENISER` as the default, only
`RN` uses it; with `FQZCOMP`, only `QS` and `QQ` do. The other data series and
the non-data-series blocks then use `GZIP`. Setting `NAME_TOKENISER` on a series
other than `RN`, or `FQZCOMP` on one other than `QS`/`QQ`, throws
`std::runtime_error`. `RANS4X16` compresses `QS`/`QQ` with order 0 and PACK
(SIMD chosen at run time), and the other blocks with order 1.
`CompressionLevel` can stay empty, and gzip then uses level 6; when set, it is a
libdeflate level in [0, 12]. `Open()` throws for a level outside that range.

### Reference

A reference-based CRAM stores only the differences of a mapped record to the
reference. `CramReference` is the interface for the reference. It gives the
bases of a window of one reference sequence. Several threads can call it at
the same time, so many readers and writers can share one reference.

- `CramReaderConfig::ReferencePath` loads a whole FASTA file into memory. The
  reader matches its sequences to the `@SQ` lines by name. Without
  `ReferencePath` and `Reference`, the reader takes the path from the
  environment variable `PBSAMOA_CRAM_REFERENCE`.
- `CramReaderConfig::Reference` and `CramWriterConfig::Reference` take a
  `CramReference`. The reader keeps an upper-cased copy of one window of it,
  and each decode task keeps its own. A window that misses reads at least
  64 Ki bases.
- `CramWriterConfig::ReferenceSequences` takes one `std::string_view` for each
  `@SQ` line, in header order.

Set one reference source only. The writer puts the MD5 of the reference bases
in each slice. It adds an `M5` tag to each `@SQ` line that has none.
`CramReferenceMd5()` computes the `M5` value. The reader checks the reference
MD5 of each single-reference slice that needs the reference and has a non-zero
MD5. A mismatch, or a missing reference, throws `std::runtime_error`. A query checks each slice one time
for each reader. Without a reference, the writer stores all the bases of each
record.

### Reader threads

`CramReaderConfig::DecompressionWorkers` (default 0) sets the size of a pbcopper
`ThreadPool`. With 0, the reader decodes all containers on the calling thread.

With workers, the sequential read works as follows:

```
 caller thread                         ThreadPool (DecompressionWorkers)
 ─────────────                         ─────────────────────────────────
 read container bytes from the file
   │  (up to 2 containers per worker)
   ├──► submit decode task ─────────────► decode container to BamRecord
   │    push future to readAhead deque      or RawRecord, with its own
   │                                        ReferenceWindow
   ▼
 pop the front future (file order) ◄──── set_value / set_exception
 give its records one by one
```

- Only the calling thread reads the file.
- Each task decodes one whole container. The task decodes to the record type
  that the caller reads (`BamRecord` or `RawRecord`). If the caller switches
  between `ReadRecord()` and `ReadRawRecord()`, the containers already queued
  keep their type, and their records are converted one by one.
- The futures stay in file order, so the records come out in file order.
- A task gives its exception to its future. `ReadRecord()` throws it in file
  order.
- The pool is the last member of the reader, so it joins its tasks first.

`Query`, `QueryRaw`, `ReadSlice`, and `ReadSliceRaw` read on the same stream
as the sequential read. The sequential read then continues unchanged. These
calls use the pool for the blocks of one slice when there are 2 or more
workers and 8 or more blocks. In a coordinate-sorted, single-reference slice,
the query decode stops after the region end when no record waits for a mate.

### Writer threads

`CramWriterConfig` sets the layout: `RecordsPerSlice{10000}` and
`SlicesPerContainer{1}`. A value of `SlicesPerContainer` <= 0 flushes only at
`Close()`. `CompressionWorkers` (default 0) sets the size of a pbcopper
`ThreadPool`. With 2 or more workers, the pool compresses the slices of a
container in parallel when it has more than one slice. Otherwise (the default
`SlicesPerContainer{1}`), it compresses the blocks of a slice in parallel when
the slice has 8 or more blocks. With workers, `std::async` writes each
container, so the next container can encode while the previous one writes.
With exactly 1 worker, only this write runs in parallel.

`EncodeContainer()` and `Serialize()` are `const` and run on any thread.
`WriteContainer()` and `WriteSerialized()` write the result in order. An
aligner can then encode on its own threads and keep only the write in a lock.

`WriteCrai` writes a samtools-compatible `.crai` file next to the output, or at
`CraiPath`.

## Error Handling

Two error policy types are defined in `<pbsamoa/io/BamRawReader.hpp>`:

- **ThrowPolicy** — `OnCorruptRecord()` is `[[noreturn]]`; throws
  `std::runtime_error` on corrupt records.
- **SkipPolicy** — `OnCorruptRecord()` silently skips; `SkippedCount()`
  returns the total number of records skipped.

No reader takes a policy at this time. The readers throw `std::runtime_error`
on a corrupt record or a corrupt block. The parallel pipelines catch an
exception on a worker thread and throw it again on the caller thread.

## Data Flow

A typical read pipeline:

```
file on disk
  → BgzfReader (decompress BGZF blocks, optionally parallel, and parse record
    boundaries from the decompressed bytes)
    → BamRawReader (chunking, record limits, queries)
      → RawRecordBatch (owns buffer, provides View(i) decoders and RecordData(i) spans)
        → RawRecordView (aligned CIGAR copy, other fields decode-on-demand)
          → BamRecord (optional: ToOwned() for mutation)
```

A typical write pipeline:

```
BamRecord or RawRecordView or RawRecord or raw span<const byte>
  → BamWriter (serialize record if needed, accumulate into BGZF blocks)
    → BgzfWriter (compress blocks with libdeflate)
      → file on disk
```

## Clipping

`BamRecord` supports coordinate-based clipping via `Clip()` (in-place) and
`Clipped()` (copy). Two clip types: `CLIP_TO_QUERY` (polymerase/ZMW
coordinates) and `CLIP_TO_REFERENCE` (genomic coordinates).

Clipping adjusts CIGAR, sequence, qualities, and auxiliary tags. Without a
clipper argument, `Clip()` uses `TagClipper::PacBioDefault()`.
Tag clipping is driven by a `TagClipper` registry that maps `TagKey` values
to `TagClipStrategy` subclasses:

- **SubstringClipStrategy** — simple substring for per-base tags (`ip`, `pw`, etc.)
- **ReverseSubstringClipStrategy** — mirrored offset for reverse-orientation tags (`ri`, `rp`)
- **PulseClipStrategy** — pulse-space mapping via the `pc` tag
- **BasemodClipStrategy** — MM/ML base-modification tags per SAMv1.6
- **PileupClipStrategy** — RLE-encoded `sa` pileup coverage tag

`TagClipper::PacBioDefault()` returns a pre-configured clipper with all
PacBio-standard strategies registered.

The clipping modules are in `core/`:

- **CigarClipping** — `ClipCigarToQuery()` and `ClipCigarToReference()` clip a
  CIGAR. They give a `ClipResult` with the query offset, the query length, and
  the new reference position.
- **Basemods** — `ParseBasemodString()` and `WriteBasemodString()` convert the
  `MM` tag. `BasemodClipStrategy` uses them. `ClipBasemodRecord()` finds which
  sites of one modification type fall before, in, and after a clip window.
  `ClipToQueryLossless()` uses it.
- **LosslessClipping** — `ClipToQueryLossless()` clips a record to a query
  range and keeps the removed flanks in a msgpack `ls` tag. The `ls` format is
  the same as in lima. `RestoreFromLossless()` uses the `ls` tag to restore
  the record, like `lima-undo`. It does not restore CIGAR or POS, so it is for
  unmapped reads.

## Metrics

The BGZF readers and writers expose `GetMetrics()` for runtime performance
monitoring: `BgzfReader`, `BgzfWriter`, `BamRawReader`, `BamRecordReader`,
`BamZmwReader`, and `BamWriter`. `SamReader`, `SamWriter`, `CramReader`, and
`CramWriter` have no metrics. `MergeBam` gives `MergeStats::Runtime`, and
`MergeReadAhead` gives `MergeReadAheadStats`.

- **BgzfMetrics** — IO throughput, decompression pool stats, SPSC queue
  produced/consumed counts (depth is their difference), stall counters,
  cumulative stage timing (nanoseconds)
- **DecodeMetrics** — decode pool stats, throughput, stall counters, timing
- **ReaderMetrics** — composite of `BgzfMetrics` + `DecodeMetrics` +
  `TotalRecordsRead`, plus the flags `ParallelBgzf` and `ParallelDecode`,
  returned by `BamRecordReader::GetMetrics()`
- **ZmwReaderMetrics** — `BamZmwReader` prefetch queue capacity, occupancy,
  throughput, and stall counters, plus `Reader`, the `ReaderMetrics` of the
  wrapped `BamRecordReader`
- **BgzfWriteMetrics** — compression pool stats, stall counters, throughput,
  timing
- **WriterMetrics** — composite of `BgzfWriteMetrics` +
  `TotalRecordsWritten`, returned by `BamWriter::GetMetrics()`

`BamRawReader::GetMetrics()` returns `BgzfMetrics` directly. `BamZmwReader`
returns `ZmwReaderMetrics`. All snapshots are
captured via cheap state reads; diff two snapshots to compute per-second
rates where appropriate.
