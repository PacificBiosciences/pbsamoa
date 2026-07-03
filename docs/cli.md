# CLI Tools

pbsamoa ships a multi-tool binary with subcommands for BAM/SAM/CRAM workflows,
indexing, and performance benchmarking.

## pbsamoa

### dump — BAM/SAM/CRAM to SAM text

Convert BAM, SAM, or CRAM input to SAM text on stdout.

```sh
pbsamoa dump input.bam
pbsamoa dump input.sam
pbsamoa dump input.cram
pbsamoa dump input.cram > output.sam
```

Options:

| Flag | Alias | Description |
| ---- | ----- | ----------- |
| `--bgzf-threads N` | | Decompression worker threads (default: auto, max 10) |
| `--format-threads N` | | SAM formatting thread pool size for BAM input (default: `max(hw_concurrency, 4)`) |
| `--reference ref.fa` | | Reference FASTA for reference-based CRAM decoding |
| `--region ref:start-end|*` | | CRAM-only region query (`*` selects unmapped) |
| `--index file.crai` | | CRAM-only CRAI path override (requires `--region`) |

Notes:

- `--region` and `--index` are only supported for CRAM input.
- If `--region` is set and `--index` is omitted, `dump` loads `<input>.crai`.

Environment:

- `PBSAMOA_METRICS=1` — Print per-second pipeline metrics and a final summary to stderr for BAM dump pipeline.

### convert — BAM/SAM to CRAM

Convert BAM or SAM input to CRAM output.

```sh
pbsamoa convert input.bam output.cram
pbsamoa convert input.sam output.cram
```

Options:

| Flag | Description |
| ---- | ----------- |
| `--records-per-slice N` | Maximum records per CRAM slice (default: `10000`) |
| `--bgzf-threads N` | BAM input decompression workers (default: auto, max 16) |
| `--convert-to-bam-record` | Opt in to decoding BAM input into `BamRecord` before CRAM write |
| `--decode-threads N` | BAM decode workers for `--convert-to-bam-record` mode (default: auto, max 16). Also enables that mode when set explicitly |
| `--compression-threads N` | CRAM block compression workers (default: auto, max 8) |
| `--block-compression METHOD` | Default CRAM block method |
| `--series-compression SERIES=METHOD` | Override method for a CRAM data series (repeatable) |
| `--write-crai` | Write samtools-compatible `<output>.crai` |

Supported `METHOD` values:

- `raw`, `gzip`, `bzip2`, `lzma`, `rans4x8`, `rans4x16`, `arith`, `fqzcomp`, `tok`
- Numeric method IDs `0` through `8`

Supported `SERIES` values:

- `BF`, `CF`, `RI`, `RL`, `AP`, `RG`, `RN`, `MF`, `NS`, `NP`, `TS`, `NF`, `TL`, `FN`, `FC`, `FP`, `DL`, `BB`, `QQ`, `BS`, `IN`, `RS`, `PD`, `HC`, `SC`, `MQ`, `BA`, `QS`

Notes:

- Output must use `.cram` extension.
- Input must use `.bam` or `.sam` extension.
- Default block compression method is `rans4x8`.
- BAM input defaults to `RawRecord` passthrough into `CramWriter`; use
  `--convert-to-bam-record` to use the decoded `BamRecord` path.
- Advanced `CramWriterConfig` fields (`SlicesPerContainer`,
  `CompressionLevel`, `UseTempFile`) are currently API-only and are not
  exposed as `convert` flags.

### sort — sort a BAM file

Sort a BAM file by coordinate, query name, or a 2-char tag, using a RAM-bounded
external merge sort (spill to temp run files, then a bounded-fan-in k-way merge —
multi-pass when the run count exceeds the open-file limit).

```sh
pbsamoa sort -o sorted.bam input.bam
pbsamoa sort --order queryname -o byname.bam input.bam
```

Options:

| Flag | Description |
| ---- | ----------- |
| `--order ORDER` | `coordinate` \| `queryname` \| `tag` (default: `coordinate`) |
| `--tag XX` | 2-char tag to sort by (required iff `--order tag`) |
| `--memory SIZE` | Peak run-buffer budget, `K`/`M`/`G` suffix (default: `768M`) |
| `--temp-dir DIR` | Directory for temporary run files (default: output directory) |
| `--threads N` | Worker threads, `0`/auto = `min(hw, 8)` |
| `--compression L` | Output BGZF level `[1,12]` (default: `6`) |

### merge — merge or concatenate BAM files

Combine several BAM files (which must share an identical `@SQ` reference list)
into one. The merged header unions their `@RG` / `@PG` / `@CO` records and an
appended `pbsamoa.merge` `@PG`. With no `--order` / `--concat`, the mode is
auto-detected from the inputs' `@HD SO`:

- all `coordinate` or all `queryname` → streaming k-way sorted merge of the
  inputs (each must already be sorted by that order); output is deterministic
  and independent of the thread count, decode timing, and `--memory` budget.
- all unsorted/unknown → byte-level concatenation.
- a mix → error (pass `--order` or `--concat`).

Performance note: `pbsamoa merge` is much faster than `samtools merge` on
sorted inputs. On a 16-input benchmark with 8,535,828 records, `samtools 1.23
merge -@20` took 8m28s, while `pbsamoa merge --threads 20 --memory 5G` took
2m16s (~3.7x faster). With `--threads 32 --memory 5G`, `pbsamoa merge`
finished in 1m03s (~8.0x faster than that samtools run).

The sorted merge reads every input in parallel and ahead of the merge: each
input is decompressed on a shared worker pool and framed into a per-input bounded
queue, while the heap consumes already-decoded records and the writer compresses
the output in parallel. `--memory` caps the total in-flight read-ahead across all
inputs (a single budget, independent of input count); each input always keeps its
head record available, so a budget below one record cannot stall the merge.

When a coordinate merge's inputs occupy strictly disjoint coordinate ranges — each
file's coordinates entirely precede the next's — the merge skips the heap entirely
and emits via verbatim BGZF block passthrough (the `--concat` machinery), reordered
by minimum coordinate, after `ProbeDisjointChain` verifies non-overlap and internal
sort order. This eliminates recompression — the dominant cost — so the common shard
/ per-chromosome / `chunk`→`merge` roundtrip runs several times faster than the heap
merge (and than `samtools merge`, which always recompresses). Overlapping, internally
unsorted, or query-name/tag merges use the heap path. The CLI summary notes
`(passthrough)` when this path is taken.

```sh
pbsamoa merge out.bam a.bam b.bam              # auto-detect
pbsamoa merge --order coordinate out.bam a.bam b.bam
pbsamoa merge --concat out.bam a.bam b.bam     # force concatenation
```

`--concat` concatenates record streams via BGZF block passthrough — compressed
blocks are copied input→output without decompress/recompress (only the single
block straddling each input's header/record boundary is rebuilt), which is the
fastest path. The output is marked `SO:unsorted` and its record count is not
reported (counting would require decompressing every block).

Options:

| Flag | Description |
| ---- | ----------- |
| `--order ORDER` | `coordinate` \| `queryname` \| `tag` (default: auto-detect) |
| `--concat` | Byte-concatenate inputs (BGZF passthrough), output `SO:unsorted` |
| `--tag XX` | 2-char tag inputs are sorted by (required iff `--order tag`) |
| `--threads N` | Default size of both CPU pools (input decode + output compress), `0`/auto = `min(hw, 8)`. Does **not** bound the per-input producer threads (one I/O-bound thread per input). |
| `--decode-threads N` | Input-decompression pool size; `0` = inherit `--threads`. Sorted merge only. |
| `--compress-threads N` | Output-compression pool size; `0` = inherit `--threads`. Sorted merge only. |
| `--memory SIZE` | Total input read-ahead budget, `K`/`M`/`G` suffix (default: `768M`); sorted merge only. Peak RAM also includes ~`NumInputs × --batch-bytes` (per-input current batches held outside the budget) plus the writer queue. |
| `--batch-bytes SIZE` | Per-input handoff batch target, `K`/`M`/`G` suffix (default: `256K`); sorted merge only. |
| `--writer-queue N` | Output writer queue depth in compressed blocks (default: `256`, must be `≥ 1`); sorted merge only. |
| `--compression L` | Output BGZF level `[1,12]` (default: `6`) |
| `--bai` | Write `<out>.bai` during the merge, built on the fly (no extra pass). Coordinate merges that take the heap path only — rejected (non-zero exit) for `--concat`, `queryname`/`tag` order, or the disjoint-chain passthrough, since none decodes records into an indexable coordinate-sorted stream. For those, run `bai-build` on the output instead. The index is byte-identical to `bai-build`. |

The thread/memory knobs (except `--compression`) apply to the sorted-merge path
only; `--concat` is single-threaded BGZF passthrough.

### bai-query — region query

Extract records overlapping a genomic region using a BAI index.
Expects a `.bai` file alongside the BAM (auto-detected as `<bam>.bai`).

```sh
pbsamoa bai-query sorted.bam chr1:1000-2000
```

Arguments: BAM path, region (`ref:start-end`, 1-based inclusive).

### bai-build — build BAI

Build a `.bai` index from a coordinate-sorted BAM file.

```sh
pbsamoa bai-build sorted.bam
# creates sorted.bam.bai
```

### zmi-build — build ZMW index

Copy a BAM file and build a `.zmi` ZMW index alongside it.

```sh
pbsamoa zmi-build input.bam output.bam
# creates output.bam and output.bam.zmi
```

### zmi-index — build `.zmi` for an existing BAM

Build a `.zmi` sidecar from an existing BAM without rewriting the BAM.
Uses parallel BGZF decompression. `--threads` defaults to
`hardware_concurrency()`.

```sh
pbsamoa zmi-index input.bam
# creates input.bam.zmi

pbsamoa zmi-index --threads 8 --quiet input.bam
```

### zmi-query — query by ZMW

Extract records for specific ZMW hole numbers.

```sh
pbsamoa zmi-query movie.bam 12345
```

### chunk — partitioned reading

Read a partition of a BAM file, divided by unique ZMWs. Records of a ZMW always
stay together (works for subreads and HiFi), and running all `TOTAL` chunks
visits every record exactly once.

```sh
pbsamoa chunk movie.bam 1 4   # chunk 1 of 4 (1-based)
pbsamoa chunk movie.bam 4 4   # chunk 4 of 4 (last)
```

By default each chunk is one contiguous ZMW range. Use `--mode scatter` for a
chaotic-but-deterministic partition that samples ZMWs across the whole file —
useful when a single chunk should be a representative sample of the dataset:

```sh
pbsamoa chunk movie.bam 1 4 --mode scatter            # spread across the file
pbsamoa chunk movie.bam 1 4 --mode scatter --tile 64  # 64 ZMWs/seek (fewer seeks)
pbsamoa chunk movie.bam 1 4 --mode scatter --seed 7   # reproducible shuffle seed
```

Scatter assigns tiles of up to `M` consecutive ZMWs (`--tile M`, default 1) to
chunks via a seeded balanced shuffle (`--seed S`, default 0). The same
`(file, CHUNK, TOTAL, M, S)` always yields the same chunk on every platform.
Larger `M` reads more ZMWs per seek (less disk seeking, coarser sampling).
Both modes require a ZMW index (`.zmi`/`.pbi`); see `zmi-build`.

### bench — benchmarks

Run throughput benchmarks on a BAM file. Executes multiple benchmarks
in sequence: `sequential_read`, `batch_read`, `pipeline_raw_read`,
`record_reader` (serial), `record_reader` (parallel), `write`, and
`region_query`.

```sh
pbsamoa bench input.bam
```

Options:

| Flag | Alias | Description |
| ---- | ----- | ----------- |
| `--bgzf-threads N` | | BGZF decompression worker threads (default: auto, max 10) |
| `--decode-threads N` | | BamRecord decode worker threads (default: 4) |
