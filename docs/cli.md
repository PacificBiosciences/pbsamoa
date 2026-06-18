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
| `--bgzf-threads N` | `-j` | Decompression worker threads (default: auto, max 10) |
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
| `--bgzf-threads N` | `-j` | BGZF decompression worker threads (default: auto, max 8) |
| `--decode-threads N` | | BamRecord decode worker threads (default: 4) |
