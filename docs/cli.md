# CLI Tools

pbsamoa ships a multi-tool binary with subcommands for common BAM/SAM
operations and BGZF utilities.

## pbsamoa

### dump — BAM/SAM to text

Convert a BAM or SAM file to human-readable SAM text on stdout.

```sh
pbsamoa dump input.bam
pbsamoa dump input.bam > output.sam
```

Options:

|         Flag         | Alias |                            Description                            |
| -------------------- | ----- | ----------------------------------------------------------------- |
| `--bgzf-threads N`   | `-j`  | BGZF decompression worker threads (default: auto, max 8)          |
| `--format-threads N` |       | SAM formatting thread pool size (default: max(hw_concurrency, 4)) |

Environment:

- `PBSAMOA_METRICS=1` — Print per-second pipeline metrics and a final summary to stderr.

### bai-query — region query

Extract records overlapping a genomic region using a BAI index.
Expects a `.bai` file alongside the BAM (auto-detected as `<bam>.bai`).

```sh
pbsamoa bai-query sorted.bam chr1:1000-2000
```

Arguments: BAM path, region (ref:start-end, 1-based inclusive).

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

### zmi-query — query by ZMW

Extract records for specific ZMW hole numbers.

```sh
pbsamoa zmi-query movie.bam 12345
```

### chunk — partitioned reading

Read a proportional partition of a BAM file, divided by unique ZMWs.

```sh
pbsamoa chunk movie.bam 1 4   # chunk 1 of 4 (1-based)
pbsamoa chunk movie.bam 4 4   # chunk 4 of 4 (last)
```

### bench — benchmarks

Run throughput benchmarks on a BAM file. Executes multiple benchmarks
in sequence: `sequential_read`, `batch_read`, `pipeline_raw_read`,
`record_reader` (serial), `record_reader` (parallel), `write`, and
`region_query`.

```sh
pbsamoa bench input.bam
```

Options:

|         Flag         | Alias |                       Description                        |
| -------------------- | ----- | -------------------------------------------------------- |
| `--bgzf-threads N`   | `-j`  | BGZF decompression worker threads (default: auto, max 8) |
| `--decode-threads N` |       | BamRecord decode worker threads (default: 4)             |
