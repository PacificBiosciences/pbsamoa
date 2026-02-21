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

### copy — BAM to BAM

Copy a BAM file.

```sh
pbsamoa copy input.bam output.bam
```

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

### parse — parse-only benchmark

Parse a BAM file without producing output. Reports throughput in MiB/s
and records/s.

```sh
pbsamoa parse input.bam
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

Run throughput benchmarks on a BAM file.

```sh
pbsamoa bench input.bam
```

### bgzf-cat — decompress BGZF

Read a BGZF file block by block and write decompressed data to stdout.

```sh
pbsamoa bgzf-cat compressed.bam > raw.bin
```

### bgzf-compress — compress to BGZF

Compress stdin to BGZF format.

```sh
cat data.bin | pbsamoa bgzf-compress data.bgzf
```
