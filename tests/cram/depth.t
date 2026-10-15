pbsamoa depth: the mean depth in bins, with samtools depth -J as the oracle.

  $ same() { oracle "$1" "$2" > expected.txt; sums "$2" "$3" > got.txt; cmp -s expected.txt got.txt; }
  $ oracle() { samtools depth -J -a "$1" | awk -v b="$2" -v OFS='\t' '{ s[$1 OFS int(($2 - 1) / b)] += $3 } END { for (k in s) if (s[k] > 0) print k, s[k] }' | sort; }
  $ sums() { awk -v b="$1" -v OFS='\t' '!/^#/ { s = $4 * ($3 - $2); if (s > 0) printf "%s\t%d\t%.0f\n", $1, $2 / b, s }' "$2" | sort; }

The records cover M, =, X, D and N, a supplementary record, MAPQ 0, and the records that
samtools depth skips: duplicate, QC-fail, secondary and unmapped. chr3 has no records:

  $ printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:30\n@SQ\tSN:chr2\tLN:12\n@SQ\tSN:chr3\tLN:5\n' > cases.sam
  $ printf 'r1\t0\tchr1\t1\t60\t5M2D3M\t*\t0\t0\t*\t*\n' >> cases.sam
  $ printf 'r3\t1024\tchr1\t1\t60\t10M\t*\t0\t0\t*\t*\n' >> cases.sam
  $ printf 'r4\t512\tchr1\t1\t60\t10M\t*\t0\t0\t*\t*\n' >> cases.sam
  $ printf 'r5\t256\tchr1\t1\t60\t10M\t*\t0\t0\t*\t*\n' >> cases.sam
  $ printf 'r2\t16\tchr1\t4\t0\t2S3=1X2N4M\t*\t0\t0\t*\t*\n' >> cases.sam
  $ printf 'r6\t2048\tchr1\t20\t60\t8M3I3M\t*\t0\t0\t*\t*\n' >> cases.sam
  $ printf 'r8\t0\tchr2\t1\t60\t4M1D2M\t*\t0\t0\t*\t*\n' >> cases.sam
  $ printf 'r7\t4\tchr2\t5\t0\t*\t*\t0\t0\t*\t*\n' >> cases.sam
  $ samtools view -b --no-PG -o cases.bam cases.sam
  $ "${PBSAMOA}" depth --bin-size 4 cases.bam cases.tsv
  pbsamoa depth: 4 of 8 records counted, 36 depth bases in 13 bins of 4 bp
  $ cat cases.tsv
  #contig\tstart\tend\tmean_depth (esc)
  chr1\t0\t4\t1.25 (esc)
  chr1\t4\t8\t1.75 (esc)
  chr1\t8\t12\t1.25 (esc)
  chr1\t12\t16\t0.25 (esc)
  chr1\t16\t20\t0.25 (esc)
  chr1\t20\t24\t1 (esc)
  chr1\t24\t28\t1 (esc)
  chr1\t28\t30\t1 (esc)
  chr2\t0\t4\t1 (esc)
  chr2\t4\t8\t0.75 (esc)
  chr2\t8\t12\t0 (esc)
  chr3\t0\t4\t0 (esc)
  chr3\t4\t5\t0 (esc)
  $ same cases.bam 4 cases.tsv
  $ awk 'END { print NR }' expected.txt
  10

The sums match samtools on the other fixtures, with a bin that does not divide the
contig lengths, and with one and with many threads:

  $ for f in diverse spec_example benchmark; do
  >   for b in 7 1000; do
  >     for t in 0 4; do
  >       "${PBSAMOA}" depth --bin-size $b --threads $t "${TESTDIR}"/../data/$f.bam $f.tsv 2>/dev/null
  >       same "${TESTDIR}"/../data/$f.bam $b $f.tsv || echo "$f $b $t differs"
  >     done
  >   done
  > done

The status stream has one phase, depth, in reads:

  $ "${PBSAMOA}" depth --bin-size 4 --status-json status.jsonl cases.bam s.tsv 2>/dev/null
  $ sed -E 's/"(t|wall_ms|started_wall_ms|pid|cpu_s|wall_s|peak_rss_bytes)":[0-9.]+/"\1":N/g; s/"(hostname|version|commandline|path)":"[^"]*"/"\1":S/g' status.jsonl
  {"command":"depth","commandline":S,"hostname":S,"num_threads":*,"phases":["depth"],"pid":N,"schema_version":1,"started_wall_ms":N,"tool":"pbsamoa","type":"run_info","version":S} (glob)
  {"name":"depth","state":"start","t":N,"type":"phase","wall_ms":N}
  {"done":8,"phase":"depth","t":N,"total":8,"type":"progress","unit":"reads","wall_ms":N}
  {"metrics":{"bin_size":4,"bins":13,"counted_records":4,"depth_bases":36,"records":8},"name":"depth","state":"end","t":N,"type":"phase","wall_ms":N}
  {"metrics":{"bin_size":4,"bins":13,"counted_records":4,"depth_bases":36,"records":8},"outputs":[{"kind":"tsv","path":S}],"phase":"FINISHED","resources":{"cpu_s":N,"num_threads":*,"peak_rss_bytes":N,"wall_s":N},"t":N,"type":"final","wall_ms":N} (glob)

With a BAI, the progress has the total of the records from the index: the mapped and
unmapped counts of all references and the unplaced reads. The metrics give it as
index_records. Without a BAI, as above, or with a BAI older than the BAM, there is no total:

  $ samtools index cases.bam
  $ "${PBSAMOA}" depth --bin-size 4 --status-json bai.jsonl cases.bam s.tsv 2>/dev/null
  $ tail -1 bai.jsonl | grep -o '"index_records":[0-9]*,"records":[0-9]*'
  "index_records":8,"records":8
  $ touch -t 200001010000 cases.bam.bai
  $ "${PBSAMOA}" depth --bin-size 4 --status-json stale.jsonl cases.bam s.tsv 2>/dev/null
  $ grep -c index_records stale.jsonl
  0
  [1]

A bin size below 1 is an error:

  $ "${PBSAMOA}" depth --bin-size 0 cases.bam bad.tsv 2>&1 | grep -o 'bin-size must be at least 1: 0'
  bin-size must be at least 1: 0
  $ test -e bad.tsv
  [1]
