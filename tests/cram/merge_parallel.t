Parallel read-ahead merge: --threads and --memory must produce a correct,
deterministic coordinate merge matching `samtools merge`.

Build two coordinate-sorted inputs sharing one header. The bodies interleave by
line so both inputs span the full coordinate range — their ranges overlap, which
forces the heap read-ahead path (not the disjoint-chain passthrough):

  $ samtools view -H "${TESTDIR}"/../data/many_records.bam > hdr.sam
  $ samtools view "${TESTDIR}"/../data/many_records.bam > body.sam
  $ awk 'NR % 2 == 1' body.sam > a_body.sam
  $ awk 'NR % 2 == 0' body.sam > b_body.sam
  $ cat hdr.sam a_body.sam | samtools sort -o a.bam 2>/dev/null
  $ cat hdr.sam b_body.sam | samtools sort -o b.bam 2>/dev/null

The overlapping inputs take the heap merge, not passthrough:

  $ "${PBSAMOA}" merge --threads 4 chk.bam a.bam b.bam 2>&1 >/dev/null | grep passthrough | wc -l | tr -d ' '
  0

Parallel merge with a tiny read-ahead budget (forces backpressure) matches samtools:

  $ "${PBSAMOA}" merge --threads 4 --memory 1K pb.bam a.bam b.bam 2>/dev/null
  $ samtools merge -f st.bam a.bam b.bam 2>/dev/null
  $ samtools view pb.bam > pb.sam
  $ samtools view st.bam > st.sam
  $ diff pb.sam st.sam

Deterministic across thread count and budget (single thread, large budget):

  $ "${PBSAMOA}" merge --threads 1 --memory 256M pb1.bam a.bam b.bam 2>/dev/null
  $ samtools view pb1.bam > pb1.sam
  $ diff pb.sam pb1.sam

All 500 records present and SO:coordinate preserved:

  $ samtools view -c pb.bam
  500
  $ samtools view -H pb.bam | grep "^@HD"
  @HD\tVN:1.6\tSO:coordinate (esc)

Granular decode/compress/batch/writer-queue knobs reproduce the samtools merge:

  $ "${PBSAMOA}" merge --decode-threads 4 --compress-threads 1 --batch-bytes 1K --writer-queue 8 g.bam a.bam b.bam 2>/dev/null
  $ samtools view g.bam > g.sam
  $ diff g.sam st.sam
  $ samtools view -c g.bam
  500

An invalid --memory value is rejected:

  $ "${PBSAMOA}" merge --memory bogus out.bam a.bam b.bam > /dev/null 2>&1
  [1]

--writer-queue 0 is rejected (queue must be >= 1):

  $ "${PBSAMOA}" merge --writer-queue 0 out.bam a.bam b.bam > /dev/null 2>&1
  [1]

A negative --decode-threads is rejected:

  $ "${PBSAMOA}" merge --decode-threads -1 out.bam a.bam b.bam > /dev/null 2>&1
  [1]
