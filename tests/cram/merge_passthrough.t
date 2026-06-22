Disjoint-chain passthrough: a coordinate merge whose inputs occupy strictly
disjoint coordinate ranges is emitted via verbatim BGZF block copy (no
recompression), reordered by minimum coordinate, and still matches `samtools merge`.

Split the sorted fixture into a low-coordinate half and a high-coordinate half so
the two inputs do not overlap:

  $ samtools view -H "${TESTDIR}"/../data/many_records.bam > hdr.sam
  $ samtools view "${TESTDIR}"/../data/many_records.bam > body.sam
  $ head -n 250 body.sam > lo_body.sam
  $ tail -n 250 body.sam > hi_body.sam
  $ cat hdr.sam lo_body.sam | samtools sort -o lo.bam 2>/dev/null
  $ cat hdr.sam hi_body.sam | samtools sort -o hi.bam 2>/dev/null

The merge reports the passthrough path, even with the high-range input given first
(it is reordered by minimum coordinate):

  $ "${PBSAMOA}" merge pb.bam hi.bam lo.bam 2>&1 >/dev/null
  pbsamoa merge: 500 records from 2 input(s) (passthrough)
  *runtime:* (glob)
  *records 500* (glob)
  *input *output* (glob)
  *no decode/compress pipeline* (glob)
  *verdict: I/O-COPY bound* (glob)

Output matches `samtools merge` and stays coordinate-sorted with every record:

  $ samtools merge -f st.bam lo.bam hi.bam 2>/dev/null
  $ samtools view pb.bam > pb.sam
  $ samtools view st.bam > st.sam
  $ diff pb.sam st.sam
  $ samtools view -c pb.bam
  500
  $ samtools view -H pb.bam | grep "^@HD"
  @HD\tVN:1.6\tSO:coordinate (esc)

A single input is a trivial disjoint chain and also passes through:

  $ "${PBSAMOA}" merge one.bam lo.bam 2>&1 >/dev/null
  pbsamoa merge: 250 records from 1 input(s) (passthrough)
  *runtime:* (glob)
  *records 250* (glob)
  *input *output* (glob)
  *no decode/compress pipeline* (glob)
  *verdict: I/O-COPY bound* (glob)
