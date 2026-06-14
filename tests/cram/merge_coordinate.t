Coordinate merge: split a file into two coordinate-sorted halves, then merge.
pbsamoa merge must match `samtools merge` record bodies and set SO:coordinate.

Build two sorted inputs sharing one header:

  $ samtools view -H "${TESTDIR}"/../data/many_records.bam > hdr.sam
  $ samtools view "${TESTDIR}"/../data/many_records.bam > body.sam
  $ head -n 250 body.sam > a_body.sam
  $ tail -n 250 body.sam > b_body.sam
  $ cat hdr.sam a_body.sam | samtools sort -o a.bam 2>/dev/null
  $ cat hdr.sam b_body.sam | samtools sort -o b.bam 2>/dev/null

Merge with pbsamoa and with samtools, compare record bodies:

  $ "${PBSAMOA}" merge pb.bam a.bam b.bam 2>/dev/null
  $ samtools merge -f st.bam a.bam b.bam 2>/dev/null
  $ samtools view pb.bam > pb.sam
  $ samtools view st.bam > st.sam
  $ diff pb.sam st.sam

All 500 records are present:

  $ samtools view -c pb.bam
  500

Header sort order and an appended @PG:

  $ samtools view -H pb.bam | grep "^@HD"
  @HD\tVN:1.6\tSO:coordinate (esc)
  $ samtools view -H pb.bam | grep -c "ID:pbsamoa.merge"
  1
