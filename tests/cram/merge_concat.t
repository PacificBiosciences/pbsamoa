Concat merge: byte-level BGZF block passthrough must reproduce `samtools cat`
record bodies, exercising the partial header/record boundary block plus the
verbatim copy of the remaining compressed blocks.

Forced --concat appends inputs verbatim; compare bodies against samtools cat:

  $ "${PBSAMOA}" merge --concat pb.bam "${TESTDIR}"/../data/many_records.bam "${TESTDIR}"/../data/many_records.bam 2>/dev/null
  $ samtools cat -o st.bam "${TESTDIR}"/../data/many_records.bam "${TESTDIR}"/../data/many_records.bam
  $ samtools view pb.bam > pb.sam
  $ samtools view st.bam > st.sam
  $ diff pb.sam st.sam

Output is a valid BAM with all records, marked SO:unsorted, one appended @PG:

  $ samtools quickcheck pb.bam && echo OK
  OK
  $ samtools view -c pb.bam
  1000
  $ samtools view -H pb.bam | grep -c "SO:unsorted"
  1
  $ samtools view -H pb.bam | grep -c "ID:pbsamoa.merge"
  1

Auto-detect: all-unsorted inputs concatenate with no --order/--concat:

  $ "${PBSAMOA}" merge auto.bam "${TESTDIR}"/../data/unsorted.bam "${TESTDIR}"/../data/unsorted.bam 2>/dev/null
  $ samtools quickcheck auto.bam && echo OK
  OK
  $ samtools view -c auto.bam
  6
  $ samtools view -H auto.bam | grep -c "SO:unsorted"
  1
