Coordinate sort: pbsamoa sort must match `samtools sort` record bodies and set
SO:coordinate in the header.

  $ "${PBSAMOA}" sort "${TESTDIR}"/../data/unsorted.bam pb.bam 2>/dev/null
  $ samtools sort -o st.bam "${TESTDIR}"/../data/unsorted.bam 2>/dev/null

Record bodies are identical (BGZF framing may differ, so compare SAM text):

  $ samtools view pb.bam > pb.sam
  $ samtools view st.bam > st.sam
  $ diff pb.sam st.sam

Header sort order is coordinate:

  $ samtools view -H pb.bam | grep "^@HD"
  @HD\tVN:1.6\tSO:coordinate (esc)

An appended @PG records the sort program:

  $ samtools view -H pb.bam | grep -c "ID:pbsamoa.sort"
  1
