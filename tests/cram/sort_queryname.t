Query-name sort: pbsamoa sort --order queryname must match `samtools sort -n`
record bodies and set SO:queryname in the header.

  $ "${PBSAMOA}" sort --order queryname "${TESTDIR}"/../data/many_records.bam pb.bam 2>/dev/null
  $ samtools sort -n -o st.bam "${TESTDIR}"/../data/many_records.bam 2>/dev/null

  $ samtools view pb.bam > pb.sam
  $ samtools view st.bam > st.sam
  $ diff pb.sam st.sam

  $ samtools view -H pb.bam | grep "^@HD"
  @HD\tVN:1.6\tSO:queryname (esc)
