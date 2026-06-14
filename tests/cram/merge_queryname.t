Query-name merge: two queryname-sorted halves merged with --order queryname must
match `samtools merge -n` record bodies and set SO:queryname.

  $ samtools view -H "${TESTDIR}"/../data/many_records.bam > hdr.sam
  $ samtools view "${TESTDIR}"/../data/many_records.bam > body.sam
  $ head -n 250 body.sam > a_body.sam
  $ tail -n 250 body.sam > b_body.sam
  $ cat hdr.sam a_body.sam | samtools sort -n -o a.bam 2>/dev/null
  $ cat hdr.sam b_body.sam | samtools sort -n -o b.bam 2>/dev/null

  $ "${PBSAMOA}" merge --order queryname pb.bam a.bam b.bam 2>/dev/null
  $ samtools merge -n -f st.bam a.bam b.bam 2>/dev/null
  $ samtools view pb.bam > pb.sam
  $ samtools view st.bam > st.sam
  $ diff pb.sam st.sam

  $ samtools view -H pb.bam | grep -o "SO:queryname"
  SO:queryname
