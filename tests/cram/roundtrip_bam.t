Full round-trip: BAM -> SAM -> BAM -> SAM should produce identical non-PG output:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/diverse.bam > step1.sam
  $ samtools view -bS step1.sam | samtools sort -o step2.bam
  $ "${PBSAMOA}" dump step2.bam > step3.sam
  $ grep -v "^@PG" step1.sam > step1_nopg.txt
  $ grep -v "^@PG" step3.sam > step3_nopg.txt
  $ diff step1_nopg.txt step3_nopg.txt
