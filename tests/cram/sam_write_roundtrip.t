Test round-trip: BAM -> SAM -> BAM -> SAM should be stable:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/spec_example.bam > first.sam
  $ samtools view -bS --no-PG first.sam > intermediate.bam
  $ "${PBSAMOA}" dump intermediate.bam > second.sam
  $ diff first.sam second.sam
