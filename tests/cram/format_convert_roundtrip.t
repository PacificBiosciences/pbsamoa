pbsamoa BAM->SAM->BAM round-trip preserves record count:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/diverse.bam > first.sam
  $ samtools view -bS first.sam | samtools sort -o second.bam
  $ "${PBSAMOA}" dump second.bam | grep -cv "^@"
  13

Unsorted BAM can be read:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/unsorted.bam | grep -cv "^@"
  3
