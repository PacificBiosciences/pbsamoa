Test SAM write from BAM input matches samtools:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/spec_example.bam > pbsamoa.sam
  $ samtools view -h --no-PG "${TESTDIR}"/../data/spec_example.bam > samtools.sam
  $ diff pbsamoa.sam samtools.sam
