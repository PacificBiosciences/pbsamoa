Read SAM written by samtools as BAM:

  $ samtools view -bS "${TESTDIR}"/../data/diverse.sam > from_samtools.bam
  $ "${PBSAMOA}" dump from_samtools.bam | grep -cv "^@"
  13

pbsamoa SAM output is valid for samtools:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/diverse.bam > pbsamoa.sam
  $ samtools view -bS pbsamoa.sam > reconverted.bam
  $ samtools view reconverted.bam | wc -l | sed 's/ //g'
  13

pbsamoa BAM->SAM->BAM round-trip preserves record count:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/diverse.bam > first.sam
  $ samtools view -bS first.sam | samtools sort -o second.bam
  $ "${PBSAMOA}" dump second.bam | grep -cv "^@"
  13

Unsorted BAM can be read:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/unsorted.bam | grep -cv "^@"
  3
