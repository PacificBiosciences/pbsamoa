Read SAM written by samtools as BAM:

  $ samtools view -bS "${TESTDIR}"/../data/diverse.sam > from_samtools.bam
  $ "${PBSAMOA}" dump from_samtools.bam | grep -cv "^@"
  13

pbsamoa SAM output is valid for samtools:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/diverse.bam > pbsamoa.sam
  $ samtools view -bS pbsamoa.sam > reconverted.bam
  $ samtools view reconverted.bam | wc -l | sed 's/ //g'
  13
