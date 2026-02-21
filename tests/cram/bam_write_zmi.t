Write BAM with pbsamoa zmi-build, then verify with samtools:

  $ "${PBSAMOA}" zmi-build "${TESTDIR}"/../data/spec_example.bam output.bam 2>/dev/null
  $ samtools view -h --no-PG output.bam > from_copy.sam
  $ samtools view -h --no-PG "${TESTDIR}"/../data/spec_example.bam > from_original.sam
  $ diff from_copy.sam from_original.sam
