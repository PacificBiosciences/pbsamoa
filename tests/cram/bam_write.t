Test BAM write: copy spec_example.bam through pbsamoa, verify samtools can read it:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/spec_example.bam > expected.sam
  $ samtools view -h --no-PG "${TESTDIR}"/../data/spec_example.bam > samtools_expected.sam
  $ diff expected.sam samtools_expected.sam

Write BAM with pbsamoa zmi-build, then verify with samtools:

  $ "${PBSAMOA}" zmi-build "${TESTDIR}"/../data/spec_example.bam output.bam
  $ samtools view -h --no-PG output.bam > from_copy.sam
  $ samtools view -h --no-PG "${TESTDIR}"/../data/spec_example.bam > from_original.sam
  $ diff from_copy.sam from_original.sam
