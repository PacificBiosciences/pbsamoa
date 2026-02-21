Test BAM write: copy spec_example.bam through pbsamoa, verify samtools can read it:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/spec_example.bam > expected.sam
  $ samtools view -h --no-PG "${TESTDIR}"/../data/spec_example.bam > samtools_expected.sam
  $ diff expected.sam samtools_expected.sam
