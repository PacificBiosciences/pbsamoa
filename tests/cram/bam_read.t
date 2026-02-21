Test BAM read: pbsamoa dump should produce the same output as samtools view -h:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/spec_example.bam > pbsamoa_out.sam
  $ samtools view -h --no-PG "${TESTDIR}"/../data/spec_example.bam > samtools_out.sam
  $ diff pbsamoa_out.sam samtools_out.sam

Test header-only BAM produces only header:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/header_only.bam | grep -c "^@"
  * (glob)

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/header_only.bam | grep -cv "^@"
  0
  [1]
