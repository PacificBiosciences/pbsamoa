Test BAM read: pbsamoa dump should produce the same output as samtools view -h:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/spec_example.bam > pbsamoa_out.sam
  $ samtools view -h --no-PG "${TESTDIR}"/../data/spec_example.bam > samtools_out.sam
  $ diff pbsamoa_out.sam samtools_out.sam
