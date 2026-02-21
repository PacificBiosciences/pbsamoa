Test SAM read round-trip (read SAM, output SAM):

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/spec_example.sam > pbsamoa_out.sam
  $ samtools view -h --no-PG "${TESTDIR}"/../data/spec_example.sam > samtools_out.sam
  $ diff pbsamoa_out.sam samtools_out.sam
