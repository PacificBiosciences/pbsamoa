pbsamoa bai-build creates an index file:

  $ cp "${TESTDIR}"/../data/diverse.bam test.bam
  $ "${PBSAMOA}" bai-build test.bam 2>/dev/null
  $ test -f test.bam.bai
