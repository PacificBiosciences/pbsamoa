pbsamoa bai-build creates an index file:

  $ cp "${TESTDIR}"/../data/diverse.bam test.bam
  $ "${PBSAMOA}" bai-build test.bam 2>/dev/null
  $ test -f test.bam.bai

Query with pbsamoa-built index, pbsamoa matches samtools:

  $ cp "${TESTDIR}"/../data/diverse.bam query.bam
  $ "${PBSAMOA}" bai-build query.bam 2>/dev/null
  $ "${PBSAMOA}" bai-query query.bam chr1:1-1000 | grep -v "^@" | sort > pbsamoa_query.txt
  $ samtools view query.bam chr1:1-1000 | sort > samtools_query.txt
  $ diff pbsamoa_query.txt samtools_query.txt

Query specific region returns subset:

  $ "${PBSAMOA}" bai-query query.bam chr1:200-310 | grep -v "^@" | wc -l | sed 's/ //g'
  2

Query chr2 returns chr2 records:

  $ "${PBSAMOA}" bai-query query.bam chr2:1-500 | grep -v "^@" | wc -l | sed 's/ //g'
  2

Empty region returns no records:

  $ "${PBSAMOA}" bai-query query.bam chr1:999-1000 | grep -v "^@" | wc -l | sed 's/ //g'
  0
