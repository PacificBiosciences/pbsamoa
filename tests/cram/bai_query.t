Test region query matches samtools:

  $ samtools index "${TESTDIR}"/../data/spec_example.bam 2>/dev/null || true
  $ "${PBSAMOA}" bai-query "${TESTDIR}"/../data/spec_example.bam ref:1-45 | grep -v "^@" | sort > pbsamoa_records.txt
  $ samtools view --no-PG "${TESTDIR}"/../data/spec_example.bam ref:1-45 | sort > samtools_records.txt
  $ diff pbsamoa_records.txt samtools_records.txt
