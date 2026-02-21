Copy then dump matches samtools on the copy (single-block BAM):

  $ "${PBSAMOA}" zmi-build "${TESTDIR}"/../data/spec_example.bam copy.bam 2>/dev/null
  $ "${PBSAMOA}" dump copy.bam | grep -v "^@PG" > pipeline.sam
  $ samtools view -h copy.bam | grep -v "^@PG" > samtools.sam
  $ diff pipeline.sam samtools.sam
