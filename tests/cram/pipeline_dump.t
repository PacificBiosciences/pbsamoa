Pipeline output matches samtools for spec_example.bam:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/spec_example.bam | grep -v "^@PG" > pipeline.sam
  $ samtools view -h "${TESTDIR}"/../data/spec_example.bam | grep -v "^@PG" > samtools.sam
  $ diff pipeline.sam samtools.sam

Pipeline output matches samtools for diverse.bam:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/diverse.bam | grep -v "^@PG" > pipeline.sam
  $ samtools view -h "${TESTDIR}"/../data/diverse.bam | grep -v "^@PG" > samtools.sam
  $ diff pipeline.sam samtools.sam

Pipeline output matches samtools for many_records.bam:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/many_records.bam | grep -v "^@PG" > pipeline.sam
  $ samtools view -h "${TESTDIR}"/../data/many_records.bam | grep -v "^@PG" > samtools.sam
  $ diff pipeline.sam samtools.sam

Pipeline output matches samtools for unsorted.bam:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/unsorted.bam | grep -v "^@PG" > pipeline.sam
  $ samtools view -h "${TESTDIR}"/../data/unsorted.bam | grep -v "^@PG" > samtools.sam
  $ diff pipeline.sam samtools.sam

Header-only BAM produces only header lines:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/header_only.bam | grep -v "^@PG" > pipeline.sam
  $ samtools view -h "${TESTDIR}"/../data/header_only.bam | grep -v "^@PG" > samtools.sam
  $ diff pipeline.sam samtools.sam
