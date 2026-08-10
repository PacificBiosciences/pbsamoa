Repeated pbsamoa invocations must keep @PG IDs unique, as required by SAMv1 §1.3.

Sort twice; the second @PG entry gets a numeric suffix and links to the first:

  $ "${PBSAMOA}" sort "${TESTDIR}"/../data/many_records.bam s1.bam > /dev/null 2>&1
  $ "${PBSAMOA}" sort s1.bam s2.bam > /dev/null 2>&1
  $ "${PBSAMOA}" dump s2.bam > s2.sam
  $ grep "^@PG" s2.sam | cut -f2
  ID:samtools
  ID:pbsamoa.sort
  ID:pbsamoa.sort.1
  $ grep "ID:pbsamoa.sort.1" s2.sam | tr '\t' '\n' | grep "^PP:"
  PP:pbsamoa.sort

Merge twice without colliding with the shared @PG entry:

  $ "${PBSAMOA}" merge m1.bam s1.bam s1.bam > /dev/null 2>&1
  $ "${PBSAMOA}" merge m2.bam m1.bam m1.bam > /dev/null 2>&1
  $ "${PBSAMOA}" dump m2.bam > m2.sam
  $ grep "^@PG" m2.sam | cut -f2
  ID:samtools
  ID:pbsamoa.sort
  ID:pbsamoa.merge
  ID:pbsamoa.merge.1

The concat merge path also keeps the ID unique:

  $ "${PBSAMOA}" merge --concat c1.bam s1.bam s1.bam > /dev/null 2>&1
  $ "${PBSAMOA}" merge --concat c2.bam c1.bam c1.bam > /dev/null 2>&1

Read the file back before counting its merge IDs:

  $ "${PBSAMOA}" dump c2.bam > /dev/null
  $ "${PBSAMOA}" dump c2.bam | grep -c "ID:pbsamoa.merge"
  2

samtools reads the header and agrees with pbsamoa on every record:

  $ samtools view -H m2.bam | grep -c "ID:pbsamoa.merge"
  2
  $ "${PBSAMOA}" dump m2.bam | grep -v "^@" > pb.sam
  $ samtools view m2.bam > st.sam
  $ diff pb.sam st.sam

pbsamoa must also read legacy files with duplicate @PG IDs. Forge one:

  $ samtools view -H s1.bam > hdr.sam
  $ sed -n '/ID:pbsamoa\.sort/p' hdr.sam | sed 's/PN:pbsamoa/PN:other_tool/' >> hdr.sam
  $ samtools reheader hdr.sam s1.bam > dup.bam
  $ samtools view -H dup.bam | grep -c "ID:pbsamoa.sort"
  2

dump keeps both lines in file order:

  $ "${PBSAMOA}" dump dup.bam | grep "^@PG" | grep -c "ID:pbsamoa.sort"
  2
  $ "${PBSAMOA}" dump dup.bam | grep "ID:pbsamoa.sort" | grep -c "PN:other_tool"
  1

sort and merge read it and use the next free suffix:

  $ "${PBSAMOA}" sort dup.bam dupsorted.bam > /dev/null 2>&1
  $ "${PBSAMOA}" dump dupsorted.bam | grep "^@PG" | cut -f2 | grep "^ID:pbsamoa.sort"
  ID:pbsamoa.sort
  ID:pbsamoa.sort
  ID:pbsamoa.sort.1
  $ "${PBSAMOA}" merge --bai dupmerged.bam dup.bam s1.bam > /dev/null 2>&1
  $ "${PBSAMOA}" dump dupmerged.bam > /dev/null
