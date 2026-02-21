Create test BAM with PacBio-style read names (3 ZMWs: 10, 20, 30):

  $ printf '@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:ref\tLN:10000\n' > input.sam
  $ printf 'movie/10/0_100\t0\tref\t101\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n' >> input.sam
  $ printf 'movie/10/100_200\t0\tref\t201\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n' >> input.sam
  $ printf 'movie/20/0_100\t0\tref\t301\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n' >> input.sam
  $ printf 'movie/30/0_100\t0\tref\t401\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n' >> input.sam
  $ printf 'movie/30/100_200\t0\tref\t501\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n' >> input.sam
  $ printf 'movie/30/200_300\t0\tref\t601\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n' >> input.sam
  $ samtools view -bS input.sam > input.bam
  $ "${PBSAMOA}" zmi-build input.bam indexed.bam 2>/dev/null

Query ZMW 30 returns 3 records:

  $ "${PBSAMOA}" zmi-query indexed.bam 30 2>/dev/null | grep -v "^@" | wc -l | sed 's/ //g'
  3

Query ZMW 20 returns 1 record:

  $ "${PBSAMOA}" zmi-query indexed.bam 20 2>/dev/null | grep -v "^@" | wc -l | sed 's/ //g'
  1

Query ZMW 10 returns 2 records:

  $ "${PBSAMOA}" zmi-query indexed.bam 10 2>/dev/null | grep -v "^@" | wc -l | sed 's/ //g'
  2

Query nonexistent ZMW returns 0 records:

  $ "${PBSAMOA}" zmi-query indexed.bam 999 2>/dev/null | grep -v "^@" | wc -l | sed 's/ //g'
  0
