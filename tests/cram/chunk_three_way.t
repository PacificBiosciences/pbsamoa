Create test BAM with 5 ZMWs of varying record counts (11 total):

  $ printf '@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:ref\tLN:10000\n' > input.sam
  $ printf 'movie/10/0_100\t0\tref\t101\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n' >> input.sam
  $ printf 'movie/20/0_100\t0\tref\t201\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n' >> input.sam
  $ printf 'movie/20/100_200\t0\tref\t301\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n' >> input.sam
  $ printf 'movie/20/200_300\t0\tref\t401\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n' >> input.sam
  $ printf 'movie/30/0_100\t0\tref\t501\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n' >> input.sam
  $ printf 'movie/30/100_200\t0\tref\t601\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n' >> input.sam
  $ printf 'movie/40/0_100\t0\tref\t701\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n' >> input.sam
  $ printf 'movie/50/0_100\t0\tref\t801\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n' >> input.sam
  $ printf 'movie/50/100_200\t0\tref\t901\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n' >> input.sam
  $ printf 'movie/50/200_300\t0\tref\t1001\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n' >> input.sam
  $ printf 'movie/50/300_400\t0\tref\t1101\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n' >> input.sam
  $ samtools view -bS input.sam > input.bam
  $ "${PBSAMOA}" zmi-build input.bam indexed.bam 2>/dev/null
  $ "${PBSAMOA}" dump indexed.bam | grep -v "^@" | sort > full.txt

Three chunks — no overlap or gaps:

  $ for i in 1 2 3; do "${PBSAMOA}" chunk indexed.bam $i 3 | grep -v "^@"; done | sort > chunked3.txt
  $ diff chunked3.txt full.txt

Three-chunk ZMW boundary preservation:

  $ "${PBSAMOA}" chunk indexed.bam 1 3 | grep -v "^@" | cut -d/ -f2 | sort -u > z1.txt
  $ "${PBSAMOA}" chunk indexed.bam 2 3 | grep -v "^@" | cut -d/ -f2 | sort -u > z2.txt
  $ "${PBSAMOA}" chunk indexed.bam 3 3 | grep -v "^@" | cut -d/ -f2 | sort -u > z3.txt
  $ cat z1.txt z2.txt z3.txt | sort | uniq -d | wc -l | sed 's/ //g'
  0
