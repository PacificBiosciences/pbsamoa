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

Three chunks cover all records with no overlap or gaps:

  $ "${PBSAMOA}" chunk indexed.bam 1 3 | grep -v "^@" | sort > chunk1.txt
  $ "${PBSAMOA}" chunk indexed.bam 2 3 | grep -v "^@" | sort > chunk2.txt
  $ "${PBSAMOA}" chunk indexed.bam 3 3 | grep -v "^@" | sort > chunk3.txt
  $ cat chunk1.txt chunk2.txt chunk3.txt | sort > chunked.txt
  $ "${PBSAMOA}" dump indexed.bam | grep -v "^@" | sort > full.txt
  $ diff chunked.txt full.txt

Single chunk reads everything:

  $ "${PBSAMOA}" chunk indexed.bam 1 1 | grep -v "^@" | wc -l | sed 's/ //g'
  6
