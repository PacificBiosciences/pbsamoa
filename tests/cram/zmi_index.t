Build a BAM with PacBio-style read names (3 ZMWs: 10, 20, 30):

  $ printf '@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:ref\tLN:10000\n' > input.sam
  $ printf 'movie/10/0_100\t0\tref\t101\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n' >> input.sam
  $ printf 'movie/10/100_200\t0\tref\t201\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n' >> input.sam
  $ printf 'movie/20/0_100\t0\tref\t301\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n' >> input.sam
  $ printf 'movie/30/0_100\t0\tref\t401\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n' >> input.sam
  $ printf 'movie/30/100_200\t0\tref\t501\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n' >> input.sam
  $ printf 'movie/30/200_300\t0\tref\t601\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n' >> input.sam
  $ samtools view -bS input.sam > input.bam

Round-trip the BAM through zmi-build to fix BGZF block layout, then drop its index:

  $ "${PBSAMOA}" zmi-build input.bam canonical.bam 2>/dev/null
  $ test -f canonical.bam.zmi
  $ mv canonical.bam.zmi canonical.bam.zmi.from-build

Run zmi-index over the same BAM and confirm the new sidecar:

  $ "${PBSAMOA}" zmi-index canonical.bam 2>/dev/null
  $ test -f canonical.bam.zmi

Output must byte-equal the zmi-build sidecar (same writer, same records, same VOs):

  $ cmp canonical.bam.zmi canonical.bam.zmi.from-build && echo MATCH
  MATCH

Determinism across thread counts:

  $ rm canonical.bam.zmi
  $ "${PBSAMOA}" zmi-index --threads 1 canonical.bam 2>/dev/null
  $ mv canonical.bam.zmi canonical.bam.zmi.t1
  $ "${PBSAMOA}" zmi-index --threads 8 canonical.bam 2>/dev/null
  $ cmp canonical.bam.zmi canonical.bam.zmi.t1 && echo MATCH
  MATCH

Missing input exits non-zero:

  $ "${PBSAMOA}" zmi-index does-not-exist.bam 2>err.log
  [1]
  $ grep -q 'input file not found' err.log && echo OK
  OK

Quiet flag suppresses stderr line:

  $ rm canonical.bam.zmi
  $ "${PBSAMOA}" zmi-index --quiet canonical.bam 2>quiet.log
  $ test ! -s quiet.log
  $ test -f canonical.bam.zmi
