Subsample unaligned HiFi reads to a target yield. The fixture has 302 records with
321748 bases in 300 ZMWs; two ZMWs have separate fwd and rev strand records.

  $ HIFI="${TESTDIR}"/../data/hifi.bam
  $ cp "${HIFI}" nopbi.bam

With a .pbi, the tool reads the lengths from the index and the BAM once. The kept
yield reaches the target and overshoots by less than one read:

  $ "${PBSAMOA}" subsample --yield 100K --seed 5 "${HIFI}" pbi.bam
  pbsamoa subsample: kept 91 of 302 records, 100536 of 321748 bases (target 100000, seed 5, 1 BAM pass with .pbi)

Without a .pbi, an extra pass over the BAM gets the lengths. The kept records are the
same:

  $ "${PBSAMOA}" subsample --yield 100K --seed 5 nopbi.bam twopass.bam
  pbsamoa subsample: kept 91 of 302 records, 100536 of 321748 bases (target 100000, seed 5, 2 BAM passes, no .pbi)
  $ samtools view pbi.bam > pbi.sam
  $ samtools view twopass.bam > twopass.sam
  $ cmp pbi.sam twopass.sam

The selection keys on (read group, ZMW), so both strands of a ZMW stay together:

  $ "${PBSAMOA}" subsample --yield 160K --seed 1 "${HIFI}" out.bam 2>/dev/null
  $ samtools view out.bam | cut -f1 | grep -E 'fwd|rev'
  m84000_240101_000000_s1/2631155/ccs/fwd
  m84000_240101_000000_s1/2631155/ccs/rev
  m84000_240101_000000_s1/33131177/ccs/fwd
  m84000_240101_000000_s1/33131177/ccs/rev

Suffixes are powers of 1000, and a larger target keeps a superset:

  $ "${PBSAMOA}" subsample --yield 0.05m --seed 5 "${HIFI}" small.bam
  pbsamoa subsample: kept 46 of 302 records, 51621 of 321748 bases (target 50000, seed 5, 1 BAM pass with .pbi)
  $ samtools view small.bam | cut -f1 | sort > small.names
  $ samtools view pbi.bam | cut -f1 | sort > pbi.names
  $ comm -23 small.names pbi.names

The header records the target and the seed. The auto seed comes from the header:

  $ "${PBSAMOA}" subsample --yield 100K "${HIFI}" auto.bam 2>/dev/null
  $ samtools view -H auto.bam | grep '^@CO'
  @CO\tSub-sampled yield=100000 seed=* (esc) (glob)
  $ samtools view -H auto.bam | grep -c '^@PG.*ID:pbsamoa.subsample'
  1

A target above the input yield keeps every record:

  $ "${PBSAMOA}" subsample --yield 1G --seed 5 "${HIFI}" all.bam
  pbsamoa subsample: WARNING: the input has 321748 bases, less than the target yield 1000000000; keeping all records
  pbsamoa subsample: kept 302 of 302 records, 321748 of 321748 bases (target 1000000000, seed 5, 1 BAM pass with .pbi)

A stale .pbi is an error, and the tool removes the partial output:

  $ samtools view -h "${HIFI}" | awk 'NR != 10' | samtools view -b --no-PG -o stale.bam -
  $ cp "${HIFI}".pbi stale.bam.pbi
  $ "${PBSAMOA}" subsample --yield 100K stale.bam stale.out.bam 2>&1 | grep -o 'the .pbi does not match BAM record [0-9]*'
  the .pbi does not match BAM record 6
  $ test -e stale.out.bam || echo removed
  removed

Read names must be PacBio names movie/zmw/...:

  $ "${PBSAMOA}" subsample --yield 1K "${TESTDIR}"/../data/benchmark.bam out.bam 2>&1 | grep -o -- '--yield needs PacBio read names.*'
  --yield needs PacBio read names movie/zmw/..., found: read00000

Invalid yields:

  $ for y in 60X abc 0 1.5.2G; do
  >   "${PBSAMOA}" subsample --yield "${y}" "${HIFI}" out.bam 2>&1 | grep -o 'invalid --yield: [^ ]*'
  > done
  invalid --yield: 60X
  invalid --yield: abc
  invalid --yield: 0
  invalid --yield: 1.5.2G
