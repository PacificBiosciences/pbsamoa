Remove the HiFi kinetics from the kept records. The fixture has no kinetics, so the test
adds them to each record: by strand (fi, fp, ri, rp), over both strands (ip, pw), and a
tag that stays (sm) between them.

  $ HIFI="${TESTDIR}"/../data/hifi.bam
  $ samtools view -h --no-PG "${HIFI}" \
  >   | awk 'BEGIN { OFS = "\t" } /^@/ { print; next } { print $0, "fi:B:C,1,2", "fp:B:C,3", "sm:B:C,9,9", "ri:B:C,4", "rp:B:C,5", "ip:B:C,6", "pw:B:C,7" }' \
  >   | samtools view -b --no-PG -o kinetics.bam -

With --strip-kinetics, the tool keeps the same records and bases:

  $ "${PBSAMOA}" subsample --yield 100K --seed 5 --strip-kinetics kinetics.bam stripped.bam
  pbsamoa subsample: kept 91 of 302 records, 100536 of 321748 bases (target 100000, seed 5, 2 BAM passes, no .pbi)
  $ "${PBSAMOA}" subsample --yield 100K --seed 5 kinetics.bam kept.bam 2>/dev/null

No kinetics tag stays. The other fields and tags, and their order, do not change:

  $ samtools view stripped.bam | grep -cE '	(fi|fp|ri|rp|ip|pw):'
  0
  [1]
  $ samtools view stripped.bam | grep -c '	sm:B:C,9,9'
  91
  $ samtools view kept.bam | sed -E 's/	(fi|fp|ri|rp|ip|pw):B:C(,[0-9]+)+//g' > expected.sam
  $ samtools view stripped.bam > stripped.sam
  $ cmp expected.sam stripped.sam

The fraction mode removes them too:

  $ "${PBSAMOA}" subsample --fraction 0.5 --seed 5 --strip-kinetics kinetics.bam half.bam 2>/dev/null
  $ samtools view half.bam | grep -cE '	(fi|fp|ri|rp|ip|pw):'
  0
  [1]
