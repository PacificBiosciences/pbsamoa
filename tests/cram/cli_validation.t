Signed thread-count options reject undocumented negative values:

  $ cp "${TESTDIR}"/../data/spec_example.bam input.bam
  $ "${PBSAMOA}" zmi-index --threads -1 input.bam > /dev/null 2>&1
  [1]
  $ "${PBSAMOA}" dump --bgzf-threads -2 input.bam > /dev/null 2>&1
  [1]
  $ "${PBSAMOA}" dump --format-threads -1 input.bam > /dev/null 2>&1
  [1]
  $ "${PBSAMOA}" bench --bgzf-threads -2 input.bam > /dev/null 2>&1
  [1]
  $ "${PBSAMOA}" bench --decode-threads -1 input.bam > /dev/null 2>&1
  [1]

Chunk help does not advertise the unused built-in thread option:

  $ "${PBSAMOA}" chunk --help | awk '/--num-threads/ { count++ } END { print count + 0 }'
  0

Convert help distinguishes auto-detection from serial execution and shows the
effective records-per-slice default:

  $ "${PBSAMOA}" convert --help | grep -c -- '-1 = auto, 0 = serial'
  3
  $ "${PBSAMOA}" convert --help | grep -c '\[10000\]'
  1

Convert accepts the documented -1 auto sentinel it advertises in --help:

  $ "${PBSAMOA}" convert --bgzf-threads -1 input.bam accepts_auto.cram > /dev/null 2>&1
