Merge error handling: bad arguments and incompatible inputs must exit non-zero.

Build a coordinate-sorted input and a second input with a different reference:

  $ samtools sort -o a.bam "${TESTDIR}"/../data/unsorted.bam 2>/dev/null
  $ printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:other\tLN:5000\n' > other.sam
  $ printf 'r1\t0\tother\t10\t30\t4M\t*\t0\t0\tACGT\t*\n' >> other.sam
  $ samtools sort -o other.bam other.sam 2>/dev/null

Incompatible @SQ reference lists:

  $ "${PBSAMOA}" merge out.bam a.bam other.bam 2>&1 | grep -c "incompatible"
  1
  $ "${PBSAMOA}" merge out.bam a.bam other.bam > /dev/null 2>&1
  [1]

Sort-order mismatch (coordinate inputs, queryname merge requested):

  $ "${PBSAMOA}" merge --order queryname out.bam a.bam > /dev/null 2>&1
  [1]

Missing input file:

  $ "${PBSAMOA}" merge out.bam a.bam no_such_file.bam > /dev/null 2>&1
  [1]

No input files (only an output operand):

  $ "${PBSAMOA}" merge out.bam > /dev/null 2>&1
  [1]
