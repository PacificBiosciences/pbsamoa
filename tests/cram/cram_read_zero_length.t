pbsamoa decodes zero-length reads (SEQ=* and QUAL=*) from CRAM.

A zero-length run must not resolve an omitted BA or QS external block.

  $ printf '@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:ref\tLN:1000\n' > zl.sam
  $ printf 'zerolen\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n' >> zl.sam
  $ printf 'withseq\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n' >> zl.sam

Round-trip through a pbsamoa-written CRAM file:

  $ "${PBSAMOA}" convert zl.sam zl.cram
  $ "${PBSAMOA}" dump zl.cram > from_cram.sam
  $ diff zl.sam from_cram.sam

Use awk for portable tab handling and require exactly 11 fields:

  $ awk -F'\t' 'NF == 11 && $1 == "zerolen" && $2 == 4 && $3 == "*" && $4 == 0 && $5 == 0 && $6 == "*" && $7 == "*" && $8 == 0 && $9 == 0 && $10 == "*" && $11 == "*"' from_cram.sam | wc -l | tr -d ' '
  1

An all-zero-length slice omits the blocks:

  $ printf '@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:ref\tLN:1000\n' > allzl.sam
  $ printf 'a\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n' >> allzl.sam
  $ printf 'b\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n' >> allzl.sam
  $ "${PBSAMOA}" convert allzl.sam allzl.cram
  $ "${PBSAMOA}" dump allzl.cram | grep -cv '^@'
  2

The samtools-written form exercises the missing-block regression directly:

  $ samtools view -C -o allzl_st.cram allzl.sam 2> /dev/null
  $ "${PBSAMOA}" dump allzl_st.cram | grep -cv '^@'
  2
