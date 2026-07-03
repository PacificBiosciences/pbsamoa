Scatter (chaotic-deterministic) chunking via `pbsamoa chunk --mode scatter`.

Build an indexed BAM with 12 ZMWs (holes 0..11) in ascending file order:

  $ printf '@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:ref\tLN:100000\n' > input.sam
  $ z=0; while [ $z -lt 12 ]; do printf 'movie/%s/0_100\t0\tref\t%s\t30\t10M\t*\t0\t0\tACGTACGTAC\t*\n' "$z" "$((z * 100 + 1))" >> input.sam; z=$((z + 1)); done
  $ samtools view -bS input.sam > input.bam 2>/dev/null
  $ "${PBSAMOA}" zmi-build input.bam indexed.bam 2>/dev/null

Contiguous chunk 1 of 4 is the first quarter of the file (holes 0,1,2):

  $ "${PBSAMOA}" chunk indexed.bam 1 4 | grep -v '^@' | sed -E 's#movie/([0-9]+)/.*#\1#' | paste -sd' ' -
  0 1 2

Scatter chunk 1 of 4 (seed 7) samples across the whole file instead:

  $ "${PBSAMOA}" chunk indexed.bam 1 4 --mode scatter --tile 1 --seed 7 | grep -v '^@' | sed -E 's#movie/([0-9]+)/.*#\1#' | paste -sd' ' -
  2 4 9

Running all 4 scatter chunks visits every record exactly once (exact partition):

  $ for c in 1 2 3 4; do "${PBSAMOA}" chunk indexed.bam $c 4 --mode scatter --tile 1 --seed 7 | grep -v '^@'; done | sed -E 's#movie/([0-9]+)/.*#\1#' | sort -n | paste -sd' ' -
  0 1 2 3 4 5 6 7 8 9 10 11

Same (chunk, total, tile, seed) is deterministic across runs:

  $ "${PBSAMOA}" chunk indexed.bam 2 4 --mode scatter --tile 1 --seed 7 > run1.txt
  $ "${PBSAMOA}" chunk indexed.bam 2 4 --mode scatter --tile 1 --seed 7 > run2.txt
  $ diff run1.txt run2.txt && echo SAME
  SAME

--tile reads M consecutive ZMWs per seek; each chunk gets a contiguous run:

  $ "${PBSAMOA}" chunk indexed.bam 4 4 --mode scatter --tile 3 --seed 7 | grep -v '^@' | sed -E 's#movie/([0-9]+)/.*#\1#' | paste -sd' ' -
  0 1 2

--tile/--seed without --mode scatter is rejected:

  $ "${PBSAMOA}" chunk indexed.bam 1 4 --tile 3 2>&1
  *pbsamoa chunk ERROR: --tile/--seed require --mode scatter (glob)
  [1]
