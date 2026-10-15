Subsample keeps the same records as samtools view --subsample.

  $ BAM="${TESTDIR}"/../data/benchmark.bam
  $ parity() {
  >   "${PBSAMOA}" subsample --fraction "$1" --seed "$2" "$3" out.bam 2>/dev/null &&
  >   samtools view out.bam > pbsamoa.rec &&
  >   samtools view --subsample "$1" --subsample-seed "$2" "$3" > samtools.rec &&
  >   cmp pbsamoa.rec samtools.rec && echo "same $(wc -l < pbsamoa.rec | tr -d ' ')"
  > }

Seed 0 is not scrambled with rand(), so these counts are the same on every platform:

  $ "${PBSAMOA}" subsample --fraction 0.1 --seed 0 "${BAM}" out.bam
  pbsamoa subsample: kept 993 of 10000 records, 198600 of 2000000 bases (fraction 0.1, seed 0)
  $ parity 0.1 0 "${BAM}"
  same 993
  $ parity 0.5 0 "${BAM}"
  same 5032
  $ parity 0.5 0 "${TESTDIR}"/../data/many_records.bam
  same 238

A non-zero seed goes through the C library rand(), as in samtools:

  $ parity 0.25 7 "${BAM}"
  same * (glob)
  $ parity 0.25 4294967295 "${TESTDIR}"/../data/diverse.bam
  same * (glob)

The auto seed is the FNV-1a hash of the input header text. The @CO line records it
like samtools >= 1.24 does:

  $ "${PBSAMOA}" subsample --fraction 0.3 "${BAM}" out.bam 2>/dev/null
  $ samtools view -H out.bam | grep '^@CO'
  @CO\tSub-sampled fraction=0.3 seed=106762583 (esc)
  $ samtools view -H out.bam | grep -c '^@PG.*ID:pbsamoa.subsample'
  1
  $ samtools view out.bam > pbsamoa.rec
  $ samtools view --subsample 0.3 --subsample-seed 106762583 "${BAM}" > samtools.rec
  $ cmp pbsamoa.rec samtools.rec

The hash uses the header text as stored, not a canonical form of it:

  $ { printf '@CO\tfirst\n'; samtools view -h --no-PG "${TESTDIR}"/../data/many_records.bam; } \
  >   | samtools view -b --no-PG -o odd.bam -
  $ "${PBSAMOA}" subsample --fraction 0.3 odd.bam odd.out.bam 2>/dev/null
  $ samtools view -H odd.out.bam | grep '^@CO	Sub'
  @CO\tSub-sampled fraction=0.3 seed=* (esc) (glob)
  $ seed=$(samtools view -H odd.out.bam | sed -n 's/^@CO.*seed=//p')
  $ samtools view odd.out.bam > pbsamoa.rec
  $ samtools view --subsample 0.3 --subsample-seed "${seed}" odd.bam > samtools.rec
  $ cmp pbsamoa.rec samtools.rec
  $ test "${seed}" != 3496665649 && echo "differs from the canonical header seed"
  differs from the canonical header seed

Records with the same name stay together:

  $ "${PBSAMOA}" subsample --fraction 0.5 --seed 0 "${TESTDIR}"/../data/diverse.bam out.bam 2>/dev/null
  $ samtools view out.bam | cut -f1 | sort | uniq -c | awk '{print $1, $2}' > kept
  $ samtools view "${TESTDIR}"/../data/diverse.bam | cut -f1 | sort | uniq -c | awk '{print $1, $2}' > all
  $ comm -23 kept all

A fraction of 1 keeps every record:

  $ "${PBSAMOA}" subsample --fraction 1 "${BAM}" out.bam
  pbsamoa subsample: kept 10000 of 10000 records, 2000000 of 2000000 bases (fraction 1, seed 106762583)

Invalid arguments:

  $ "${PBSAMOA}" subsample "${BAM}" out.bam 2>&1 | grep -o -- '--fraction must be in (0, 1\]: 0'
  --fraction must be in (0, 1]: 0
  $ "${PBSAMOA}" subsample --fraction 0.5 --yield 1G "${BAM}" out.bam 2>&1 | grep -o 'give --fraction or --yield, not both'
  give --fraction or --yield, not both
  $ "${PBSAMOA}" subsample --fraction 1.5 "${BAM}" out.bam 2>&1 | grep -o -- '--fraction must be in (0, 1\]: 1.5'
  --fraction must be in (0, 1]: 1.5
  $ "${PBSAMOA}" subsample --fraction 0.5 --seed abc "${BAM}" out.bam 2>&1 | grep -o 'invalid --seed: abc'
  invalid --seed: abc
  $ "${PBSAMOA}" subsample --fraction 0.5 "${BAM}" 2>&1 | grep -o 'subsample requires exactly two arguments.*'
  subsample requires exactly two arguments: <input.bam> <output.bam>
