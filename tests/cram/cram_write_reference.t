Write reference-based CRAM with --reference: a mapped record stores only its differences to the
reference (substitutions, clips, an insertion, a deletion, a lowercase reference stretch, a base
past the reference end), and each @SQ line gets M5.

  $ printf ">ref\nACGTACGTACGTacgtacgtACGTNCGTAC\n>ref2\nTTTTGGGGCCCCAAAA\n" > ref.fa
  $ printf "@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:ref\tLN:30\n@SQ\tSN:ref2\tLN:16\n" > in.sam
  $ printf "match\t0\tref\t1\t60\t30M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTNCGTAC\t??????????????????????????????\n" >> in.sam
  $ printf "subst\t16\tref\t3\t60\t2S8M2I5M3D6M\t*\t0\t0\tTTGTACATACNNGTACGGRNCGT\t???????????????????????\n" >> in.sam
  $ printf "eqx\t0\tref\t11\t60\t4=1X5=\t*\t0\t0\tGTACTTACGT\t??????????\n" >> in.sam
  $ printf "pastend\t0\tref2\t13\t60\t8M\t*\t0\t0\tAAAAGGGG\t????????\n" >> in.sam
  $ printf "unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t????\n" >> in.sam
  $ "${PBSAMOA}" convert --reference ref.fa in.sam out.cram
  $ samtools view -H out.cram | grep '^@SQ'
  @SQ\tSN:ref\tLN:30\tM5:e6db79b744cd9475e4325b08ac87fabe (esc)
  @SQ\tSN:ref2\tLN:16\tM5:aca9e0e1019366b28c22c661b728ae6a (esc)

samtools reads the records back; it gives M for = and X:

  $ samtools view --input-fmt-option decode_md=0 -T ref.fa out.cram | cut -f 1-11
  match\t0\tref\t1\t60\t30M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTNCGTAC\t?????????????????????????????? (esc)
  subst\t16\tref\t3\t60\t2S8M2I5M3D6M\t*\t0\t0\tTTGTACATACNNGTACGGRNCGT\t??????????????????????? (esc)
  eqx\t0\tref\t11\t60\t10M\t*\t0\t0\tGTACTTACGT\t?????????? (esc)
  pastend\t0\tref2\t13\t60\t8M\t*\t0\t0\tAAAAGGGG\t???????? (esc)
  unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t???? (esc)

pbsamoa reads them back too:

  $ "${PBSAMOA}" dump --no-header --reference ref.fa out.cram | cut -f 1-11 > pbsamoa.sam
  $ samtools view --input-fmt-option decode_md=0 -T ref.fa out.cram | cut -f 1-11 > samtools.sam
  $ diff pbsamoa.sam samtools.sam

The reference-free output of the same input stores all bases:

  $ "${PBSAMOA}" convert in.sam free.cram
  $ samtools view free.cram | cut -f 1-11 | diff - samtools.sam
