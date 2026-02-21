Read reference-based CRAM using explicit reference FASTA:

  $ printf ">ref\nACGTACGTACGTACGTACGT\n" > ref.fa
  $ printf "@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:ref\tLN:20\nread1\t0\tref\t1\t60\t10M\t*\t0\t0\tACGTTCGTAC\tFFFFFFFFFF\n" > ref.sam
  $ samtools view -bS ref.sam > ref.bam
  $ samtools view -T ref.fa --output-fmt cram,use_rans=0,use_bzip2=0,use_lzma=0,use_fqz=0,use_tok=0,use_arith=0 -o ref.cram ref.bam
  $ "${PBSAMOA}" dump --reference ref.fa ref.cram > pbsamoa_ref.sam
  $ grep -v "^@PG" pbsamoa_ref.sam > pbsamoa_ref_nopg.sam
  $ sed -E 's/\tMD:Z:[^\t]*//g; s/\tNM:i:[^\t]*//g' pbsamoa_ref_nopg.sam > pbsamoa_ref_norm.sam
  $ samtools view -h --no-PG -T ref.fa ref.cram > samtools_ref.sam
  $ grep -v "^@PG" samtools_ref.sam > samtools_ref_nopg.sam
  $ sed -E 's/\tMD:Z:[^\t]*//g; s/\tNM:i:[^\t]*//g' samtools_ref_nopg.sam > samtools_ref_norm.sam
  $ diff pbsamoa_ref_norm.sam samtools_ref_norm.sam

Read reference-based CRAM with wrong FASTA fails with MD5 mismatch:

  $ printf ">ref\nTTTTTTTTTTTTTTTTTTTT\n" > wrong.fa
  $ { "${PBSAMOA}" dump --reference wrong.fa ref.cram > pbsamoa_wrong_ref.sam 2> pbsamoa_wrong_ref.err; dump_status=$?; [ "${dump_status}" -ne 0 ]; } 2> /dev/null
  $ grep "reference MD5 mismatch" pbsamoa_wrong_ref.err
  *reference MD5 mismatch* (glob)

Read reference-free CRAM without providing a FASTA:

  $ samtools view --output-fmt cram,no_ref,use_rans=0,use_bzip2=0,use_lzma=0,use_fqz=0,use_tok=0,use_arith=0 -o ref_no_ref.cram ref.bam
  $ "${PBSAMOA}" dump ref_no_ref.cram > pbsamoa_ref_no_ref.sam
  $ grep -v "^@PG" pbsamoa_ref_no_ref.sam > pbsamoa_ref_no_ref_nopg.sam
  $ awk 'BEGIN{OFS="\t"} /^@/{print;next} {$7="*";$8=0;$9=0; print}' pbsamoa_ref_no_ref_nopg.sam > pbsamoa_ref_no_ref_norm.sam
  $ samtools view -h --no-PG ref_no_ref.cram > samtools_ref_no_ref.sam
  $ grep -v "^@PG" samtools_ref_no_ref.sam > samtools_ref_no_ref_nopg.sam
  $ awk 'BEGIN{OFS="\t"} /^@/{print;next} {$7="*";$8=0;$9=0; print}' samtools_ref_no_ref_nopg.sam > samtools_ref_no_ref_norm.sam
  $ diff pbsamoa_ref_no_ref_norm.sam samtools_ref_no_ref_norm.sam
