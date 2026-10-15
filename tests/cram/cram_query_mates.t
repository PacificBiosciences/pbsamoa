Mates attached by samtools (NF), in full reads and in a CRAI query that stops before the mate:

  $ awk 'BEGIN{srand(7); s=""; for(i=0;i<2000;i++) s=s substr("ACGT",int(rand()*4)+1,1); print ">ref"; print s}' > ref.fa
  $ samtools faidx ref.fa
  $ seqat() { samtools faidx ref.fa ref:$1-$(($1+$2-1)) | tail -n +2 | tr -d '\n'; }
  $ Q=FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
  $ { printf '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:ref\tLN:2000\n'
  >   printf 'r1\t99\tref\t101\t60\t50M\t=\t1501\t1450\t%s\t%s\n' "$(seqat 101 50)" "$Q"
  >   printf 'r2\t99\tref\t201\t60\t50M\t=\t301\t150\t%s\t%s\n' "$(seqat 201 50)" "$Q"
  >   printf 'r2\t147\tref\t301\t60\t50M\t=\t201\t-150\t%s\t%s\n' "$(seqat 301 50)" "$Q"
  >   printf 'r1\t147\tref\t1501\t60\t50M\t=\t101\t-1450\t%s\t%s\n' "$(seqat 1501 50)" "$Q"; } > pairs.sam
  $ samtools view -C -T ref.fa -o pairs.cram pairs.sam
  $ samtools index pairs.cram
  $ "${PBSAMOA}" dump --reference ref.fa pairs.cram | grep -v "^@" | cut -f1-9 > pbsamoa_all.txt
  $ samtools view -T ref.fa pairs.cram | cut -f1-9 > samtools_all.txt
  $ diff pbsamoa_all.txt samtools_all.txt
  $ "${PBSAMOA}" dump --reference ref.fa --region ref:100-150 pairs.cram | grep -v "^@" | cut -f1-9
  r1\t99\tref\t101\t60\t50M\t=\t1501\t1450 (esc)
