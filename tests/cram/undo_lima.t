Undo lima demultiplexing: restore the pre-clip reads and drop the barcode bookkeeping.
The fixture is a lima 2.5.1 output (three CCS reads, `ls` undo stack present).

  $ "${PBSAMOA}" undo "${TESTDIR}"/../data/lima_demux.bam undo.bam
  pbsamoa undo: 3 records, 3 restored

Every read grows by its clipped leading + trailing barcode (bx):

  $ samtools view "${TESTDIR}"/../data/lima_demux.bam | awk '{print $1, length($10)}'
  m64013e_220721_210631/0/ccs 2483
  m64013e_220721_210631/1/ccs 2261
  m64013e_220721_210631/4/ccs 3324
  $ samtools view undo.bam | awk '{print $1, length($10), length($11)}'
  m64013e_220721_210631/0/ccs 2513 2513
  m64013e_220721_210631/1/ccs 2293 2293
  m64013e_220721_210631/4/ccs 3356 3356

The clipped flanks (bl/bt) are back at the read ends:

  $ samtools view "${TESTDIR}"/../data/lima_demux.bam | head -1 | grep -o 'b[lt]:Z:[ACGT]*'
  bl:Z:TGCATGTCATGAGTAT
  bt:Z:ATACTCATGACATG
  $ samtools view undo.bam | head -1 | cut -f10 | cut -c1-16
  TGCATGTCATGAGTAT
  $ samtools view undo.bam | head -1 | cut -f10 | rev | cut -c1-14 | rev
  ATACTCATGACATG

Demultiplexing tags, the CCS query window, and the barcode context bits are gone:

  $ samtools view undo.bam | grep -c -E '\t(bc|bq|bx|bl|bt|ql|qt|ls|qs|qe|cx):'
  0
  [1]
  $ samtools view undo.bam | cut -f12- | tr '\t' '\n' | grep '^RG:' | sort -u
  RG:Z:c8c82efd

Header: barcoded read groups collapse onto their base ID without Barcode* fields, the
lima @PG is dropped, and a pbsamoa.undo @PG is chained:

  $ samtools view -H undo.bam | grep '^@RG' | cut -f2
  ID:34327a9c
  ID:c8c82efd
  ID:d1c2e366
  $ samtools view -H undo.bam | grep -c 'Barcode'
  0
  [1]
  $ samtools view -H --no-PG undo.bam | grep '^@PG' | cut -f2
  ID:ccs
  ID:zmwfilter
  ID:pbsamoa.undo
