CRAI-backed region queries in pbsamoa dump:

  $ cat > query_region.sam <<'EOF'
  > @HD	VN:1.6	SO:unknown
  > @SQ	SN:chr1	LN:1000
  > region_mapped_1	0	chr1	11	60	5M	*	0	0	ACGTA	IIIII
  > region_mapped_2	0	chr1	19	60	5M	*	0	0	TTGCA	IIIII
  > region_mapped_3	0	chr1	41	60	5M	*	0	0	GGGGG	IIIII
  > region_unmapped_1	4	*	0	0	*	*	0	0	CCCCC	IIIII
  > EOF
  $ samtools view -bS query_region.sam > query_region.bam
  $ "${PBSAMOA}" convert --write-crai query_region.bam query_region.cram
  $ "${PBSAMOA}" dump --region chr1:11-20 query_region.cram | grep -v "^@" | cut -f1
  region_mapped_1
  region_mapped_2
  $ "${PBSAMOA}" dump --index query_region.cram.crai --region '*' query_region.cram | grep -v "^@" | cut -f1
  region_unmapped_1
  $ "${PBSAMOA}" dump --region chr1:1-5 query_region.cram | awk 'BEGIN{n=0} !/^@/{++n} END{print n}'
  0
