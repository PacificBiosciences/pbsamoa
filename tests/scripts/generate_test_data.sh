#!/usr/bin/env sh
set -eu

OUTDIR="${1:-.}"

# 1. Coordinate-sorted BAM with diverse CIGAR, tags, and quality
cat > "${OUTDIR}/diverse.sam" << 'EOF'
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
@SQ	SN:chr2	LN:500
@RG	ID:rg1	SM:sample1	PL:ILLUMINA
@PG	ID:samtools	PN:samtools	VN:1.17
soft_clip	0	chr1	1	30	2S8M	*	0	0	NNACGTACGT	IIIIIIIIII	NM:i:0	RG:Z:rg1
insertion	0	chr1	10	30	5M2I3M	*	0	0	ACGTAXXCGT	IIIIIIIIII	NM:i:2
deletion	0	chr1	20	30	5M3D5M	*	0	0	ACGTAACGTA	IIIIIIIIII	NM:i:3
skip	0	chr1	50	30	3M100N3M	*	0	0	ACGACG	IIIIII	NM:i:0	XS:A:+
hard_clip	0	chr1	100	30	3H7M	*	0	0	ACGTACG	IIIIIII	NM:i:0
paired_r1	99	chr1	200	30	10M	=	300	110	ACGTACGTAC	IIIIIIIIII	NM:i:0	RG:Z:rg1
paired_r2	147	chr1	300	30	10M	=	200	-110	ACGTACGTAC	IIIIIIIIII	NM:i:1
unmapped	4	*	0	0	*	*	0	0	ACGT	IIII
secondary	256	chr1	400	0	5M	*	0	0	ACGTA	IIIII	NM:i:0
supplementary	2048	chr2	1	30	5M	*	0	0	ACGTA	IIIII	SA:Z:chr1,100,+,5M,30,0;
chr2_read	0	chr2	100	30	5M	*	0	0	ACGTA	IIIII	NM:i:0
no_qual	0	chr1	500	30	5M	*	0	0	ACGTA	*
array_tag	0	chr1	600	30	5M	*	0	0	ACGTA	IIIII	XA:B:i,1,2,3,4,5
EOF

samtools view -bS "${OUTDIR}/diverse.sam" | samtools sort -o "${OUTDIR}/diverse.bam"
samtools index "${OUTDIR}/diverse.bam"

# 2. Large file for performance testing (~10k records, 200bp reads)
python3 -c "
import random
random.seed(42)
print('@HD\tVN:1.6\tSO:coordinate')
print('@SQ\tSN:chr1\tLN:10000000')
bases = 'ACGT'
for i in range(10000):
    seq = ''.join(random.choice(bases) for _ in range(200))
    qual = 'I' * 200
    pos = i * 500 + 1
    nm = random.randint(0, 5)
    print(f'read{i:05d}\t0\tchr1\t{pos}\t30\t200M\t*\t0\t0\t{seq}\t{qual}\tNM:i:{nm}\tRG:Z:bench')
" | samtools view -bS - | samtools sort -o "${OUTDIR}/benchmark.bam"
samtools index "${OUTDIR}/benchmark.bam"

# 3. Header-only BAM
printf '@HD\tVN:1.6\n@SQ\tSN:empty\tLN:1000\n' | \
  samtools view -bS - > "${OUTDIR}/empty_records.bam"
samtools index "${OUTDIR}/empty_records.bam"

# 4. Unsorted BAM
cat > "${OUTDIR}/unsorted.sam" << 'EOF'
@HD	VN:1.6
@SQ	SN:chr1	LN:1000
read_b	0	chr1	500	30	5M	*	0	0	ACGTA	IIIII
read_a	0	chr1	100	30	5M	*	0	0	ACGTA	IIIII
read_c	0	chr1	300	30	5M	*	0	0	ACGTA	IIIII
EOF
samtools view -bS "${OUTDIR}/unsorted.sam" > "${OUTDIR}/unsorted.bam"

# 5. Unaligned HiFi BAM with a .pbi (300 ZMWs, two of them split by strand)
python3 -c "
import random
random.seed(7)
movie = 'm84000_240101_000000_s1'
print('@HD\tVN:1.6\tSO:unknown\tpb:5.0.0')
print(f'@RG\tID:a1b2c3d4\tPL:PACBIO\tDS:READTYPE=CCS;BINDINGKIT=102-739-100;SEQUENCINGKIT=102-118-800;BASECALLERVERSION=5.0;FRAMERATEHZ=100.000000\tPU:{movie}\tPM:REVIO')
zmws = sorted(random.sample(range(1, 50000000), 300))
for i, zmw in enumerate(zmws):
    suffixes = ['ccs/fwd', 'ccs/rev'] if i in (10, 200) else ['ccs']
    for suffix in suffixes:
        length = random.randint(200, 2000)
        seq = ''.join(random.choice('ACGT') for _ in range(length))
        qual = ''.join(chr(33 + random.randint(20, 40)) for _ in range(length))
        print(f'{movie}/{zmw}/{suffix}\t4\t*\t0\t255\t*\t*\t0\t0\t{seq}\t{qual}\tRG:Z:a1b2c3d4\tzm:i:{zmw}\trq:f:0.999')
" | samtools view -b --no-PG -o "${OUTDIR}/hifi.bam" -
python3 "$(dirname "$0")/make_pbi.py" "${OUTDIR}/hifi.bam"

echo "Test data generated in ${OUTDIR}"
