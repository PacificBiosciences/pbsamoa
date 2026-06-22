Merge --bai builds the BAI on the fly. The streamed index must be byte-identical to
running bai-build on the finished output, and must query like samtools.

Split a coordinate-sorted file into two overlapping halves (alternating records) so the
merge takes the record-decoding heap path, not the disjoint-chain passthrough:

  $ samtools view -H "${TESTDIR}"/../data/many_records.bam > hdr.sam
  $ samtools view "${TESTDIR}"/../data/many_records.bam > body.sam
  $ awk 'NR%2==1' body.sam > a_body.sam
  $ awk 'NR%2==0' body.sam > b_body.sam
  $ cat hdr.sam a_body.sam | samtools sort -o a.bam 2>/dev/null
  $ cat hdr.sam b_body.sam | samtools sort -o b.bam 2>/dev/null

Merge with --bai; the sidecar index is written next to the output:

  $ "${PBSAMOA}" merge pb.bam a.bam b.bam --bai 2>/dev/null
  $ test -f pb.bam.bai && echo "bai exists"
  bai exists

The on-the-fly index is byte-identical to building it from the finished file:

  $ cp pb.bam.bai streamed.bai
  $ "${PBSAMOA}" bai-build pb.bam 2>/dev/null
  $ cmp streamed.bai pb.bam.bai && echo "identical"
  identical

The streamed index queries the same records samtools returns:

  $ "${PBSAMOA}" bai-query pb.bam chr1:1-1000 2>/dev/null | grep -v "^@" | sort > pb_q.txt
  $ samtools view pb.bam chr1:1-1000 | sort > st_q.txt
  $ diff pb_q.txt st_q.txt

--bai is rejected for an unsorted concat (no valid index exists), and no sidecar is left:

  $ "${PBSAMOA}" merge cc.bam a.bam b.bam --concat --bai 2>&1 | grep -c "bai"
  1
  $ "${PBSAMOA}" merge cc.bam a.bam b.bam --concat --bai > /dev/null 2>&1
  [1]
  $ test -e cc.bam.bai && echo "leaked" || echo "no sidecar"
  no sidecar

--bai is rejected for a queryname merge:

  $ "${PBSAMOA}" merge qn.bam a.bam b.bam --order queryname --bai > /dev/null 2>&1
  [1]
