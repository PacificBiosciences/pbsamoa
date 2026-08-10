This test checks a regression in the parallel BGZF read path. The path must not deadlock
when one decompressed batch produces more records than the output queue can hold.

The consumer thread pushed records into a fixed-capacity SPSC queue (OUTPUT_QUEUE_CAPACITY
in src/Bgzf.cpp). The thread signalled the reader only after it finished a whole batch.
When the queue filled, the consumer thread spun forever on the full queue. At the same
time, the reader thread waited on an empty-queue condition variable that no thread would
ever notify. Both threads stalled permanently.

The fixture must contain more records than one output queue in a batch; the command below
pins that requirement:

  $ samtools view -c "${TESTDIR}"/../data/benchmark.bam
  10000

With default options, which use parallel BGZF, the tool must terminate and match samtools:

  $ "${PBSAMOA}" dump "${TESTDIR}"/../data/benchmark.bam | grep -v "^@PG" > parallel.sam
  $ samtools view -h "${TESTDIR}"/../data/benchmark.bam | grep -v "^@PG" > samtools.sam
  $ diff parallel.sam samtools.sam

An explicit worker count follows the same code path. This test checks that path directly,
because the default worker count can change:

  $ "${PBSAMOA}" dump --bgzf-threads 4 "${TESTDIR}"/../data/benchmark.bam | grep -v "^@PG" > parallel4.sam
  $ diff parallel4.sam samtools.sam

The parallel path produces output that is byte-identical to the serial path, so
backpressure never reorders or drops records:

  $ "${PBSAMOA}" dump --bgzf-threads 0 "${TESTDIR}"/../data/benchmark.bam | grep -v "^@PG" > serial.sam
  $ diff serial.sam parallel.sam
  $ grep -vc "^@" parallel.sam
  10000
