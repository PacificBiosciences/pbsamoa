Forced spill under a constrained file-descriptor budget: a tiny --memory budget
forces an external merge over many run files, and a low open-file limit forces a
bounded multi-pass merge (the run count exceeds the fds available to open them at
once). The merged output must still equal the whole-file coordinate sort.

  $ ( ulimit -n 128; "${PBSAMOA}" sort --memory 1K "${TESTDIR}"/../data/many_records.bam spill.bam ) 2>/dev/null
  $ "${PBSAMOA}" sort "${TESTDIR}"/../data/many_records.bam mem.bam 2>/dev/null

Spilled and in-memory sorts produce identical record bodies (determinism):

  $ samtools view spill.bam > spill.sam
  $ samtools view mem.bam > mem.sam
  $ diff spill.sam mem.sam

And both match the samtools oracle:

  $ samtools sort -o st.bam "${TESTDIR}"/../data/many_records.bam 2>/dev/null
  $ samtools view st.bam > st.sam
  $ diff spill.sam st.sam

No temporary run files are left behind:

  $ ls spill.bam.run*.bam 2>/dev/null | wc -l | tr -d ' '
  0
