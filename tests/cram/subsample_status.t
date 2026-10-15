The live status stream of subsample (--status-json): one JSON object per line, as in the
tool status stream of pbhuman. The fixture has 302 records in 300 ZMWs.

  $ HIFI="${TESTDIR}"/../data/hifi.bam
  $ cp "${HIFI}" nopbi.bam
  $ lines() { sed -E 's/"(t|wall_ms|started_wall_ms|pid|cpu_s|wall_s|peak_rss_bytes)":[0-9.]+/"\1":N/g; s/"(hostname|version|commandline|path)":"[^"]*"/"\1":S/g' "$1"; }

With a .pbi, the yield mode selects the ZMWs from the index (no progress), then writes.
The write phase counts reads, with the total from the .pbi:

  $ "${PBSAMOA}" subsample --yield 100K --seed 5 --status-json status/pbi.jsonl "${HIFI}" pbi.bam 2>/dev/null
  $ lines status/pbi.jsonl
  {"command":"subsample","commandline":S,"hostname":S,"num_threads":*,"phases":["select","write"],"pid":N,"schema_version":1,"started_wall_ms":N,"tool":"pbsamoa","type":"run_info","version":S} (glob)
  {"name":"select","state":"start","t":N,"type":"phase","wall_ms":N}
  {"metrics":{"bases":321748,"records":302,"target_bases":100000},"name":"select","state":"end","t":N,"type":"phase","wall_ms":N}
  {"name":"write","state":"start","t":N,"type":"phase","wall_ms":N}
  {"done":302,"phase":"write","t":N,"total":302,"type":"progress","unit":"reads","wall_ms":N}
  {"metrics":{"bases":321748,"kept_bases":100536,"kept_records":91,"records":302},"name":"write","state":"end","t":N,"type":"phase","wall_ms":N}
  {"metrics":{"bases":321748,"kept_bases":100536,"kept_records":91,"records":302,"seed":5,"target_bases":100000},"outputs":[{"kind":"bam","path":S}],"phase":"FINISHED","resources":{"cpu_s":N,"num_threads":*,"peak_rss_bytes":N,"wall_s":N},"t":N,"type":"final","wall_ms":N} (glob)

Without a .pbi, a count phase reads the BAM first. Its total is not known:

  $ "${PBSAMOA}" subsample --yield 100K --seed 5 --status-json nopbi.jsonl nopbi.bam out.bam 2>/dev/null
  $ lines nopbi.jsonl | grep -o '"phases":[^]]*]'
  "phases":["count","write"]
  $ lines nopbi.jsonl | grep '"type":"progress"'
  {"done":302,"phase":"count","t":N,"total":302,"type":"progress","unit":"reads","wall_ms":N}
  {"done":302,"phase":"write","t":N,"total":302,"type":"progress","unit":"reads","wall_ms":N}

The fraction mode has only the write phase:

  $ "${PBSAMOA}" subsample --fraction 0.5 --seed 5 --status-json fraction.jsonl "${HIFI}" half.bam 2>/dev/null
  $ lines fraction.jsonl | grep -o '"phases":[^]]*]'
  "phases":["write"]
  $ lines fraction.jsonl | tail -1 | grep -o '"phase":"[A-Z]*"'
  "phase":"FINISHED"

A failure ends the stream with a FAILED line and its cause. The phase that failed has no
end line:

  $ samtools view -h "${HIFI}" | awk 'NR != 10' | samtools view -b --no-PG -o stale.bam -
  $ cp "${HIFI}".pbi stale.bam.pbi
  $ "${PBSAMOA}" subsample --yield 100K --status-json stale.jsonl stale.bam stale.out.bam 2>/dev/null
  [1]
  $ lines stale.jsonl | grep -o '"name":"[a-z]*","state":"[a-z]*"'
  "name":"select","state":"start"
  "name":"select","state":"end"
  "name":"write","state":"start"
  $ lines stale.jsonl | tail -1 | grep -o '"cause":"[^"]*"\|"phase":"[A-Z]*"'
  "cause":"[pbsamoa] the .pbi does not match BAM record 6; rebuild it with pbindex"
  "phase":"FAILED"

Only the runs with --status-json wrote a stream:

  $ ls *.jsonl status
  fraction.jsonl
  nopbi.jsonl
  stale.jsonl
  
  status:
  pbi.jsonl
