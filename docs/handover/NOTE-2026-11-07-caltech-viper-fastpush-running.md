# NOTE 2026-11-07 Caltech -> viper: fastpush running (12:10 PST)

- USER ran `go.sh 3` at ~12:10 PST: old copies stopped by exact PID (fast.sh, its pcp + 8 HG streams, the RSG rsync + ssh); fastpush.sh 3 running in tmux on login3 with 4 connections (existing master + cmk2/3/4) x 4 = 16 streams.
- Order: RSG -> `raven_gap_1105/rsgE/rst/rsg.01295.rst.part` (re-sent from scratch), then HG -> `raven_gap_1105/hg/rst/hegiant.00139.rst.part` (byte ranges overwrite the partial .part).
- First reading: RSG .part 0.61 GB allocated after ~30 s (~20 MB/s, early; not yet a measured rate). The old rsync partial is on Raven as `rsg.01295.rst.TRUNCATED_rsync_partial` (23.7 GB): yours to delete after the .part checks out.
- Small states hgb/rsb: you already checked and renamed them (thanks).
- Caltech will post the per-file MB/s lines that fastpush prints when RSG and HG complete.
