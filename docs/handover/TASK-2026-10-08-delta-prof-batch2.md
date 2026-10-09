# TASK viper -> Delta: prof batch2 = RE-ISSUE of prof batch1 (user 10-09: queue it on Delta after all)

Run exactly TASK-2026-10-08-delta-prof-batch1.md (its WITHDRAWN header is lifted by this re-issue) with the scripts in
docs/handover/delta-prof-1009/ (prof.sbatch, prof_cuda.sh, prof_group.py, labels_98835d99.tsv): one 1-node x 4 A100 job,
30 min limit (<= 4 GPU-h), AG Car B 60 cycles + BSG reduced 4x1 30 cycles, kernel timer + implicit_timers arms.
- DeltaAI has the same task (TASK-2026-10-09-deltaai-prof on bsg-files-1009): first result wins. If
  NOTE-2026-10-09-deltaai-prof.md appears on fork/bsg-files-1009 before your job starts, cancel it (0 GPU-h) and say so.
- NOTE back: NOTE-2026-10-09-delta-prof-batch2.md on this branch (same content as batch1 asked) + LEDGER line "prof2".
