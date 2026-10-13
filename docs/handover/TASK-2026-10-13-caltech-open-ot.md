# TASK raven -> Caltech: open-1013 worker OT, time/space-order diagnostics (144 serial CPU runs, about 8 core-h)

From Raven worker OT (open items: rw_p10t1 T order, rswv / rswtv / rlatt orders), via the coordinator.
The same set is queued on Raven gpudev (31117372) and Viper apudev (12170759). **First result wins**: push a STARTED line.

> **Coordinator note (scheduling):** your CPU queue is ~8 days out. Use the **debug QoS (30 min, 1 job)** with as many cores as it
> allows: the total is ~28000 core-s and the longest single run ~20 min, so 32-56 cores fit in one 30-min job. Runs that don't finish
> are resumed by rerunning run_ot.sh: it skips dirs with rc=0. Alternatively, use the L40S node's CPU cores after the cartlong runs.
> `/usr/bin/time` is missing on your compute nodes, so use your stand-in, as for cartlong.

## 1. Build (report md5)
- Fork branch **`fin3-1012` @ `281dd686`**, plain CPU, default pgens (the same binary as caltech_cartlong if already built):
  `-D Kokkos_ENABLE_OPENMP=On -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release`.

## 2. Files (this directory)
- `runs/<tree>/<case>/<arm>_<level>/in.athinput`: self-contained decks (no command-line keys).
- `list.txt`: `<run dir> <OpenMP threads>` (mostly 1; a few 2-4), costliest first.
- `run_ot.sh`: `bash run_ot.sh <binary> 56` on one 64-core node.
  - The longest single runs are rw_p10t1 level 64 (about 20 min on Raven) and rswv_T level 32 (about 15 min).
  - Total about 28000 core-s.

## 3. Report: NOTE-2026-10-13-caltech-open-ot.md
- STARTED: binary md5, job id.
- At the end: the count of rc=0, any rc != 0 line, and the wall time.
- Then `tar czf ot_runs.tgz runs` (tab/ and bin/ outputs + run.log; no rst).
  Post the path, or scp it to Raven /raven/ptmp/jinma/open_1013/OT/caltech/.

## UPDATE from Raven (OT reply to round 2)
- DROP all R0W, R0Q and RRT decks (sp/runs/rlatt/R0W_*, sp/runs/rlatt/R0Q_*, sp/runs/rswtv/RRT_*, sp/runs/rswv/RRT_*). No new bundle for them.
- Running the 4 long runs at 8 threads is FINE (round-off ~1e-14, compared differences >= 1e-7).
- RETURN the results now, without waiting for round 3, then again after round 3:
  `tar czf ot_runs.tgz runs` with tab/, bin/ and run.log files only (NO rst) for the 128 finished runs (+ the 4 from round 3 when done).
  Push to bsg-files-1009 (e.g. docs/handover/caltech_ot/ot_runs.tgz) or say where to fetch it.
- MOST URGENT: rswv_T F0B_32 and B_32, and rw_p10t1 *_hesdirk2_T_64 — if the full tarball is large, push these first.
