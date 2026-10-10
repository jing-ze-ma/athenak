# NOTE Caltech -> Raven: capacity and turnaround (2026-10-10 ~07:00 PDT)

Caltech can take Raven's SHORT test tasks and turn them around fast. A standing poller checks bsg-files-1009 every 15 min
for `TASK-*-caltech-*.md` without a matching NOTE (see TASK-2026-10-10-caltech-test-runner). Each task gets a "picked up"
NOTE line when it starts, the full NOTE when it is done, and a row in docs/handover/caltech-runner-1010/LEDGER.md.

## Capacity right now
- **CPU (expansion): FAST.** About 5000 cores are idle, so batteries of ~100-2000 np-1 runs start immediately and take
  minutes to an hour. The plm gates rerun took 10 min of wall time and 3.7 CPU-h.
- **H200: SLOW today.** The partition is full. Our 2-H200 vgdfuse smoke (job 4316032) is pending on priority, with an
  estimated start of 10-13 23:20 (Slurm's estimate; backfill is often sooner). The vgdfuse-scaling jobs j1 and j2 wait
  behind it. Put GPU-timing tasks where an A100/MI300A is free, or accept a delay of a day or more here.
- Builds (CPU or CUDA, incremental) take about 10-20 min.

## In progress
- TASK-2026-10-10-caltech-beam-baseline: picked up. 90 CPU runs with athena_cpu_built_in_pgens_e50ec7d6 (f684f0c0).
- TASK-2026-10-10-caltech-vgdfuse-scaling: binaries ready; the jobs are queued behind the H200 queue (above).

## Not started (user hold)
No long production runs (BSG, He giant, AG Car, accretor) until the updated radiation scheme AND the user's word.
