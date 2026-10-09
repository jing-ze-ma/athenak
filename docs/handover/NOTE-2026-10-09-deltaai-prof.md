# NOTE DeltaAI: prof job CANCELLED (Caltech answered first)

Answer to TASK-2026-10-09-deltaai-prof (547b106f). Caltech pushed the results (NOTE-2026-10-09-caltech-prof,
7a20aeec) first, so DeltaAI job **3349371** (radprof, 1 node x 4 GH200, 30 min, ghx4-interactive) was scancelled
while still **PENDING** (queued ~2 h 15 min; never ran; 0 GPU-h). The user had chosen to keep it queued until a
first result appeared.

- Prepared and kept: `/work/nvme/bivj/jma20/prof_1009` (bundle, dry-run checked: all command-line keys exist in the
  inputs), simple-kernel-timer built in `bsg_hrdet_1009/kt/profiling/simple-kernel-timer`.
- Note for a re-run here: Caltech found that current kokkos-tools master prints the kernel table to stdout instead
  of writing `.dat` files (`KOKKOS_TOOLS_TIMER_BINARY=1` restores them); our checkout is the same master (1e98316).
- Interactive queue this afternoon: 30-min 1-node jobs waited 1-2 h+ (fairshare 0.093, down from 0.105 this
  morning), versus minutes this morning.
- Nothing queued on DeltaAI.
