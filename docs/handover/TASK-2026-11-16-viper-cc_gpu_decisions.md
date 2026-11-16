# TASK 2026-11-16 Orion -> viper (USER): clean scheme GPU design decisions (Raven A100, second task)

**Raven A100, second task, independent of stage 1; commit origin/clean-cc 32af6349e**
(push 10 = 4680ce500 + the wavefront-parallel FV ray sweep; outputs bitwise identical to 4680ce500,
where the references were generated, so they are valid). Package:
`docs/handover/task_1116_cc_gpu_decisions/`. The user wants the open design choices of the clean
scheme decided by GPU runs. Every arm is timed and scored for accuracy on the same
production-like decks. Run it after or alongside stage 1 (TASK-2026-11-16-viper-cc_gpu1-a100); it
needs nothing from it.

## Code and build
Fork `jing-ze-ma/athenak`, **origin/clean-cc = 32af6349e**: the same binary as stage 1
(`build_cuda_raven.sh` of task_1116_cc_gpu1, NMAX 432). It has the keys these arms need:
clean_line_solver (pcr | hybrid | thomas; default hybrid), clean_dissip_power (2 | 3; default 2),
FV-ray rotation for clean, and clean_rays = lld. The temporary keys are read only when they are
named in the deck. The clean FV ray sweep is wavefront-parallel on this commit: one launch per
dependency level over all cells x directions.

## Arms (runs.txt; every MeshBlock fits one A100-40GB: 1 block per GPU)
Decks:
- **env sp**: the rhse envelope, sp 548 x 64 x 64, 1 block, 1 GPU (~30 GB). Accuracy: max v/c_s
  (well-balanced, should stay ~1e-7).
- **env cs**: the rhse envelope on the cubed sphere, 548 x 64^2 per panel, 6 blocks, 6 GPUs.
  Accuracy as env sp.
- **thin1 / thin1l**: cs thin sphere, 64 x 8^2 and 64 x 16^2 per panel, 1 GPU (small).
  Accuracy: L1(E) against the steady exact solution.

All the clean arms use the production GL 10 x 12 direction set unless the arm says otherwise.

| group | arm (tag suffix) | keys | decks |
|---|---|---|---|
| G5L | clean_pcr / clean_hybrid / clean_thomas | clean_line_solver | env sp, env cs |
| G6 | base | GL 10x12, rotation off, A^2, fv1 rays | env sp, env cs, thin1, thin1l |
| G6 | a3 | clean_dissip_power = 3 | same |
| G6 | lld | clean_rays = lld | same |
| G6 | rot | vet_rotate_avg 4 + vet_rotate_advance + vet_gd_rotate_xfade, vet_gd_rotate_dt 6.5e-3 (envelope, ~5 cycles) / 2.8125 (thin, 5 cycles) | same |
| G6 | geo4 | vet_gd_geo = 4 (geodesic) | same |
| G6 | geo4rot | geo4 + rot | same |

Rotation with K = 4 direction sets means K FV sweeps per ray build (the union of the sets).
timing.py reports that cost as "FV ray sweep: ... s per build, ... passes per build"
(cln_fv time and launches per cln_gd_chi call; launches = passes x wavefront levels). Compare rot with base, and geo4rot with geo4.

## Run
```
export KOKKOS_TOOLS_LIBS=<kokkos-tools>/libkp_kernel_timer.so KP_READER=<kokkos-tools>/kp_reader
LAUNCH='srun -n %N --gpus-per-task=1 <your Raven A100 binding>' KEEPGOING=1 TIMEOUT=30m \
  ./run_all.sh <build_a100_n432>/src/athena G5L G6
python3 timing.py > timing.txt
# memory: base and rot of env sp and env cs again with nlim 6 and the space-time-stack tool
KOKKOS_TOOLS_LIBS=<kokkos-tools>/libkp_space_time_stack.so NLIM=6 RUNS=$PWD/runs_mem LAUNCH=... \
  ./run_all.sh <build_a100_n432>/src/athena G6
```
The envelope decks run nlim 25 (s/cycle over cycles 6-25); thin1/thin1l run nlim 100. Fresh
starts only. G5L/G6 are report-only (no gates), but send the log of any run that fails (FATAL,
NaN, rc != 0) or prints CLEAN WARNING.

## What to send back (one table per deck: rows = arms)
s/cycle, the initialisation time, the top kernels (timing.py), FV sweep cost and launches per ray
build, memory high-water mark per GPU, and the accuracy metric (max v/c_s or L1(E)). Also the
modules, the binary md5, and the run-dir paths.
cpu_values.txt holds Orion's CPU values of the thin-sphere arms (4680ce500) for comparison.

Reply: NOTE-2026-11-16-viper-orion-cc_gpu_decisions.md on fork/bsg-files-1009.
