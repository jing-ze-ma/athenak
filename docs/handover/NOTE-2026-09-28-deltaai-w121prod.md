# NOTE 2026-09-28 (DeltaAI -> viper, Caltech): WASP-121b fresh start is ALSO queued on DeltaAI (user 09-28)

User 09-28: queue the WASP-121b 1x and 10x fresh start on DeltaAI too, alongside Caltech
(NOTE-2026-09-28-caltech-w121prod.md, jobs 3608646/3608647) and viper (NOTE-2026-09-29-viper-w121prod.md, jobs
12018387/12018389). All three copies run the same setup; the user decides which to keep.

## Jobs

- **1x 3257296, 10x 3257297**: partition ghx4, 1 node x 2 GH200 each, 24 h limit, submitted 2026-09-28 18:40 CDT.
  Slurm's start estimate is 09-29 ~22:13 CDT (it moved between 09-29 19:00 and 09-30 13:30).
- Expected speed (2 GH200, fresh-start dt): ~39 rot/h (1x), ~26 rot/h (10x); 300 rot in ~8 h / ~12 h.
- Once one of the other copies of an arm is running, the DeltaAI copy is stopped with
  `touch /u/jma20/ATHENAK/prod_claim/CANCEL_<1x|10x>` (the job exits within seconds of starting) or `scancel`.
  Tell the user, or push a NOTE here.

## Setup (same as Caltech)

- Inputs: `docs/handover/caltech-2026-09-26/inputs/w121prod_0928/w121prod_{1x,10x}.athinput` from 1832f493, unchanged
  (lhllc, 1x maxit 16, 10x maxit 24 + dtmax 0.25, xstep 8, flux_hst_floor, closed wall, top sponge on / bottom off,
  tlim 300 rot; hst+log 10/rot, bin 2 rot, rst 0.5 rot); ICs from `bench-2026-09-28/inputs`. MeshBlocks 16^2
  (24), as on the other machines.
- Binary: src 6d690e09 (the Caltech source), CUDA 12.9, Kokkos HOPPER90 + ARMV9_GRACE, cray-mpich 9.0.1 GPU-aware,
  md5 01a6cc93. Built by the incremental worktree script `build_inc_deltaai.sh` (on DeltaAI, not in the repo).
- Launch: the AthenaK docs' DeltaAI recipe (all GPUs visible, `KOKKOS_MAP_DEVICE_ID_BY=mpi_rank`,
  `SLURM_CPU_BIND=cores`). NOT `--gpus-per-task=1 --gpu-bind=closest`: its cgroups block CUDA IPC, and GPU-aware
  Cray MPICH then hangs, or runs 2.8-3.3x slower with IPC off.
- Each job first runs a 50-cycle smoke (ck_impl_verbose) and stops on FATAL or NOT-CONVERGED, as Caltech does.

## Checks done on DeltaAI before queueing

- GPU test suite at 6d690e09: 41/42 pass; the only failure is `test_nr_geneos_repro_gpu`, whose shipped input has a
  relative ck path and which asserts one .hst while the input now writes two (repo test still needs fixing). With
  absolute ck paths and all .hst compared: 3 runs bitwise identical.
- Smoke on the exact production inputs, both arms (jobs 3257349/3257350): rc 0, 0 FATAL, 0 NOT-CONVERGED, no NaN,
  rst/bin/log/hst written. `lhllc` really is active: dt at cycle 20 is 1x 15.54181 s vs 15.54189 s (hllc), the same
  pair viper found; 10x 16.20543 s vs 16.20547 s.
- The code at 6d690e09 is ~1 % slower per cycle than the bench binary 11c9a5be (1x 11.87 vs 11.71 ms, 10x 15.01 vs
  14.87 ms at 2 GPUs).

## Other DeltaAI results (09-28)

- Cross-cluster timing: branch `bench-results-deltaai` (0d804760), `bench-2026-09-28/RESULTS_deltaai.md`. 1x 18.19 /
  11.71 / 8.12, 10x 22.68 / 14.87 / 10.53 ms/cycle at 1 / 2 / 4 GH200; 1.4-1.7x faster than MI300A.
- MeshBlock split at 2 GPUs (bench inputs): per panel 1x1 (32^2, 6 MBs) is 6.9 % faster than 2x2 (16^2); 4x4 (8^2)
  is 23-28 % slower. Not bitwise across splits (window dt differs by 0.1-0.3 %). User: keep 16^2 for this spin-up;
  32^2 is an option for later 2-GPU runs that are not cross-machine checks. 4 GPUs not probed.
- Data: DeltaAI's copy of ckdata10 had cia/ray/sw_flux symlinks pointing to viper paths; repointed to the local
  exo_fms_ck tables. Other machines copying data dirs should check with `find <data> -xtype l`.
