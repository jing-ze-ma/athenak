# NOTE 2026-11-16 Orion -> viper (USER): clean scheme GPU stage 1 is coming

The user asked to start GPU tests of the clean cell-centred scheme now. A full TASK
(TASK-2026-11-16-viper-cc_gpu1.md + docs/handover/task_1116_cc_gpu1/) follows within a few hours,
once Orion has finished a CUDA compile check of the clean code (it has never been compiled for CUDA).

Code: fork branch clean-cc, head 18309ae7a (push 7). The TASK may name a newer commit with CUDA
portability fixes; Orion will push it to the fork first.

Planned content, so you can prepare an allocation and a CUDA build environment:
- Build: Kokkos CUDA (Ampere80; Hopper90 if you use H200), GPU-aware MPI, -ffp-contract=off,
  device --fmad=false, built-in pgens; two binaries with Athena_M1_VGD_NMAX 432 and 4096.
- G1 CPU parity on 1 GPU (small Cartesian IDORT decks + one linear wave), G2 4 GPUs / 4 blocks,
  G3 spherical-polar wedges (well-balanced envelope, thin pulse, 3-D with lateral blocks),
  G4 cubed sphere (6 and 24 blocks), G5 timing of one production-size block, clean vs production,
  with a kernel breakdown and the device memory high-water mark.
- Gates: scores against Orion CPU references (shipped in the package), 0 CLEAN WARNING, 0 dt_retry,
  energy / momentum / L_z round-off; G5 report only.

Known: the clean FV ray sweep is not yet GPU-parallel (one serial thread per block and direction),
so G5 will be slow; that is part of what we want to measure. No action needed until the TASK arrives.
