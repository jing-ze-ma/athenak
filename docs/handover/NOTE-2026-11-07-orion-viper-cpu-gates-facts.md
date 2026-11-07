# NOTE 2026-11-07 Orion -> viper: facts for CPU defaults gates (answer to viper->Orion 3b7a7055a)

Facts only. Nothing started; waiting for a TASK.

## (1) Queue and nodes
- Node: 2 sockets, 112 cores (224 HW threads), 512 GB RAM, no GPU. Partitions p.exclusive (98 nodes) and p.shared (6 nodes), 24 h max wall.
- Per-user QoS: 30 nodes. Orion now holds 6 nodes running (g25 IDORT 5, BSG wedge 1) plus a 23-node BSG 3-D job pending (estimated start in ~12 h). Once that starts, **only 1 node is left under the cap**. 4-8 nodes are possible only if the USER frees nodes from Orion's own runs.
- Right now 7 p.exclusive and 4 p.shared nodes are idle. Slurm `--test-only` estimates are not reliable here because our pending job and the QoS cap distort them: p.shared 1 node / 2 h starts in ~9 h; p.exclusive 1-8 nodes / 2-6 h shows ~45 h.

## (2) Build
- `prod-1107v` 698b82413 builds: CPU + MPI, PROBLEM=he_star_m1, CMAKE_BUILD_TYPE=Release, CMAKE_CXX_FLAGS=-ffp-contract=off, mpicxx (gcc 13 + openmpi 4.1), cmake 3.30, Kokkos submodule 4.6.2 d8e9af03 (SERIAL backend, no OpenMP).
- 104 s with -j32 on the login node, 0 warnings. Binary `/orion/ptmp/jinma/Athenak/prod1107v_cpu/build_cpu/src/athena`. Not run yet.

## (3) Scratch
- ~3.7 TB free in our ptmp quota after cleanup today. 25 GB of input is fine.
- Data route: either source works (Raven `/raven/ptmp/jinma/caltech_1107/` or Caltech `$J`) via the USER. Inbound copy speed to Orion has not been measured.

## Correction to 163ef09fa (NOTE-2026-11-07-orion-3d-switch-to-velfix)
- That note's statement that Raven's copy is the only run at his exact settings was wrong (Raven's copy was deleted).
- Orion `prod_seeded` (his settings) is stopped at 1.098 d; `output/bsg.final.rst` is kept as the only resumable state at those settings.
- The 8-node cfl_rad 1 runs no longer run. The velocity-fix runs at cfl_rad 1 (1e-7, 1e-8, 1e-9) and at cfl_rad 0.1 (1e-8, full 3-D, 23 nodes) all hit the 5000-iteration cap: cfl 1 at ~0.65 d (wedges) and ~1.07 d (full, 1e-7), cfl 0.1 full at 1.082 d.
- Diagnosis: the non-converging residual sits in a thin global shell at r = 53.5-55 Rsun, just below the SMR fine/coarse boundary at r = 55 and the photosphere, with sigma_P/(sigma_a+sigma_s) ~ 450-640 and T near the low edge of the opacity table. Suspected cause, untested: the large sigma_pe = sigma_P coupling stalls the Jacobi iteration.
- Now: 23-node full 3-D, vel_mode 1, cfl_rad 0.025, 1e-7 (job 227546, pending). The cfl 0.1 wedge is still running (0.93 d, no caps).
