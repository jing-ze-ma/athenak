# NOTE 2026-10-08 (Delta): AthenaK bring-up on NCSA Delta -- results

Answer to TASK-2026-10-08-delta-bringup. Report only; nothing longer than the smokes was run.
Scripts: `docs/handover/delta-2026-10-08/` (build_delta.sh, smoke_delta.sh, compare_smoke.py).
Delta work dir: `/work/nvme/bivj/jma20/delta_1008/`.

**Summary:** the build works, smoke A passes (~1e-13 of viper and Raven), smoke B matches viper to ~1e-11 (10x
looser than the ~1e-12 band; Raven B reference asked for in NOTE-2026-10-08-delta-raven-b-question).
CPU test `test_hydro_fofc_cs_cpu` FAILS with the Cray compiler wrapper because of FMA contraction (traced below).

## 1. Machine
- Delta now runs the HPE Cray PE (RHEL 9.6). Login node: 128 cores, 250 GB.
- Account `bivj-delta-gpu` (project "performance-portabl..."): **100 GPU-hours deposited, 99 left**. No CPU
  account, so CPU work goes on the login node or a charged GPU node.
- Partitions (sinfo 15:00 CDT):
  - gpuA100x4 (gpua002-100: 4 x A100-SXM4-40GB, 64 cores, 250 GB; max 2 days; billing GPU=1000, i.e. about
    1 SU per A100-hour).
  - gpuA100x4-interactive (max 1 h, 1 job, billing 2x).
  - gpuA100x4-preempt; gpuA100x8 (6 nodes); gpuA40x4 (98); gpuH200x8 (8); gpuMI100x8 (1); cpu (136).
  - Queue: about 1980 jobs pending on gpuA100x4. Of those, about 750 are waiting on Dependency, about 240 are
    JobHeldUser, and about 920 are waiting on Priority.
  - **The smoke job (20 min, exclusive node) waited 19 s.** Short jobs backfill fast. Long jobs were not tested.
- Filesystems (1 TB quota each): `/u/jma20` (100 GB), `/projects/bivj`, `/work/nvme/bivj`, `/work/hdd/bivj`.
  - **/work is shared with DeltaAI**: `/work/nvme/bivj/jma20` already holds DeltaAI's runs (hegiant_1008, bsg_*,
    ...). The Delta session keeps to `delta_1008/` there.
  - Purge policy: not checked.

## 2. Build (build_delta.sh, incremental: persistent worktree + build dir per target)
- Modules: the defaults (PrgEnv-gnu/8.7.0 with gcc-native/14 = gcc 14.2.1, cray-mpich/9.1.0,
  craype-x86-milan, craype-accel-nvidia80), then `module unload cudatoolkit; module load cuda/12.9`.
  - The C++ compiler is `kokkos/bin/nvcc_wrapper` with `NVCC_WRAPPER_DEFAULT_COMPILER=CC` (the Cray wrapper,
    which brings in MPI and GTL).
  - **Why not the default cudatoolkit/13.2:** in-tree Kokkos 4.6.2 does not compile with CUDA 13.
    `Kokkos_Cuda_Instance.hpp` calls `cudaMemAdvise(..., int device)`, and CUDA 13 removed that overload.
  - Side effect: cray-mpich's GTL (`libmpi_gtl_cuda`) is linked against libcudart.so.13, so the binary loads
    libcudart 12 and 13 side by side. This worked in the smokes (GPU-aware MPI, 4 ranks).
- he_gpu: `PROBLEM=he_star_m1`, MPI, CUDA, AMPERE80 + ZEN3.
  - Commit 081dd60b. Full build 845 s on the login node with -j32.
  - `bin/athena_he_gpu_081dd60b`, md5 `c605f0e83e9ba7fa253c340caddcfd01`.

## 3. CPU tests (login node, commit 081dd60b)
- `rad_m1/test_rad_m1_slab_cpu.py`: **SKIPPED**. Four of its five data files (ic_m1_V3edd_pgen.txt,
  arad_V3edd.txt, rosseland_he_x0.0_z0.02.txt, planck_he_x0.0_z0.02_ferg+tops.txt) are not on Delta.
- `hydro/test_hydro_fofc_cs_cpu.py`: **FAILS** with the default compiler (CXX=CC, the Cray wrapper):
  `FOFC left 309 energy-floor events on the cubed sphere (control run: 232917)`.
  - It also fails at 67ab56ba, so it does not come from 148a7859.
- **Trace:** the Cray wrapper silently adds `-march=znver3` (craype-x86-milan), which turns on FMA contraction.
  - Same CC + `CXXFLAGS=-ffp-contract=off`: **0** energy floors in the fofc=true run, so it **passes**
    (control 239660).
  - Plain `g++` (gcc-toolset-14, no -march): identical to that, also passes.
  - The FMA build gives 264 efloor events in cycles 0-5. They appear in both the fofc=true and fofc=false runs
    (FOFC count 0 at that point). It gives 45 more at cycle 46.
  - Total mass at t=0.02 differs by 0.2 % between the FMA and no-FMA builds; both runs floor heavily.
- **Reading:** the floor decisions in this near-vacuum rarefaction are knife-edge, and the test asserts exactly 0.
  - FOFC's floor-TEST pass reproduces GnomonicEquiangleRaiseVel's arithmetic on the trial state. With
    contraction on, the compiler can fuse the two copies differently, so the prediction and the real pass can
    disagree.
  - nvcc contracts FMAs by default on device code, so GPU cs runs may see the same thing. Not tested.
  - AG Car runs on the sp wedge, not the cubed sphere, so it is unaffected.
- Options, for the user / viper to decide:
  - Build CPU tests on Cray with `-ffp-contract=off`, or with plain g++.
  - Make the test tolerate a small residual floor count.
  - Make the floor-test and RaiseVel share one non-inlined routine, or set `#pragma STDC FP_CONTRACT OFF` in
    gnomonic_raisevel.hpp.
- No code was changed.

## 4. GPU smoke, AG Car (1 node, 4 x A100-40GB, 4 ranks, one 480x64x64 block per GPU; job 22759294, gpua092)
Settings: `KOKKOS_MAP_DEVICE_ID_BY=mpi_rank`, `MPICH_GPU_SUPPORT_ENABLED=1`, `srun --cpu-bind=cores`,
`time/nlim=10`. The agcar-files-1008 MD5SUMS check passed.

| | A Delta | rel. viper 12131991 | rel. Raven 31006215 | B Delta | rel. viper 12130898 |
|---|---|---|---|---|---|
| t | 1.2525591814774953e+04 | 1.3e-13 | 1.4e-13 | 2.4186988198547679e+03 | 1.9e-11 |
| mass | 6.1192879615057494e+32 | 4.9e-13 | 5.5e-13 | 1.7891809043987407e+31 | 5.8e-12 |
| tot-E | 1.6552578844686083e+47 | 4.2e-13 | – | 1.0633773260178835e+46 | 4.8e-12 |
| IC T-check / he_ic_balance | 7.64e-14 / 3.23845e-07 | | | 5.88e-14 / 1.07898e-07 | |
| Picard mean / max | 13.4 / 17 (viper 13.5/16) | | | 6.9 / 24 (= viper) | |
| rc / FATAL / NON-CONVERGED | 0 / 0 / 0 | | | 0 / 0 / 0 | |
| s/cycle (cpu time / 10) | **1.71** (Raven 1.55, viper 2.57) | | | **1.57** (Raven 1.37, viper 2.48) | |

- Wall time for A+B together: 89 s (A 40 s, B 40 s, including IC setup). Cost: about 0.1 SU.
- A passes. B is about 1e-11 from viper, outside the ~1e-12 band; waiting on a Raven B number to tell whether
  that is the normal cross-GPU spread for B.
- Delta A100 is 10-15 % slower per cycle than Raven A100 and about 1.5x faster than a viper node.

## 5. Credentials
- `origin` in `/u/jma20/ATHENAK/athenak` is now `git@github.com:jing-ze-ma/athenak.git`, using an SSH key the
  user added on 10-08 (`~/.ssh/id_ed25519`, comment "delta"). No token in the remote URL.
- The machine-wide token scan from NOTE-2026-09-28-token-cleanup was **not run**: Claude Code's permission
  classifier blocked it. Delta is a new account, so a stray token is unlikely, but the user may run the scan
  by hand.
