# TASK for the Delta session: bring up AthenaK on NCSA Delta (2026-10-08)

You are a Claude Code session on **NCSA Delta** (the user's new account; NOT DeltaAI, which is a separate GH200
machine the user already uses). The user develops AthenaK on viper (MPCDF, AMD MI300A, HIP), Raven (MPCDF, 4 x A100
per node, CUDA), DeltaAI (GH200), Caltech (H100/H200) and Orion (CPU). Goal here: build AthenaK, prove it runs
correctly against a viper/Raven reference, measure speed and queue behaviour. **No production runs** until the user
decides. Read this file to the end before doing anything.

## 0. Code and data
```bash
git clone --depth 1 --branch rt-integration --single-branch --recurse-submodules --shallow-submodules \
    https://github.com/jing-ze-ma/athenak.git athenak
```
- Expect rt-integration at 03590601 or newer. Kokkos is the in-tree submodule; do not use a site Kokkos.
- Read `CLAUDE.md` (build, tests, style, architecture), `docs/handover/NOTE-2026-09-28-token-cleanup.md` (do the
  one-time credential cleanup) and `docs/handover/NOTE-2026-09-28-incremental-builds.md` (build incrementally).
- Smoke files (AG Car, 1.7 MB) live on branch `he-ic-eint-from-t`, dir `docs/handover/agcar-files-1008/`
  (ICs, TOPS tables, inputs with `@AGCAR_DIR@`, `SETUP.sh`, `MD5SUMS`, viper references in `smoke_ref_viper/`).
  Get them without switching branch, e.g. `git fetch --depth 1 origin he-ic-eint-from-t` then
  `git archive FETCH_HEAD docs/handover/agcar-files-1008 | tar -x -C <work dir>`; run `SETUP.sh <dir>` there.
- Pushing needs GitHub credentials. Until the user sets up an SSH key on Delta, do not push; write results to files
  and tell the user, who relays them.

## 1. Find out the machine first (verify, do not assume)
From general knowledge (CHECK with `sinfo -s`, `module avail`, `sacctmgr show assoc user=$USER
format=account,partition,qos`, `nvidia-smi` / `rocm-smi` on a compute node, `accounts` or the NCSA Delta docs):
- GPU partitions like `gpuA100x4` (4 x A100 40 GB, AMD Milan host), `gpuA100x8`, `gpuA40x4`, maybe `gpuH200x8`,
  `gpuMI100x8`, plus `*-interactive` variants; CPU partition `cpu`. Accounts look like `<proj>-delta-gpu` /
  `-delta-cpu`. Filesystems: home (small), `/projects/<proj>`, `/work/hdd|nvme/<proj>` (check which exist, quotas,
  purge policy). Record SU charging per GPU-hour and job limits (max time, max nodes).

## 2. Build
- CUDA build for A100: `-D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_AMPERE80=On -D Kokkos_ARCH_ZEN3=On
  -D CMAKE_CXX_COMPILER=$PWD/kokkos/bin/nvcc_wrapper -D Athena_ENABLE_MPI=ON` (A40: AMPERE86). Pick gcc + CUDA 12.x +
  a CUDA-aware MPI from the modules; record the exact module list in a `build_delta.sh`.
- Targets: (a) `PROBLEM=he_star_m1` (the AG Car / He giant / BSG problem generator), (b) a plain build without
  PROBLEM for the test suite. A100 needs the scratch fix de528da5 (already in rt-integration).
- Compile big CUDA builds in a compute/interactive job, not on the login node, if the login node is weak.
- CPU check: `cd tst && python run_test_suite.py --test test_suite/hydro/test_hydro_fofc_cs_cpu.py` and the
  `rad_m1` slab test (needs ATHENAK_M1_DATA; skip if absent and say so).

## 3. GPU smoke vs reference (1 node, 4 A100, 4 ranks = one 480x64x64 MeshBlock per GPU)
Run `time/nlim=10` for AG Car B (`agcar_shakeB_ge.athinput`) and A (`agcar_shakeA_ge.athinput`), env as Raven:
one rank per GPU, `KOKKOS_MAP_DEVICE_ID_BY=mpi_rank` (or Kokkos' default mapping; check), `srun --cpu-bind=cores`.
Compare the hst after 10 cycles with `smoke_ref_viper/smoke{A,B}.hydro.hst` (viper 12131991 A / 12130898 B) and
the Raven A100 numbers:

| | A (viper 12131991) | A (Raven A100 31006215) | B (viper 12130898) |
|---|---|---|---|
| t | 1.2525591814773279e+04 | 1.2525591814773266e+04 | 2.4186988199012844e+03 |
| mass (col 3) | 6.1192879615027439e+32 | 6.1192879615024124e+32 | 1.7891809043883109e+31 |
| tot-E (col 7) | 1.6552578844679151e+47 | – | 1.0633773260127629e+46 |
| Picard mean / max | 13.5 / 16 | 13.6 / 16 | 6.9 / 24 |

Pass: rc 0, 0 FATAL, 0 NON-CONVERGED, t/mass/tot-E within ~1e-12 relative of Raven/viper (different GPU and
compiler: not bitwise), the two IC T-check lines ("IC column T(rho,eint)/T_col" ~6e-14, "he_ic_balance cells"
3.24e-7 A / 1.08e-7 B). Record s/cycle (Raven node: 1.37 B / 1.55 A s/cycle; viper node 2.5).

## 4. Report (stop after this)
Write `NOTE-2026-10-08-delta-bringup.md` (results file; push only once credentials exist): machine facts
(partitions, account, limits, SU rate, filesystems), `build_delta.sh`, binary md5, test results, smoke table vs the
references, s/cycle, how long the smoke jobs waited in the queue. Do not start anything longer than the smokes.
