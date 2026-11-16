# TASK 2026-11-16 Orion -> viper (USER): clean scheme GPU stage 1

First GPU run of the clean cell-centred implicit grey M1 scheme on **MI300A**: (1) build with
Kokkos HIP for gfx942 APU, (2) a short correctness check against Orion CPU references (G1-G4), and
(3) **G5 timing, the main point of this task**. We want clean vs the production scheme on
production-size blocks (548 x 64 x 64) on 1, 2, 4 and 8 MI300A, with a kernel breakdown and the
memory high-water mark. Every real 3-D production run is on GPUs, so G5 decides whether clean can
replace the production scheme.

Self-contained package: `docs/handover/task_1116_cc_gpu1/` (decks, runs.txt, run_all.sh,
compare.py with the Orion CPU reference in ref/, timing.py, build_hip.sh). Copy it somewhere
writable and run from there.
This TASK replaces the plan in NOTE-2026-11-16-orion-viper-cc_gpu_heads_up.md, which said CUDA
and 18309ae7a.

## Code
Fork `jing-ze-ma/athenak`, **origin/clean-cc = 2a8002df6** (push 8). This is the restyled clean code,
with the line solver now a hybrid PCR-Thomas team solve, fused BiCGStab reductions, ghost copies
of the shell only, production boundary conditions, esrc, a preflight check of the production
keys, and clean_layout = cell as the default.

Compile check on Orion: Orion has no ROCm, so we compiled 2a8002df6 with Kokkos CUDA instead
(nvcc 12.1, Ampere80, `--fmad=false`). That catches the usual device-lambda and host-call
errors. Result: 0 errors. The only warning in the clean files is nvcc 20013-D at rad_m1_clean_solver.cpp:224-225: CleanDots indexes a captured std::array inside a device lambda. That is a host-only constexpr operator under nvcc. It does not affect HIP, because clang treats constexpr functions as host+device. Orion has a one-line fix (local branch clean-cc-cuda 430488de4) for later CUDA builds; it is not needed for this task. The clean sources call no CUDA or HIP API directly. They use no NCCL
or RCCL and no atomics. The new line solver uses TeamPolicy with team scratch, so check the
team-size and scratch limits on HIP (team size <= 1024, LDS 64 KB per CU). The HIP-specific code
in rad_m1 is production code that your MI300A builds already compile.

## Build (MI300A, built-in pgens, no PROBLEM)
Run `./build_hip.sh <worktree at 2a8002df6>`. It is your MI300A he_star_m1 recipe plus `-ffp-contract=off`:
```
module purge; module load gcc/16 rocm/7.2 openmpi_gpu/5.0 cmake/4.0
cmake -S <wt> -B <wt>/build_hip_n432 -D Kokkos_ENABLE_HIP=On -D Kokkos_ARCH_AMD_GFX942_APU=On \
      -D CMAKE_CXX_COMPILER=$(which hipcc) -D Kokkos_ENABLE_IMPL_HIP_MALLOC_ASYNC=OFF \
      -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_FLAGS=-ffp-contract=off
# second binary: the same plus -D Athena_M1_VGD_NMAX=4096 (build_hip_n4096)
make -C <build> -j16 VERBOSE=1
```
- `-ffp-contract=off`: clang applies it to host AND device code (the HIP device default is
  fast-honor-pragmas). The fmad-1031 CMake logic is CUDA-only and does nothing here. The script
  counts the compile lines that carry the flag.
- NMAX 4096 only sizes the per-thread `Real buf[M1_VGD_NMAX]` of VetGdWall, a kernel nothing calls
  on this commit, plus a run-time cap on the direction count. Run everything with n432. Rerun
  G1 with n4096 (it should give the same answers).
- GPU-aware MPI on ROCm: your usual openmpi_gpu/5.0 (UCX) runtime and the runtime environment of
  your MI300A he_star_m1 runs. The production NCCL option (`Athena_M1_VGD_NCCL_DIR`, CUDA-only)
  stays unset.

## G5 timing (main deliverable; report only, no gates)
Deck: the radiative-hydrostatic envelope `rhse`: grey scattering kappa_s 500, gravity, gas
feedback, implicit radiation, vet_gd rays with GL 10 x 12 directions (the production set),
well-balancing off. Every MeshBlock is a production block, 548 x 64 x 64 (2.2 M cells).
Arms: **clean** (the commit's clean defaults), **prod** (the production scheme on the same deck:
scheme and clean keys removed), and **prodrot** (prod with production rotation: vet_rotate_avg 4,
vet_rotate_advance, vet_gd_rotate_xfade, vet_gd_rotate_dt 6.5e-3, ~5 cycles). nlim 25; timing.py
reports s/cycle over cycles 6-25 (cycles 0-5 are warm-up and initialisation).

| layout | mesh | blocks | MI300A | blocks/GPU | arms |
|---|---|---|---|---|---|
| sp1x2 | sp 548 x 64 x 128 | 2 | 1 | 2 | clean, prod, prodrot |
| sp1x4 | sp 548 x 128 x 128 | 4 | 1 | 4 | clean, prod |
| sp4x2 | sp 548 x 128 x 256 | 8 | 4 | 2 | clean, prod |
| sp8x2 | sp 548 x 256 x 256 | 16 | 8 | 2 | clean, prod |
| cs2x3 | cs 548 x 64^2 per panel | 6 | 2 | 3 | clean, prod |
| cs4x6 | cs 548 x 128^2 per panel | 24 | 4 | 6 | clean, prod |
| cs8x3 | cs 548 x 128^2 per panel (RSG E layout) | 24 | 8 | 3 | clean, prod, prodrot |

The sp wedges keep dtheta = dphi = 0.4/64 (theta reflecting about the equator, phi periodic).
```
export KOKKOS_TOOLS_LIBS=<kokkos-tools>/libkp_kernel_timer.so KP_READER=<kokkos-tools>/kp_reader
LAUNCH='srun -n %N --gpus-per-task=1 <your MI300A binding>' KEEPGOING=1 TIMEOUT=30m \
  ./run_all.sh <build_hip_n432>/src/athena G5
python3 timing.py > timing.txt
# memory: the clean and prod arms of sp1x2 and cs8x3 again with nlim 6 and the space-time-stack tool
KOKKOS_TOOLS_LIBS=<kokkos-tools>/libkp_space_time_stack.so NLIM=6 RUNS=$PWD/runs_mem \
  LAUNCH=... ./run_all.sh <build_hip_n432>/src/athena G5    # (stop after those four, or run all)
```
Start with sp1x2 (all three arms) and send its numbers as soon as you have them; then do the
rest. If a run is not at cycle 25 after 30 min, TIMEOUT stops it. timing.py still reports the
cycles it reached. Please send:
1. timing.txt: s/cycle per arm and layout, the clean/prod ratio, the initialisation time;
2. the kernel table per clean run (timing.py prints the top 15 from kp_reader, summed over
   ranks). Kernels of interest: cln_fv (FV ray sweep; timing.py also gives its time and passes
   per ray build), the line solver (clc_thomas / the hybrid PCR-Thomas kernel), clc_pcsetup,
   clc_apply, cln_dot / fused reductions, cln_hcopy / cln_hback / cln_ghost* (ghost copies);
   and the same table for a prod run;
3. MPI share, from your usual profiler if it is cheap to get (rocprofv3 --kernel-trace --stats
   is also fine instead of kokkos-tools; name what you used);
4. the memory high-water mark per GPU (space-time-stack "high water" for the HIP space and the
   host; MI300A memory is unified, so say which numbers you report).

## G1-G4 correctness (compact; 1-8 min per run on the CPU)
```
LAUNCH=... ./run_all.sh <build_hip_n432>/src/athena G1 G2 G3 G4
python3 compare.py > gates.md                       # tables + GATE FAILURES or "all gates passed"
LAUNCH=... RUNS=$PWD/runs_n4096 ./run_all.sh <build_hip_n4096>/src/athena G1
python3 compare.py --runs runs_n4096 G1 > gates_n4096.md
```
- **G1, 1 GPU, 1 block:** Cartesian IDORT decks nltek1 x0.5, hcyl x0.5, cb x0.5, shd x0.25, and the
  radiation-wave pair rwp/rwm (+-amp 1e-3, P 3.33, tau 10, N32, scored on the odd part).
- **G2, 4 GPUs / 4 blocks (x1 split):** nltek1, hcyl, shd, rw. cb is left out: vet_sc refuses cb x0.5
  on any multi-block layout ("per-ray band wider than the band"). This limitation already
  exists on 18309ae7a (push 7), so it is not caused by the second-order SC interpolation.
- **G3, spherical polar:** rhse_sp 2-D and 3-D wedge with well-balancing (force_reference = clean,
  wellbalance_static(+_reconst)), the 3-D wedge on 2x2 lateral blocks (4 GPUs), the thin pulse
  wedge pulsew, and the thin 3-D FV-ray wedge spthin.
- **G4, cubed sphere:** the rhse_cs WB envelope on 6 blocks (1 GPU) and 24 blocks (4 GPUs), cs3 thin1,
  cs3 feed1.
- Gates (compare.py):
  - rc 0, no FATAL or NaN, 0 CLEAN WARNING, 0 DT_RETRY.
  - The same number of CLEAN lines, with each line's `sum` and `pic` within +-1 of the CPU.
  - G1 score within max(1e-6, 2 x CPU thread spread).
  - G2 4-vs-1 field maxrel <= max(1e-7, 10 x the CPU's own 4-vs-1 value).
  - WB max v/c_s <= max(1e-7, 2 x CPU).
  - dE, dP, dL_z <= max(1e-13, 10 x |CPU|).

  Each widened gate uses a CPU number stored in ref/ref.json. Only nltek moves noticeably with
  the reduction order on the CPU, so only its gates are wider. All decks have every clean key
  written out; the explicit decks reproduce the defaults bitwise on Orion.

## Orion CPU reference (ref/ref.json, final dumps in ref/dumps/)
Commit 2a8002df6, CPU build (gcc/13, openmpi/4.1, MPI + OpenMP, Release, `-ffp-contract=off`), md5
e5cc6fa226c993572666ccb79f10dce8, 4 OpenMP threads, ranks as in runs.txt.

| run | score | CLEAN lines / sum / pic | CPU conservation, v/c_s | CPU 4-vs-1 maxrel | CPU wall s |
|---|---|---|---|---|---|
| g1_nltek1_x0.5 | L1 9.857942e-05 | 59 / 451 / 151 | - | - | 20 |
| g1_hcyl_x0.5 | metric 1.200974e-02 | 64 / 563 / 169 | - | - | 68 |
| g1_cb_x0.5 | metric 4.599888e-01 | 14 / 614 / 31 | - | - | 20 |
| g1_shd_x0.25 | metric 1.180539e-02 | 4 / 134 / 8 | - | - | 13 |
| g1_rwp | l1_rho 3.648449e-03, l1_E 3.582953e-03 | 127 / 432 / 684 | dE -1.2e-15, dP 2.4e-15 | - | 5 |
| g1_rwm | - | 127 / 432 / 684 | - | - | 6 |
| g2_nltek1_x0.5 | L1 9.889859e-05 | 59 / 476 / 150 | - | 4.1e-06 | 18 |
| g2_hcyl_x0.5 | metric 1.200974e-02 | 64 / 610 / 169 | - | 2.7e-10 | 39 |
| g2_shd_x0.25 | metric 1.180539e-02 | 4 / 150 / 8 | - | 2.5e-10 | 10 |
| g2_rwp | l1_rho 3.648450e-03, l1_E 3.582953e-03 | 127 / 432 / 685 | dE 3.1e-16, dP 1.1e-14 | 6.6e-08 | 12 |
| g2_rwm | - | 127 / 432 / 685 | - | 6.6e-08 | 11 |
| g3_rhse_sp2 | - | 1319 / 460 / 3317 | dE -1.3e-08, vcs_maxall 1.7e-07 | - | 17 |
| g3_rhse_sp3 | - | 341 / 200 / 876 | dE -1.1e-08, dLz 2.8e-15, vcs_maxall 1.1e-07 | - | 46 |
| g3_rhse_sp3_mb4 | - | 341 / 216 / 875 | dE -9.2e-09, dLz -9.5e-15, vcs_maxall 4.0e-07 | 2.7e-08 | 31 |
| g3_pulsew | - | 179 / 2091 / 1384 | - | - | 10 |
| g3_spthin | - | 39 / 82 / 87 | - | - | 11 |
| g4_rhse_cs6 | - | 341 / 159 / 869 | dE 1.4e-08, vcs_maxall 7.7e-08 | - | 206 |
| g4_rhse_cs24 | - | 341 / 160 / 869 | dE 1.4e-08, vcs_maxall 3.0e-07 | - | 102 |
| g4_cs_thin1 | - | 199 / 117 / 346 | - | - | 128 |
| g4_cs_feed1 | - | 19 / 250 / 202 | - | - | 40 |

CPU thread spread (1 and 16 vs 4 OpenMP threads): g1_nltek1_x0.5 score 2.4e-04 field 7.6e-07 dsum 2; g1_hcyl_x0.5 score 0.0e+00 field 7.2e-16 dsum 0; g1_cb_x0.5 score 1.6e-15 field 6.7e-16 dsum 0; g1_shd_x0.25 score 0.0e+00 field 4.1e-16 dsum 0; g1_rwp score 2.5e-11 field 9.7e-11 dsum 0; g1_rwm score 0.0e+00 field 9.0e-11 dsum 0

All CPU runs: rc 0, 0 CLEAN WARNING, 0 DT_RETRY. Running compare.py on Orion against the 16-thread runs gives "all gates passed".

## Rules
- Fresh starts only. Stop on FATAL, NaN or rc != 0 (G1-G4) and send the last 40 log lines;
  do not adapt decks or keys. grep the logs for `CLEAN WARNING` (compare.py counts them).
- G5: KEEPGOING=1, so that one slow or failed arm does not block the others.

## What to send back
1. Build: ROCm/hipcc/module versions, the Kokkos configure summary, compile errors and warnings in
   the clean files (`build_hip_n*.clean_diag.txt`), the -ffp-contract=off line count, md5 of both
   binaries.
2. G5: timing.txt, the kernel tables, the memory high-water marks, the MPI share (first sp1x2,
   then the rest).
3. gates.md and gates_n4096.md.
4. The run-dir paths on viper, `md5sum runs/*/bin/*.bin runs/*/tab/*.tab | sort`, and the run.log of
   any failed run.

Reply: NOTE-2026-11-16-viper-orion-cc_gpu1.md on fork/bsg-files-1009 (an early NOTE with the sp1x2
numbers is welcome).
