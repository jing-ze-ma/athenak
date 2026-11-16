# TASK 2026-11-16 Orion -> viper (USER): clean scheme GPU stage 1 on Raven A100 (amended)

**Replaces TASK-2026-11-16-viper-cc_gpu1 (MI300A plan): Raven A100.** The user moved the GPU tests
to Raven A100 (NVIDIA, CUDA, Ampere80). Please ignore the MI300A TASK. The package
`docs/handover/task_1116_cc_gpu1/` has been updated in place: CUDA build script, G5 layouts sized
for A100-40GB, and the commit. G1-G4 are unchanged.

What is asked: (1) the CUDA build; (2) a short correctness check against Orion CPU references
(G1-G4); (3) **G5 timing, the main point of this task**. G5 compares clean with the production
scheme on production-size blocks (548 x 64 x 64) on 1, 2, 4, 6 and 8 A100, with a kernel breakdown
and the memory high-water mark. Every real 3-D production run is on GPUs, so G5 decides whether
clean can replace the production scheme.

## Code
Fork `jing-ze-ma/athenak`, **origin/clean-cc = 32af6349e** (push 10). It is 4680ce500 plus the
GPU-parallel wavefront FV ray sweep: one launch per dependency level over all cells x directions,
instead of one serial thread per block and direction. The regression is bitwise identical to
4680ce500 on every output file, so the shipped references (made on 4680ce500) are valid for
32af6349e. 4680ce500 contains push 8 (2a8002df6), the CUDA fix of dae5c4133 (CleanDots copies its std::array
slot indices into ints before the device lambda), FV-ray rotation, the temporary keys
clean_dissip_power and clean_line_solver (read only when named), and clean_rays = lld. The defaults
the decks use are unchanged: hybrid line solver, A^2, sc rays on Cartesian, fv1 on sp/cs. 

Orion compiled 2a8002df6 with Kokkos CUDA (nvcc 12.1, Ampere80, `--fmad=false`): 0 errors. The
only warning in the clean files was the CleanDots one, which the fix removes. The clean code uses
no NCCL and no atomics. The line solver is a TeamPolicy kernel ("clc_line") with team scratch of
104 x nl Reals per line (~456 KB at nl = 548), so it runs from scratch level 1 (global memory).

## Build (Raven A100, built-in pgens, no PROBLEM)
Run `./build_cuda_raven.sh <worktree at origin/clean-cc 32af6349e>` (NMAX 432, the production GPU default). It uses
your Raven A100 production recipe:
```
module purge; module load gcc/13 cuda/12.6 openmpi_gpu/5.0 cmake
cmake -S <wt> -B <wt>/build_a100_n432 -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_AMPERE80=On \
      -D Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC=OFF -D CMAKE_CXX_COMPILER=<wt>/kokkos/bin/nvcc_wrapper \
      -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_FLAGS=-ffp-contract=off
make -C <wt>/build_a100_n432 -j16
```
- `-ffp-contract=off` makes the fmad-1031 CMake logic append `--fmad=false` for the device code.
  The script checks both the configure line and flags.make.
- GPU-aware MPI: openmpi_gpu/5.0 (CUDA-aware) with your usual Raven runtime settings.
- NMAX 4096 is optional (`./build_cuda_raven.sh <wt> 4096`). On this commit it only sizes
  `Real buf[M1_VGD_NMAX]` in VetGdWall (32 KB per thread), a kernel nothing calls, plus a run-time
  cap on the direction count. No run here needs more than 432 directions. If you build it, rerun G1.

## G5 timing (main deliverable; report only, no gates)
Deck: the radiative-hydrostatic envelope `rhse`: grey scattering kappa_s 500, gravity, gas
feedback, implicit radiation, vet_gd rays with GL 10 x 12 directions (the production set),
well-balancing off. Arms: **clean** (the commit's clean defaults), **prod** (the production scheme
on the same deck: scheme and clean keys removed), and **prodrot** (prod with production rotation:
vet_rotate_avg 4, vet_rotate_advance, vet_gd_rotate_xfade, vet_gd_rotate_dt 6.5e-3, ~5 cycles).
nlim 25; timing.py reports s/cycle over cycles 6-25.

**Memory.** Orion measured the CPU RSS slope with ghosts included: **clean ~11.5 kB per cell =
~30 GB per 548 x 64 x 64 block; prod ~4.6 kB per cell = ~12 GB**. So clean fits one production block
per A100-40GB. To see several blocks per GPU we use half and quarter blocks. If sp1x1_clean runs
out of memory, report it and use sp1x1h. The measured device high-water mark is one of the main
results.

| layout | mesh | block | blocks | A100 | blocks/GPU | clean GB/GPU (est.) | arms |
|---|---|---|---|---|---|---|---|
| sp1x1 | sp 548 x 64 x 64 | 548x64x64 | 1 | 1 | 1 | ~30 | clean, prod, prodrot |
| sp1x1h | sp 548 x 32 x 64 | 548x32x64 | 1 | 1 | 1 | ~15 | clean, prod |
| sp2x1 | sp 548 x 64 x 128 | 548x64x64 | 2 | 2 | 1 | ~30 | clean, prod |
| sp4x1 | sp 548 x 128 x 128 | 548x64x64 | 4 | 4 | 1 | ~30 | clean, prod |
| sp8x1 | sp 548 x 128 x 256 | 548x64x64 | 8 | 8 | 1 | ~30 | clean, prod |
| sp4x2h | sp 548 x 128 x 128 | 548x32x64 | 8 | 4 | 2 | ~30 | clean, prod |
| cs6x1 | cs 548 x 64^2 per panel | 548x64x64 | 6 | 6 | 1 | ~30 | clean, prod |
| cs8x3q | cs 548 x 64^2 per panel | 548x32x32 | 24 | 8 | 3 | ~23 | clean, prod, prodrot |

The sp wedges have dtheta = dphi = (pi/6)/64 (theta reflecting about the equator, phi periodic
with a period that is a multiple of 2 pi/12, as GL nphi 12 requires). cs8x3q is the RSG E
layout (24 blocks, 3 per GPU) with quarter blocks.
```
export KOKKOS_TOOLS_LIBS=<kokkos-tools>/libkp_kernel_timer.so KP_READER=<kokkos-tools>/kp_reader
LAUNCH='srun -n %N --gpus-per-task=1 <your Raven A100 binding>' KEEPGOING=1 TIMEOUT=30m \
  ./run_all.sh <build_a100_n432>/src/athena G5
python3 timing.py > timing.txt
# memory: clean and prod of sp1x1 and cs8x3q again with nlim 6 and the space-time-stack tool
KOKKOS_TOOLS_LIBS=<kokkos-tools>/libkp_space_time_stack.so NLIM=6 RUNS=$PWD/runs_mem \
  LAUNCH=... ./run_all.sh <build_a100_n432>/src/athena G5    # (stop after those four, or run all)
```
nsys instead of kokkos-tools is fine (`nsys profile -t cuda,nvtx,mpi` with kokkos-tools'
libkp_nvtx_connector.so for the kernel names). nvidia-smi `--query-gpu=memory.used` polling
also works for the high-water mark. Name what you used. Start with sp1x1 (all three arms) and
send its numbers as soon as you have them; then do the rest. TIMEOUT stops any run that is not
at cycle 25 after 30 min, and timing.py still reports the cycles it reached. Please send:
1. timing.txt: s/cycle per arm and layout, the clean/prod ratio, the initialisation time;
2. the kernel table per clean run (top 15, summed over ranks) and for a prod run. Kernels of
   interest: cln_fv (wavefront FV ray sweep; timing.py also gives its time and launches per ray build),
   clc_line (hybrid PCR-Thomas line solver), clc_pcsetup, clc_apply, cln_dot (fused
   reductions), cln_hcopy / cln_hback / cln_ghost* (ghost copies);
3. the MPI share, if cheap (nsys gives it);
4. the device memory high-water mark per GPU for clean and prod.

## G1-G4 correctness (compact; 1-8 min per run on the CPU)
```
LAUNCH=... ./run_all.sh <build_a100_n432>/src/athena G1 G2 G3 G4
python3 compare.py > gates.md                       # tables + GATE FAILURES or "all gates passed"
# optional, only if you build n4096:
LAUNCH=... RUNS=$PWD/runs_n4096 ./run_all.sh <build_a100_n4096>/src/athena G1
python3 compare.py --runs runs_n4096 G1 > gates_n4096.md
```
- **G1, 1 GPU, 1 block:** Cartesian IDORT decks nltek1 x0.5, hcyl x0.5, cb x0.5, shd x0.25, and the
  radiation-wave pair rwp/rwm (+-amp 1e-3, P 3.33, tau 10, N32, scored on the odd part).
- **G2, 4 GPUs / 4 blocks (x1 split):** nltek1, hcyl, cb, shd, rw. (cb x0.5 on several blocks runs
  on 4680ce500; up to 2a8002df6 vet_sc refused it with "per-ray band wider than the band".)
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
Commit 4680ce500 (outputs bitwise identical to 32af6349e), CPU build (gcc/13, openmpi/4.1, MPI + OpenMP, Release, `-ffp-contract=off`), md5
1326bb7410d558b5cec09a38b6f87cd5, 4 OpenMP threads, ranks as in runs.txt. All decks have every clean key written out
(taken from the parameter dump); on Orion the explicit decks reproduce the implicit defaults bitwise.

| run | score | CLEAN lines / sum / pic | CPU conservation, v/c_s | CPU 4-vs-1 maxrel | CPU wall s |
|---|---|---|---|---|---|
| g1_nltek1_x0.5 | L1 9.857942e-05 | 59 / 451 / 151 | - | - | 22 |
| g1_hcyl_x0.5 | metric 1.200974e-02 | 64 / 563 / 169 | - | - | 80 |
| g1_cb_x0.5 | metric 4.599888e-01 | 14 / 614 / 31 | - | - | 23 |
| g1_shd_x0.25 | metric 1.180539e-02 | 4 / 134 / 8 | - | - | 15 |
| g1_rwp | l1_rho 3.648449e-03, l1_E 3.582953e-03 | 127 / 432 / 684 | dE -1.2e-15, dP 2.4e-15 | - | 8 |
| g1_rwm | - | 127 / 432 / 684 | - | - | 6 |
| g2_nltek1_x0.5 | L1 9.889859e-05 | 59 / 476 / 150 | - | 4.1e-06 | 23 |
| g2_hcyl_x0.5 | metric 1.200974e-02 | 64 / 610 / 169 | - | 2.7e-10 | 46 |
| g2_cb_x0.5 | metric 4.599888e-01 | 14 / 697 / 31 | - | 3.6e-07 | 16 |
| g2_shd_x0.25 | metric 1.180539e-02 | 4 / 150 / 8 | - | 2.5e-10 | 17 |
| g2_rwp | l1_rho 3.648450e-03, l1_E 3.582953e-03 | 127 / 432 / 685 | dE 3.1e-16, dP 1.1e-14 | 6.6e-08 | 18 |
| g2_rwm | - | 127 / 432 / 685 | - | 6.6e-08 | 19 |
| g3_rhse_sp2 | - | 1322 / 2700 / 5286 | dE 1.2e-05, vcs_maxall 7.4e-06 | - | 28 |
| g3_rhse_sp3 | - | 347 / 962 / 1446 | dE 1.5e-03, dLz -1.6e-08, vcs_maxall 1.1e-03 | - | 81 |
| g3_rhse_sp3_mb4 | - | 347 / 1072 / 1441 | dE 1.5e-03, dLz -1.6e-08, vcs_maxall 4.3e-03 | 1.1e-07 | 54 |
| g3_pulsew | - | 179 / 2091 / 1384 | - | - | 10 |
| g3_spthin | - | 39 / 82 / 87 | - | - | 11 |
| g4_rhse_cs6 | - | 377 / 2298 / 1811 | dE 2.3e-04, vcs_maxall 1.2e-03 | - | 798 |
| g4_rhse_cs24 | - | 377 / 2269 / 1816 | dE 2.6e-04, vcs_maxall 3.8e-03 | - | 251 |
| g4_cs_thin1 | - | 199 / 117 / 346 | - | - | 145 |
| g4_cs_feed1 | - | 19 / 250 / 202 | - | - | 42 |

CPU thread spread (1 and 16 vs 4 OpenMP threads): g1_nltek1_x0.5 score 2.4e-04 field 7.6e-07 dsum 2; g1_hcyl_x0.5 score 0.0e+00 field 7.2e-16 dsum 0; g1_cb_x0.5 score 1.6e-15 field 6.7e-16 dsum 0; g1_shd_x0.25 score 0.0e+00 field 4.1e-16 dsum 0; g1_rwp score 2.5e-11 field 9.7e-11 dsum 0; g1_rwm score 0.0e+00 field 9.0e-11 dsum 0

All CPU runs: rc 0, 0 CLEAN WARNING, 0 DT_RETRY. Running compare.py on Orion against the 16-thread runs gives "all gates passed".

## Rules
- Fresh starts only. Stop on FATAL, NaN or rc != 0 (G1-G4) and send the last 40 log lines;
  do not adapt decks or keys. grep the logs for `CLEAN WARNING` (compare.py counts them).
- G5: KEEPGOING=1, so that one slow or failed arm does not block the others.

## What to send back
1. Build: modules (gcc/cuda/openmpi_gpu), the Kokkos configure summary, the `fmad-1031` configure
   line and `--fmad=false` in flags.make, compile errors and warnings in the clean files
   (`build_a100_n*.clean_diag.txt`), binary md5(s).
2. G5: timing.txt, the kernel tables, the memory high-water marks, the MPI share (first sp1x2,
   then the rest).
3. gates.md (and gates_n4096.md if built).
4. The run-dir paths on viper, `md5sum runs/*/bin/*.bin runs/*/tab/*.tab | sort`, and the run.log of
   any failed run.

Reply: NOTE-2026-11-16-viper-orion-cc_gpu1.md on fork/bsg-files-1009 (an early NOTE with the sp1x1
numbers is welcome).
