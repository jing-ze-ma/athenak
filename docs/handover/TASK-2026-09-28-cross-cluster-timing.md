# TASK 2026-09-28: cross-cluster timing benchmark, WASP-121b 1x and 10x

The same benchmark runs on three machines: viper (MI300A, HIP), Caltech (H200, CUDA) and NCSA DeltaAI
(GH200, CUDA with an ARM Grace host). The goal is one table of ms/cycle and wall per simulated second.
The package is in `docs/handover/bench-2026-09-28/`:

| file | what it is |
|---|---|
| `inputs/w121_bench_1x.athinput`, `inputs/w121_bench_10x.athinput` | the production inputs `w121prod_{1x,10x}` with the changes listed below |
| `inputs/ic_w121_1x.txt`, `inputs/ic_w121.txt` | IC profiles (1x, 10x); `run_bench.sh` copies them into the run dir |
| `run_bench.sh` | one fresh-start run of 2000 cycles: `run_bench.sh <binary> <nranks> <1x\|10x> <new run dir>` |
| `ana_bench.py` | reads the logs and prints one line per run |
| `bench_viper.sub` | the viper job script (reference for the other machines) |
| `RESULTS_viper.md` | the viper rows |

## Pinned commit

**11c9a5be** on fork/rt-integration (`11c9a5be7d5cb8fc00a144f48f20897282df92d9`), with Kokkos at the in-tree
submodule d8e9af03 (4.6.02). Build it from `git archive 11c9a5be`, not from a working tree. The docs commit
that adds this package comes later than 11c9a5be. Only `docs/` differs, so the binary is the same.

## Inputs: what changed from production

The source is `/viper/ptmp2/jinma/w121prod_0927/w121prod_{1x,10x}.athinput` on viper. Every changed line is
marked `bench-2026-09-28`:

- **Data paths** are the placeholders `ATHENAK_CK_DATA` (1x tables, the repo `data/exo_fms_ck`) and
  `ATHENAK_CK_DATA10` (10x `ckdata10`). They follow the Caltech handover convention. `run_bench.sh`
  substitutes `CKDATA` and `CKDATA10`.
- **IC profile** is read from the run dir.
- **Newton fixes** come from `/viper/ptmp2/jinma/cknewton_0928/RESULTS.md`:
  - 1x: `ck_impl_maxit = 16`.
  - 10x: `ck_impl_maxit = 24` and `ck_impl_dtmax = 0.25`.
- `ck_impl_verbose = false`.
- `time/nlim = 2000` (fresh start) and `time/ndiag = 8`.
- **Outputs:** hst only. There are no bin, rst or log outputs during timing.

Grid: C32 cubed sphere, 6 x 2 x 2 = **24 MeshBlocks** (nx1 76 for 1x, 74 for 10x; 16 x 16 per block).
This means at most 24 ranks.

## Build per machine

All three builds use MPI, Release and `-D PROBLEM=deep_hot_jupiter_rt`.

| machine | recipe | Kokkos arch | compiler / modules |
|---|---|---|---|
| viper | `/viper/ptmp2/jinma/bench_0928/build.sh` (a copy of `w121prod_0927/build.sh` at 11c9a5be) | `Kokkos_ENABLE_HIP=On`, `Kokkos_ARCH_AMD_GFX942_APU=On` | gcc/14 rocm/6.3 (hipcc 6.3.4) openmpi_gpu/5.0 |
| Caltech | `docs/handover/caltech-2026-09-26/scripts/build_caltech.sh bench gpu 11c9a5be` | `Kokkos_ENABLE_CUDA=On`, `Kokkos_ARCH_HOPPER90=On`, nvcc_wrapper | gcc/13.2.0, cuda/12.9.0, hpcx/2.17.1 (CUDA-aware) |
| DeltaAI | per `TASK-2026-09-28-deltaai-bringup.md` section 2 | `Kokkos_ARCH_HOPPER90=On`, `Kokkos_ARCH_ARMV9_GRACE=On`, nvcc_wrapper | to be recorded by that session |

Record the md5 of the binary. `run_bench.sh` also writes it into every log.

## What to run

- **GPU counts:** 1, 2 and 4 GPUs, with 1 MPI rank per GPU.
  - 4 GPUs is 1 node on Caltech and DeltaAI.
  - On viper, 4 GPUs is 2 nodes x 2 (`-p apu -N 2`).
- **Arms:** 1x and 10x.
- **Repeats:** at least 2. Interleave the arms in one job (1x, 10x, 1x, 10x), with one job per GPU count,
  so that each job uses one node set and one binary.
- **Smoke first:** run 50 cycles with `EXTRA="time/nlim=50 time/ndiag=1 problem/ck_impl_verbose=true"`.
  Check that there is no FATAL and count the `NOT-CONVERGED` lines.
- **Analysis:** `python3 ana_bench.py <run>/bench.log ...`. It prints:
  - **ms/cyc:** the median over the 8-cycle windows in cycles 1000-2000. This window skips the start-up
    transient.
  - **wall/sim-s:** wall seconds per simulated second, over the same windows. Report it because dt differs
    slightly between machines.
  - Both are also given over the whole window.

## Machine notes

- **viper:**
  - Use `bench_viper.sub`. Use apudev for 1 and 2 GPUs (under 15 min) and apu for 4 GPUs.
  - Set `HSA_XNACK=1`, `HSA_NO_SCRATCH_RECLAIM=1` and `SBATCH_EXPORT=NONE` (pass variables with `--export`).
- **Caltech:**
  - Use `srun --mpi=pmix` (the default pmi2 starts singleton ranks).
  - **Time on H200 only** (`--gres=gpu:nvidia_h200:<n>`).
  - Exclude **hpc-sm-01-09 and hpc-sm-02-16**.
  - An H100 row is allowed only if it is labelled H100.
  - Use `MACHINE=caltech` in `run_bench.sh`. It loads the build modules and sets the launcher.
- **DeltaAI:** fill the `deltaai)` header in `run_bench.sh` (modules, launcher, data dirs) in your results
  branch. Do not edit it on rt-integration.

## Rules for the other sessions

- Run timing jobs only: no production runs and no code changes.
- Use exactly 11c9a5be.
- Report results as a results file `docs/handover/bench-2026-09-28/RESULTS_<machine>.md`. Push it to a
  **new branch `bench-results-<machine>`**, not to rt-integration, or relay it through the user.
- Report measured numbers only, with the job ids and log paths.

## Results table template

| machine | GPU | ranks | 1x ms/cycle | 1x wall/sim-s | 10x ms/cycle | 10x wall/sim-s | commit | compiler / modules |
|---|---|---|---|---|---|---|---|---|
| viper | MI300A | 1 | | | | | 11c9a5be | gcc/14 rocm/6.3 openmpi_gpu/5.0 |
| Caltech | H200 | 1 | | | | | 11c9a5be | gcc/13.2.0 cuda/12.9.0 hpcx/2.17.1 |
| DeltaAI | GH200 | 1 | | | | | 11c9a5be | |

Give the median of the repeats, with the spread of the repeats in brackets (min-max). List the job ids
under the table.
