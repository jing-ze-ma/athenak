# TASK 2026-09-29: cross-cluster timing benchmark, 3-D He box with implicit M1 + vet_sc

> **HOLD LIFTED (user 09-29): run the benchmark.** The earlier hold said the M1 He box does not convect; that was a
> measurement error (a plane-mean flux that cannot see this very inefficient FeCZ, and cm/s read as km/s). The box does
> convect: v1'/v_MLT 1.2-1.5 saturated, r(v1',T') +0.22..+0.47, same as the two-stream box_w8
> (/viper/ptmp2/jinma/heconv_0929/RESULTS.md). Resolution checked: no doubling needed. The package is unchanged.

For the Caltech (H200) and DeltaAI (GH200) sessions; viper (MI300A) has run its arm. Same structure
and rules as `TASK-2026-09-28-cross-cluster-timing.md` (WASP-121b). The goal is one table of
ms/cycle and wall per simulated second for the implicit-M1 He box. **Timing only: no production, no
code changes.**

Package: `docs/handover/bench-2026-09-29-hebox/` (README.md has the input differences, the window
and the smoke; RESULTS_viper.md has the viper rows).

| file | what it is |
|---|---|
| `inputs/hebox_bench.athinput` | box_convection, 84 x 104 x 104, 4 MeshBlocks, lhllc + plm, polytropic WB, etotgrav, implicit M1, closure vet_sc (full tensor), hesdirk2, **cfl 0.3**, ke-dt fix keys, hst only, nlim 1200 |
| `run_bench_hebox.sh` | one fresh-start run: `run_bench_hebox.sh <binary> <1\|2\|4> <new run dir>` |
| `ana_bench_hebox.py` | one line per run: ms/cycle, wall/sim-s, Picard and Krylov counts, NC, fallbacks |
| `bench_viper_hebox.sub` | the viper job script (reference) |

## Pinned commit

**401875f0** on fork/rt-integration (`401875f01beb8cd69c215eb096bbebd020e0886b`), Kokkos = in-tree
submodule. Build `-D PROBLEM=box_convection`, MPI on, Release. The docs commit adding this package
is later; only `docs/` differs.

## Data (fresh start, no restart)

`HE_BOX_DATA` = the `he_box/` directory of `athenak_data_2026-09-26.tar.gz` (you already have it
from the 09-26 handover): `ic_m1_V3edd_pgen.txt m1_rad_ic_V3edd.txt rosseland_he_x0.0_z0.02.txt
planck_he_x0.0_z0.02_ferg+tops.txt arad_V3edd.txt`. `run_bench_hebox.sh` substitutes the
placeholder and checks that the five files exist. Set `HE_BOX_DATA=<dir>` if the machine header
default is wrong.

## Build per machine (incremental, `NOTE-2026-09-28-incremental-builds.md`)

| machine | recipe | Kokkos arch | compiler / modules |
|---|---|---|---|
| viper | `/viper/ptmp2/jinma/builds/build_inc_viper.sh box_gpu72 401875f0` | `Kokkos_ENABLE_HIP=On`, `Kokkos_ARCH_AMD_GFX942_APU=On` | gcc/16 rocm/7.2 (hipcc 7.2.4) openmpi_gpu/5.0 |
| Caltech | `docs/handover/caltech-2026-09-26/scripts/build_inc.sh hebox gpu 401875f0 box_convection` (on a compute node, see its header) | `Kokkos_ENABLE_CUDA=On`, `Kokkos_ARCH_HOPPER90=On`, nvcc_wrapper | gcc/13.2.0, cuda/12.9.0, hpcx/2.17.1 |
| DeltaAI | your incremental script, new target `box_gpu` (PROBLEM=box_convection) | `Kokkos_ARCH_HOPPER90=On`, `Kokkos_ARCH_ARMV9_GRACE=On`, nvcc_wrapper | as recorded for bench-2026-09-28 |

Record the binary md5 (`run_bench_hebox.sh` also writes it into every log).

## GPU counts: 1, 2, 4 (1 MPI rank per GPU)

The box has **4 MeshBlocks** (84 x 52 x 52), so 4 ranks is the maximum and 1 / 2 / 4 are the legal
counts (the script refuses others). vet_sc on several blocks/ranks needs a uniform mesh, periodic
x2/x3 and a ray reach per layer no wider than a block: here 1-2 cells vs 52, fine; the 64 rays are
not split over ranks (`vet_mb_agroup` = 1 default). Set no `vet_mb_*` key. 4 GPUs is 1 node on
Caltech and DeltaAI (viper: 2 nodes x 2).

## What to run

1. **Smoke** (1 job, 2 GPUs): `EXTRA="time/nlim=200 time/ndiag=1 rad_m1/implicit_picard_log=20"`,
   then `python3 ana_bench_hebox.py <run>/bench.log --c0 50 --c1 200`. Pass = rc 0, FATAL 0, nan 0,
   NC 0, fb 0. Viper smoke: 102 ms/cycle at 2 GPUs, Picard 6.76 / solve.
2. **Timing:** one job per GPU count (1, 2, 4), **2 repeats** back to back in the job
   (`for r in 1 2; do run_bench_hebox.sh $BIN $n $OUT/n${n}_r$r; done`), same binary everywhere.
   Each run is 1200 cycles; viper wall per run incl. start-up: see RESULTS_viper.md.
3. **Analysis:** `python3 ana_bench_hebox.py $OUT/*/bench.log` (window cycles 400-1200).
   Report ms/cyc (median of the 10-cycle windows), wall/sim-s, dt, pic/solve, kry/lin, NC, fb.
   The dt (~0.161) should agree across machines to a few 1e-4; the counts to a few %; if they
   differ more, say so.

## Machine notes

- **Caltech:** `MACHINE=caltech`; `srun --mpi=pmix`; **time on H200 only**
  (`--gres=gpu:nvidia_h200:<n>`), exclude hpc-sm-01-09 and hpc-sm-02-16; an H100 row only if labelled.
- **DeltaAI:** fill the `deltaai)` header of `run_bench_hebox.sh` (modules, launcher, HE_BOX_DATA)
  in your results branch only (copy it from your bench-2026-09-28 header).
- **viper:** `bench_viper_hebox.sub`, apudev for 1 and 2 GPUs, apu 2 nodes for 4;
  HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1, SBATCH_EXPORT=NONE.

## Rules

- Timing jobs only, exactly 401875f0, no code changes, no production.
- Results file `docs/handover/bench-2026-09-29-hebox/RESULTS_<machine>.md`, pushed to a **new branch
  `bench-results-<machine>`** (not rt-integration), or relayed through the user.
- Measured numbers only, with job ids and log paths. The rsolver note of 09-29 (lhllc for dhj)
  does not apply: the box already uses lhllc.

## Results table template

| machine | GPU | GPUs (ranks) | ms/cycle | wall s / sim s | dt | Picard/solve | Krylov it/linear solve | NC / fb | commit | compiler / modules |
|---|---|---|---|---|---|---|---|---|---|---|
| viper | MI300A | 1 (1) | | | | | | | 401875f0 | gcc/16 rocm/7.2 openmpi_gpu/5.0 |
| Caltech | H200 | 1 (1) | | | | | | | 401875f0 | gcc/13.2.0 cuda/12.9.0 hpcx/2.17.1 |
| DeltaAI | GH200 | 1 (1) | | | | | | | 401875f0 | |

Give both repeats (r1 / r2) or the median with the spread; list the job ids under the table.
