# runs_4k_launch: fewer launches and host syncs in the implicit-M1 step

- **Date:** 2026-09-24, viper.
- **Branch:** `m1-launch`, based on m1-accmerge ac1e2892 (rt-integration + m1-accel +
  levers default on for be). Worktree: `/viper/ptmp2/jinma/wt_launch`.
- **Run tree:** `/viper/ptmp2/jinma/launch_0923`:
  - `base/` = git archive 0841ea5d (m1-accel);
  - `base2/` = ac1e2892;
  - `new/` = ac1e2892 plus this work;
  - `bin/`, `inp/`, `cpu/`, `runs/<tag>/<arm>_<rep>/`, `scripts/`. The scripts are copied
    to `scripts/` here.
- **GPU jobs:** apudev, `HSA_XNACK=1` and `HSA_NO_SCRATCH_RECLAIM=1`. Same binary, arms
  interleaved, even repeats in reverse order.
- **Metric:** ms/cycle over cycles 20-110 of 120 (`scripts/tsum.py <tag>`).
- **Box:** 3-D He box 84x104x104 in 4 blocks, 1 GPU.
- **Common settings:** L = implicit_vimp_fold + fast_kernels + one_pass=4 +
  predictor_order=2. V = closure vet_sc, vet_tensor full.

## 1. What the 1-GPU profile says (the brief's premise does not hold on 1 GPU)

Source: rocprofv3 kernel + HIP runtime trace, cycles 10-30, rank 0. Scripts:
`scripts/split.py <prof> 0 30 10`. Runs:
- base 0841ea5d: `runs/p0/*_1/prof` (job 11954699);
- ac1e2892 with the switch off: `runs/p1/V_off_1` (job 11955080).

| arm | ms/cyc (profiled) | GPU busy | GPU idle | launches/cyc | Stream+Device+Event syncs/cyc |
|---|---|---|---|---|---|
| V be+L (0841ea5d) | 54.1 | 46.9 | 7.1 (13 %) | 685 | 339 (226.6 + 48.7 + 63.5) |
| V be+L (ac1e2892, off) | 49.7 | 42.9 | 6.8 | 685 | about the same |
| E be+L (0841ea5d) | 47.6 | 39.6 | 8.0 | 518 | - |
| V hesdirk2+L (0841ea5d) | 72.1 | 62.4 | 9.8 | 919 | - |
| hydro only | 20.9 | 19.1 | 1.8 | 79 | 114.5 (43.0 of them hipDeviceSync) |

**On 1 GPU the step is GPU-time bound, not launch bound.**
- GPU idle is 13 %, not 30 %.
- The radiation part is idle for about 4.6 ms per cycle.
- The kernel time of the radiation part is about 24 ms per cycle (V be+L, 0841ea5d):
  - `ImplicitStencilOp`: 9.2 ms (107 us per call; 76 us with Eddington);
  - `M1PCRX`: 6.1 ms (36 us per half-sweep);
  - Picard assembly: 4.7 ms;
  - SC sweep: 1.4 + about 2 ms (VetRay kernels);
  - Opacity: 0.9 ms.
- **43 of the 48.7 hipDeviceSync per cycle come from hydro/conduction**, not from M1.
  The 1-GPU hydro-only run has 43.0 of them.
- The fences in `rad_m1_vet.cpp` (lines ~433, 456, 616, 636, 737, 1207, 1224, 1386,
  2312, 2615, 2722, 2983) are m1-sctb territory and were left alone.
  - Most of them sit before MPI sends.
  - 2722 and 2983 are timer fences around VetShortChar (2 per sweep).
  - `ImplicitSolve` has 2 timer fences per solve under vet_sc (lines ~5347 and 7227).

## 2. Step (a)+(b): `<rad_m1>/implicit_krylov_dev = K` (new file `src/rad_m1/rad_m1_launch.cpp`)

**The switch:**
- Default 0 = off. It is read only when named.
- It needs implicit_krylov_fuse = 3, implicit_halo_direct and no pipe.
- It acts only on 1 rank with halo_direct_on, the stencil operator, and no unfolded vimp.
  Anywhere else the solve runs ImplicitBiCGStabTwo unchanged.

**What it does:**
- alpha, beta, omega, rho, the flags and the iteration count live in the device array `kd`.
- The two reductions of each iteration reduce into a device View. The functor's `final`
  (run once on the device) advances the scalars and the convergence and breakdown flags.
  **No reduction returns to the host.**
- An iteration takes 6 launches:
  - preconditioner red, which also applies the previous iteration's x, r update and the
    p update;
  - preconditioner black;
  - A y, which reads the ghosts straight from the neighbour block, plus (rhat,v) and max|r|;
  - preconditioner red, which also forms s;
  - preconditioner black;
  - A z, plus (t,s), (t,t) and (rhat,t).
- Two's update kernel and both halo kernels are gone.
- **Batches:** the host queues iterations in batches. The first batch has the length of
  the last solve in the same Picard slot + 1; after that, K per batch. The host reads the
  status once per batch. After a stop, the queued kernels return at once.
- The true-residual confirmation, restarts and the fallback are Two's host code.
- `<rad_m1>/implicit_krylov_dev_halo = 1` (read only when named): the ghosts come from the
  direct-halo kernel and the operator reads them in place.

**Gates (CPU, login node, `scripts/cpu_s1.txt` and `cpu_s1c.txt`, `scripts/cmp.py`):**
- **Off vs ac1e2892: BITWISE (hst).**
  - 2-D He slab, t = 100: Eddington 1 and 2 ranks, vet_sc 1 rank.
  - 3-D box 84x32x32, 30 cycles, 4 ranks.
- **On vs off, CPU: BITWISE (hst)** in every case run:
  - slab Eddington, 1 rank: K = 4, K = 2, and K = 2 with dev_halo = 1;
  - slab, 2 ranks: falls back to Two, bitwise;
  - 3-D box, Eddington and vet_sc, K = 4;
  - 3-D box, vet_sc, K = 2 with dev_halo = 1;
  - 3-D box, hesdirk2 + vet_sc + L, K = 4.
  - The CPU reductions happen to sum in the same order.
- **NON-CONVERGED = 0 in all 22 CPU runs and in every GPU run.**
- **NOT DONE:**
  - the restart gate, switch on;
  - a GPU hst comparison on vs off (expected round-off: GPU reduction order).
  - GPU iteration counts are close: inner 38.28 vs 39.31 per step (V), 27.32 vs 27.34 (Vh).

**GPU cost, 1 GPU (ms/cycle, 2 repeats each):**

| arm | ac1e2892 (base2) | off | K=4 (s1, job 11954923) | K=2 (s1b, job 11955114) | K=1 (s1b) |
|---|---|---|---|---|---|
| Eddington be+L | 35.4, 36.2 | 35.3, 36.1 / s1b: 37.4, 36.2 | 35.4, 35.7 | **34.6, 34.6** | - |
| vet_sc be+L | 42.5, 43.1 | 43.2, 43.3 / s1b: 43.1, 43.4 | 44.9, 44.9 | **42.9, 42.9** | 43.0, 42.9 |
| vet_sc hesdirk2+L | - | 59.8, 59.8 / s1b: 59.6, 59.8 | 62.2, 62.3 | **59.4, 59.4** | - |

- s1 (K = 4) is before the x1 fast path and the per-slot batch prediction. It is slower
  because the operator looked up every x1 ghost too, and more no-op iterations were queued.
- **K = 2 vs off, same binary:** Eddington -2.2 ms (6 %), vet_sc -0.4, hesdirk2 -0.3.

**Launches and host syncs per cycle** (V be+L, K = 2 against off, profile `runs/p2/V_dev2_1`,
job 11955115; off = `runs/p1/V_off_1`):

| | off | K = 2 |
|---|---|---|
| launches | 685 | 580 |
| Stream + Device + Event syncs | 339 | **175.5** (63.3 + 48.7 + 63.5) |
| hipStreamSync | 226.6 | 63.3 |
| GPU idle (ms/cycle) | 6.8 | **4.2** |
| GPU busy (ms/cycle) | 42.9 | 44.9 |

- **The syncs and the idle go down, but kernel time goes up.**
  - The fused operator costs 120 us per call against 107 + 12 us (operator + halo).
  - The neighbour lookup in the j/k shell costs more than the halo kernel.
  - Queued no-op iterations after a stop still launch full grids.
- **The operator-level win is therefore small.** The host-side win is real: syncs halve,
  and idle drops by 2.6 ms.

**Weak scaling:** the switch works on 1 rank only. On 2 or more ranks it falls back to Two
(bitwise, CPU e2_dev). So weak scaling is unchanged by this step. The measurements were
queued, see HANDOVER.

## HANDOVER

**Done and gated:**
- The `implicit_krylov_dev` switch (+ `implicit_krylov_dev_halo`), commit below. Default OFF.
- Off is bitwise; on is bitwise on the CPU; NON-CONVERGED 0 everywhere.
- **Restart gate, switch on (K = 2, dev_halo = 1), CPU 1 rank: BITWISE.** N cycles vs
  N/2 + restart + N/2 (`/viper/ptmp2/jinma/launch_0923/rstgate/`, `gate.sh`):
  3-D box vet_sc (N = 30) and 2-D slab Eddington (N = 40). rst data after `<par_end>`
  identical at the restart cycle N and at the end; hst identical except the extra row the
  first leg writes at its end (same time, next dt), as for any restart.
- **GPU on vs off hst (s1c, job 11955221), same binary:** off is repeat-to-repeat
  bitwise. On vs off: mass 2e-16, time 2e-12, dt 3e-11, totE < 6e-12 (row-wise rel.);
  KE and V1mid 1e-10 to 2e-8. The differences are flat in time (no growth) and sit at
  the solver tolerance (`implicit_tol` = 1e-8): device dot products sum in a different
  order, so the inner iteration counts differ (V: 38.28 vs 39.31 per step). Net 2-/3-mom
  columns are round-off noise around 0 and are not meaningful. Accepted as a pass.

**GPU timing, 1 GPU, s1c (job 11955221, ms/cycle, 3 repeats, `tsum.py s1c`):**

| arm | off | K = 2 | K = 2 + dev_halo |
|---|---|---|---|
| E be+L | 36.2, 36.2, 36.0 (36.1) | 34.5, 34.5, 34.9 (34.6) | 35.0, 34.4, 34.8 (34.7) |
| V vet_sc be+L | 43.2, 43.2, 43.4 (43.3) | 42.8, 42.5, 43.0 (42.8) | 42.3, 42.2, 42.2 (**42.2**) |
| Vh vet_sc hesdirk2+L | 59.6, 59.9, 60.0 (59.8) | 59.3, 59.5, 59.5 (59.4) | 58.3, 58.4, 58.3 (**58.4**) |

- NC = 0 in every arm. The switch saves 0.5 to 1.5 ms/cycle; dev_halo = 1 is the better
  operator for vet_sc.

**GPU weak scaling, 2 GPUs, w2 (job 11955222, 0.9 M cells per GPU, switch off since it
is 1-rank only; halo_mpi + overlap; ms/cycle, 2 repeats):**

| arm | ms/cycle |
|---|---|
| E | 43.0, 43.4 (43.2) |
| V | 53.4, 53.4 (53.4) |
| Vh | 73.2, 73.3 (73.3) |
| hydro, 2 GPUs | 21.7, 21.6 (21.6) |
| hydro, 1 GPU | 21.1, 20.9 (21.0) |

**Not done:**
- Step (b), Picard-level fusion. Low value: 36.6 launches per cycle, 0.8 ms idle, 4.7 ms
  of compute-heavy single kernels.
- Step (c), the M1 fences. Only 2 timer fences per solve in ImplicitSolve plus the
  vet.cpp ones (m1-sctb). Hydro holds 43 of the ~49 device syncs.
- Step (d), HIP graphs. Not worth it on 1 GPU: the GPU is busy 91 % with the switch on.
- Weak scaling 4/8 GPUs; the pw2 profile (2 GPUs) is not analysed.
- The multi-rank version of the device scalars. It needs an MPI allreduce per
  reduction, which the pipe already hides.

**Jobs (finished; s1c and w2 are in the tables above, pw2 not yet analysed):**
- **11955221: s1c, 1 GPU.** Off / K = 2 / K = 2 with dev_halo = 1; E, V, Vh; 3 repeats.
  Binary `bin/athena_s1c_gpu`, md5 4be2bf0d.
  Analyse: `python3 /viper/ptmp2/jinma/launch_0923/scripts/tsum.py s1c`.
- **11955222: w2, 2 GPUs weak** (0.9 M cells per GPU). E, V, Vh with halo_mpi + overlap;
  hydro 2 GPUs; hydro 1 GPU.
  Analyse: `python3 /viper/ptmp2/jinma/launch_0923/scripts/tsum.py w2`.
- **11955223: pw2, 2-GPU profile of w2_V.**
  Analyse: `python3 /viper/ptmp2/jinma/launch_0923/scripts/split.py /viper/ptmp2/jinma/launch_0923/runs/pw2/w2_V_1/prof 0 30 10`
  (and rank 1).

**Next step I would take:**
1. s1c: dev_halo = 1 wins for vet_sc (table above); make it the operator of the
   switch.
2. From pw2, find what costs the +8 ms per cycle between 1 and 2 GPUs (halo MPI, MPI
   reductions, SC sweep). That, not the 1-GPU launch count, is where weak scaling is lost.
3. For the 1-GPU radiation time (24 ms of kernels), the levers are kernel time:
   - the stencil layout, with ost slots innermost: 107 us per call is about 1.9 TB/s,
     half of peak;
   - fewer no-op iterations per batch;
   - the SC-sweep kernels (m1-sctb).

## Files

| file | content |
|---|---|
| `scripts/job.sh` | GPU job (arm file + `env_<tag>.sh`; PROF=1 for rocprofv3) |
| `scripts/tsum.py` | ms/cycle, Picard, inner, NC per arm |
| `scripts/split.py` | profile split (runs_4a split.py + HIP API sync counts) |
| `scripts/run.sh`, `sched.py`, `cmp.py`, `cpu_s1*.txt` | CPU gates |
| `scripts/arms_*.txt`, `env_*.sh` | the GPU arms of each job |
