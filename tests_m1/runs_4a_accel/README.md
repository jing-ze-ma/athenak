# runs_4a_accel: cutting the cost of `time_scheme = hesdirk2` + `implicit_vimp`

- **Branch.** `m1-accel`, branched from rt-integration d0c59f7c; worktree
  `/viper/ptmp2/jinma/wt_accel`. Commit: see the end of this file.
- **Build and run tree.** `/viper/ptmp2/jinma/accel_0923`:
  - `base/` = `git archive d0c59f7c`, `new/` = the worktree snapshot;
  - binaries: GPU `bin_base_gpu` (md5 dd3d1ae4) and `bin_v6_gpu` (a007d9b1, the final
    code); CPU `new/build_boxcpu` (e9579410) and `new/build_mpi` (341ae27b);
  - `runs/` GPU runs, `cpu/` 2-D He slab, `rw/` radwave, `gpu_g0/`, `gpu_gk/`.
- **Spatial and coupling setting of every arm:** `implicit_enthalpy = plm`,
  `implicit_vimp = true`, `time2_enth_vel = central` (the default), the fast-path defaults
  of d0c59f7c (stencil, rbgs_fwd, krylov_fuse 3, halo_direct, predictor step, EW 1e-2).
- **GPU environment of every job:** `HSA_XNACK=1`, `HSA_NO_SCRATCH_RECLAIM=1`, apudev.

## Summary

- **Recommended set L:** `implicit_vimp_fold = true`, `implicit_fast_kernels = true`,
  `implicit_one_pass = 4` (or 8) and `implicit_predictor_order = 2`.
- **1 GPU:** hesdirk2 + vimp drops from 75.2 to 60.6 ms/cycle (Eddington) and from 87.4 to
  71.8 (vet_sc).
- **Against be at d0c59f7c:** 1.17x (Eddington) and 1.25x (vet_sc), down from 1.46x and 1.52x.
  With `implicit_one_pass = 8`: 1.16x and 1.21x. On 2 GPUs: 1.20x and 1.24x.
- **Against be with the same levers:** 1.39x (Eddington) and 1.42x (vet_sc). The levers make be
  faster too: 51.6 -> 43.5 and 57.5 -> 50.7 ms/cycle.
- **Gates:**
  - PASS: G0 (CPU and GPU), restart, He slab, stiff boundedness.
  - Borderline: radwave G1 at two stiff corners (1.89 and 1.86), a solver-tolerance effect.

## The levers (all `<rad_m1>`, all default OFF and read only when named)

| key | what it does | exact? |
|---|---|---|
| `implicit_vimp_fold = true` | The implicit_vimp operator part goes into the stored stencil: the x2/x3 +-1 terms into slots 3-6, the +-2 neighbours into new slots 19-24 (`ost` 25 slots). `ImplicitStencilOp` no longer calls `M1VimpRow` (10 separate iw reads per cell per apply). Needs `implicit_op_stencil`. | round-off |
| `implicit_fast_kernels = true` | (i) closure = eddington: the stencil build and the right-hand side skip the off-diagonal Eddington terms. (ii) hesdirk2 + vet_sc with `time2_vet_extrap = false`: the tensor save of `Time2VetExtrapolate` is a plain copy (was a 1.2 ms MDRange reduction). | (ii) bitwise; (i) bitwise on the CPU, round-off on the GPU (below) |
| `implicit_one_pass = N` | A solve that starts from the predictor takes its first Picard pass to the full linear tolerance (no Eisenstat-Walker loosening), and is accepted after that pass when `res_0 qe/(1-qe) < implicit_tol`, `qe = implicit_one_pass_safety` (default 3) x the larger of the last two MEASURED contractions `q = res_1/res_0` of this kind of solve (be / stage 1 / stage 2). Every solve that takes a second pass measures q; at least every N-th solve of a kind is made to. The q state travels in the restart file (new marked header `M1ONEP01`). | within solver tolerance |
| `implicit_predictor_order = 2` | The predictor extrapolates the increment RATE linearly in time: `x0 = start + dt (g1 + h dt1)`, `h = (g1 - g2)/dt2`, stored per cell as ipred channels 3-4 (E, T); ipred / ipred2 get 5 channels and travel in the restart file as before. Under hesdirk2 a backward-Euler step (first step, fallback) leaves ipred alone. | within solver tolerance |
| `implicit_op_split_red = true` | stencil apply as par_for + a separate reduction. DIAGNOSTIC: no gain (68.8 vs 69.2 ms, one run), not recommended. | round-off |

**Why (i) is not bitwise on the GPU.** `M1EddOff = 0.5 (3 chi - 1) n_a n_b` with chi = 1/3 is
exactly 0 on the CPU, but the GPU contracts `3 chi - 1` into an FMA and gets -5.6e-17.
So at d0c59f7c the GPU stencil of an Eddington run carries non-zero edge coefficients of
relative size 1e-17: `st_edges` is true and every apply reads all 19 slots. That is why (i)
saves 7.4 ms/cycle on the GPU. The GPU check (job 11954110, `gpu_gk/`, 40 cycles): vet_sc
hst and rt_profile are identical; Eddington KE differs by 2.4e-10 relative and the bin data
by <= 2.6e-7 relative (a near-zero field).

**Why one_pass measures q instead of assuming it.** The first version accepted after pass 0
with a FIXED assumed q1 = 5e-4 (the box has q ~1e-4). In the stiff radwave corners the real
q is ~2e-3 with the d0c59f7c solver, and the fixed version raised e1024 there by 50-70 %.

**What was measured to decide the levers** (1 GPU, 3-D box, bin_v2..v5; `runs/S_*`, `T_*`, `U_*`):
- Profile of d0c59f7c (job 11952845, `runs/P0_*/split.txt`, cycles 5-40):
  - `ImplicitStencilOp` takes 120 us per call with vimp. The stencil alone takes 76 us.
  - PCR takes 36.5 us per half-sweep.
  - hesdirk2 has 4.4 Picard passes per cycle against 2.4 under be. Each pass costs about 3 ms of assembly, vimp build and T solve.
- Steady-state Picard logs (`runs/S_Eh2_1`, `runs/S2_Eh2ew0_1`):
  - With EW 1e-2, pass 0 takes 3-5 inner iterations and pass 1 takes 22-26. Pass 1 exists only to finish the inexact pass 0.
  - With pass 0 exact (EW off), `res_1 = 2-3e-10` against `res_0 = 2-4e-6`. So the Picard contraction is q ~ 1e-4 in the box (vet_sc 1.3e-4, be 3.5e-4).
- Options that did not help (1 run each):
  - `implicit_lin_ew_max` 3e-3 / 3e-2 / 0.1;
  - `implicit_precond = rbgs`;
  - `implicit_pcr_team = 64`;
  - `implicit_precond_float` (-0.9 ms, within noise).

## GPU cost (3-D He box 84x104x104, 4 blocks, 120 cycles, ms/cycle over cycles 20-120)

Same binary (bin_v6_gpu, md5 a007d9b1), arms interleaved, repeat 2 reversed, 3 repeats.
L = `implicit_vimp_fold + implicit_fast_kernels + implicit_one_pass=4 + implicit_predictor_order=2`;
L8 = the same with `implicit_one_pass=8`. Picard = passes per solve, inner = Krylov iterations
per solve (2 solves per step under hesdirk2). NON-CONVERGED 0 in every arm.

### 1 GPU (jobs 11954111 Eddington, 11954112 vet_sc full)

| arm | Eddington ms/cycle (a, b, c) | mean | Picard | inner | vet_sc ms/cycle (a, b, c) | mean | Picard | inner |
|---|---|---|---|---|---|---|---|---|
| be | 52.4, 51.3, 51.2 | 51.6 | 2.125 | 38.3 | 57.2, 57.1, 58.1 | 57.5 | 2.733 | 39.9 |
| hesdirk2 | 74.8, 75.9, 75.0 | 75.2 | 2.067 | 29.8 | 86.5, 87.8, 87.8 | 87.4 | 2.623 | 31.5 |
| be + L | 43.1, 44.3, 43.0 | 43.5 | 1.350 | 35.8 | 49.8, 51.8, 50.5 | 50.7 | 1.358 | 39.3 |
| hesdirk2 + L | 61.0, 60.4, 60.4 | 60.6 | 1.318 | 25.7 | 70.5, 72.3, 72.5 | 71.8 | 1.339 | 27.3 |
| hesdirk2 + L8 | 61.0, 59.3, 59.4 | 59.9 | 1.201 | 25.4 | 69.2, 70.7, 69.1 | 69.7 | 1.218 | 26.9 |

| ratio | Eddington | vet_sc |
|---|---|---|
| hesdirk2 / be, d0c59f7c | 1.46 | 1.52 |
| hesdirk2 + L / be (d0c59f7c) | **1.17** | **1.25** |
| hesdirk2 + L8 / be (d0c59f7c) | **1.16** | **1.21** |
| hesdirk2 + L / be + L | 1.39 | 1.42 |

### Per lever, cumulative (1 GPU, jobs 11954113 Eddington, 11954114 vet_sc)

| hesdirk2 + | Eddington (a, b, c) | mean | step | vet_sc (a, b, c) | mean | step |
|---|---|---|---|---|---|---|
| nothing | 77.4, 77.0, 76.4 | 76.9 | | 88.6, 89.0, 89.5 | 89.0 | |
| vimp_fold | 75.1, 76.3, 75.9 | 75.8 | -1.1 | 86.1, 87.4, 87.6 | 87.0 | -2.0 |
| + fast_kernels | 68.8, 68.4, 67.9 | 68.4 | -7.4 | 84.2, 84.9, 86.2 | 85.1 | -1.9 |
| + one_pass=4 | 65.0, 65.1, 65.1 | 65.1 | -3.3 | 80.3, 74.9, 75.5 | 76.9 | -8.2 |
| + predictor_order=2 | 60.5, 61.2, 59.7 | 60.5 | -4.6 | 70.9, 71.8, 75.8 | 72.8 | -4.1 |

(vet_sc repeats are noisier: one_pass a = 80.3, predictor b/c spread 5 ms.)

### 2 GPUs (jobs 11954115 Eddington, 11954116 vet_sc; 2 ranks, 2 blocks each)

| arm | Eddington ms/cycle (a, b, c) | mean | vet_sc ms/cycle (a, b, c) | mean |
|---|---|---|---|---|
| be | 50.9, 51.4, 51.5 | 51.3 | 58.4, 58.3, 58.6 | 58.4 |
| hesdirk2 | 77.8, 76.5, 77.4 | 77.2 | 88.4, 89.4, 88.1 | 88.6 |
| be + L | 42.9, 43.3, 42.5 | 42.9 | 50.7, 50.8, 51.3 | 50.9 |
| hesdirk2 + L | 61.8, 61.6, 61.7 | 61.7 | 75.3, 71.0, 70.8 | 72.4 |

Ratios, 2 GPUs:

| ratio | Eddington | vet_sc |
|---|---|---|
| hesdirk2 / be, d0c59f7c | 1.50 | 1.52 |
| hesdirk2 + L / be (d0c59f7c) | **1.20** | **1.24** |
| hesdirk2 + L / be + L | 1.44 | 1.42 |

The box is too small for 2 GPUs to help: be takes 51.3 ms against 51.6 on one GPU. The
halo code was left alone, as the brief asks. NON-CONVERGED 0 in every arm.

## Gates

| gate | result | evidence |
|---|---|---|
| G0 CPU, switches off, bitwise vs d0c59f7c | PASS: radwave 6 cases (4 be + 2 hesdirk2, Eddington and vet_sc), 204 tab/hst files identical; 2-D He slab be / hesdirk2 / hesdirk2+vet_sc 200 s, 36 bin/hst/rst files identical | `rw/g0_{base,new}`, `cpu/g0_*_{base,new}` |
| G0 GPU, switches off, bitwise vs d0c59f7c | PASS with bin_v5_gpu (md5 75a0e8e5, the final code apart from the fast_kernels vet-copy path and the od_skip -> fast_kernels rename): 3-D box, 40 cycles, be / hesdirk2 / hesdirk2+vet_sc, 39 bin/hst/rt_profile files identical (job 11953988). Re-checked with the final bin_v6_gpu (md5 a007d9b1, job 11954212): 39 files identical, 0 DIFF | `gpu_g0/`, `log.out.11953988`, `log.out.11954212` |
| radwave G1, x1 64x4, 12 cases x nt 128..2048, tol 1e-11 | **BORDERLINE**, see below: the median order is >= 1.90 in 10 of 12 cases (Eddington and vet_sc alike). At (10,1e3) it is 1.89 and at (100,1e3) 1.86, against 2.00 / 1.99 at d0c59f7c. That difference is a solver-tolerance effect: at tol 1e-12 both runs give 2.00 / 1.99. d0c59f7c itself has 1.89 at (100,10) | `RESULTS_g1_edd.txt`, `RESULTS_g1_vet.txt`, `RESULTS_tol_check.txt` |
| stiff P=100 tau_lambda=1e5 (40 periods) and the vimp stiff arms (10 periods) | PASS (bounded): P=100 tau=1e5 end amplitudes are the same as d0c59f7c (0.719 / 0.800 / 0.801 against exact 0.799), N32 nt24 NC 1 in both; all P=100 arms NC 0, 0 fallbacks. P=1 tau=1e5 at the unreachable tol 1e-11 (the known stall): NC 531 against 164 at d0c59f7c. At tol 1e-10 and 1e-9, both have NC 0, 0 fallbacks and end 0.795 (exact 0.799) | `RESULTS_stiff.txt` |
| He slab 2-D, levers on | PASS: NON-CONVERGED 0 and 0 stage fallbacks in every run. 200 s runs: Eddington, vet_sc, 2 ranks and be; one pass on 1853 of 2481 solves, Picard 1.25/solve. 1000 s runs: Eddington and vet_sc. The forced fallback (`time2_dbg_fail=50`) takes 1 fallback and runs through. At t=1000, KE1 = 1.0858e26 with and without levers; KE2 = 2.3307e23 in both. At 200 s the hst differs from the lever-off run by <= 3e-7 relative in KE and V1max | `cpu/Z_*` (final binary), `cpu/M_*_1000` and `cpu/X_h2_1000` (binary before the od_skip -> fast_kernels rename, same code path) |
| restart, levers on | PASS: bitwise | below |

### Radwave G1 (A = d0c59f7c, B = levers: fast_kernels + one_pass=4 + predictor_order=2)

`implicit_vimp_fold` cannot run here: the radwave grid is x1-periodic and the stencil path
(which the fold needs) refuses a periodic x1. It is covered by the slab and box runs instead.

Entries: median successive-difference order (lowest single order) e128 e1024 (`med_table.py`,
error against the Richardson value of the two finest nt, as `time_table.py`).

```
Eddington
     P     tau | A: med p (min) e128 e1024      | B: med p (min) e128 e1024     
   0.1     0.1 | 2.00 (1.99) 4.98e-04 7.82e-06 | 2.00 (1.99) 4.98e-04 7.82e-06
   0.1      10 | 1.99 (1.98) 4.14e-04 6.65e-06 | 1.99 (1.98) 4.14e-04 6.65e-06
   0.1    1000 | 1.99 (1.98) 4.21e-04 6.72e-06 | 1.99 (1.98) 4.21e-04 6.74e-06
     1     0.1 | 1.99 (1.98) 4.42e-04 7.09e-06 | 1.99 (1.98) 4.42e-04 7.09e-06
     1      10 | 1.99 (1.99) 3.86e-04 6.15e-06 | 1.99 (1.99) 3.86e-04 6.15e-06
     1    1000 | 1.92 (1.85) 4.01e-04 7.43e-06 | 1.90 (1.88) 4.01e-04 7.24e-06
    10     0.1 | 1.97 (1.96) 3.91e-04 6.48e-06 | 1.97 (1.96) 3.91e-04 6.48e-06
    10      10 | 1.99 (1.98) 4.02e-04 6.44e-06 | 1.99 (1.98) 4.02e-04 6.44e-06
    10    1000 | 2.00 (1.99) 1.26e-03 1.67e-05 | 1.89 (1.86) 1.26e-03 2.35e-05
   100     0.1 | 1.99 (1.99) 3.81e-04 6.06e-06 | 1.99 (1.99) 3.81e-04 6.06e-06
   100      10 | 1.89 (1.79) 4.97e-04 8.35e-06 | 2.03 (1.94) 5.27e-04 7.55e-06
   100    1000 | 1.99 (1.96) 6.94e-04 9.52e-06 | 1.86 (1.70) 6.95e-04 1.52e-05
vet_sc (full)
     P     tau | A: med p (min) e128 e1024      | B: med p (min) e128 e1024     
   0.1     0.1 | 1.95 (1.91) 5.11e-04 8.86e-06 | 1.95 (1.91) 5.11e-04 8.86e-06
   0.1      10 | 1.92 (1.87) 4.33e-04 7.95e-06 | 1.92 (1.87) 4.33e-04 7.95e-06
   0.1    1000 | 1.99 (1.98) 4.21e-04 6.69e-06 | 1.99 (1.98) 4.21e-04 6.71e-06
     1     0.1 | 1.98 (1.98) 4.43e-04 7.21e-06 | 1.98 (1.98) 4.43e-04 7.21e-06
     1      10 | 1.96 (1.94) 3.96e-04 6.68e-06 | 1.96 (1.94) 3.96e-04 6.68e-06
     1    1000 | 1.93 (1.85) 4.01e-04 7.38e-06 | 1.91 (1.88) 4.00e-04 7.15e-06
    10     0.1 | 1.97 (1.95) 3.92e-04 6.51e-06 | 1.97 (1.95) 3.92e-04 6.51e-06
    10      10 | 1.96 (1.95) 4.11e-04 6.91e-06 | 1.97 (1.95) 4.11e-04 6.91e-06
    10    1000 | 2.00 (2.00) 1.25e-03 1.65e-05 | 1.89 (1.87) 1.26e-03 2.34e-05
   100     0.1 | 1.98 (1.98) 3.79e-04 6.11e-06 | 1.98 (1.98) 3.79e-04 6.11e-06
   100      10 | 1.95 (1.77) 5.33e-04 9.11e-06 | 1.99 (1.90) 5.79e-04 8.75e-06
   100    1000 | 1.99 (1.97) 6.93e-04 9.40e-06 | 1.86 (1.69) 6.95e-04 1.50e-05
```

- **Where it changes.** Ten cases have the same e128 and e1024 to 3 digits. The changes sit at
  (10,1e3) and (100,1e3): e1024 rises from 1.67e-5 to 2.35e-5 and from 9.5e-6 to 1.52e-5, and the
  median order drops from 2.00 / 1.99 to 1.89 / 1.86. e128 is unchanged.
- **Each lever alone does it.** `implicit_predictor_order=2` alone and `implicit_one_pass=4` alone
  each give the same change; `implicit_fast_kernels` alone is identical to A (`rw/iso`).
- **It is a solver-tolerance effect** (`RESULTS_tol_check.txt`). At tol 1e-12 the two agree:
  e1024 1.67e-5 vs 1.68e-5 and 9.5e-6 vs 1.06e-5, with medians 2.00 / 1.99 in both. With a
  wave of amplitude 1e-5 and tol 1e-11 relative to E, a coherent per-step solver error of
  ~1e-11 amounts to ~1e-5 of the wave after 1000+ steps. Which point inside the tolerance
  ball the solve stops at depends on the starting iterate and the pass count.
- **Status against the brief's ">= ~1.9 in the 12 cases": borderline.** 1.89 and 1.86 at the
  stiffest two corners at the gate tolerance, 2.00 / 1.99 at tol 1e-12.

### Restart (CPU, 2-D slab, levers on)

Restarted from the t=100 file and run to 200 s, compared with the continuous run:
- **Runs.** hesdirk2 Eddington (`cpu/R_h2`), hesdirk2 vet_sc (`cpu/R_hv`) and be (`cpu/R_be`).
- **Result.** In all three, the t=200 bin (hydro_w, m1) and rst files are identical, and all
  100 overlapping hst lines are identical.
- **The restart carries the lever state.** The restarted hesdirk2 runs took 0 BE steps, so the
  slope and ipred2 (5 channels) were read. The one-pass q state was read from M1ONEP01.

Final binary: `RESULTS_restart.txt` (`cpu/Z_*` continuous, `cpu/RZ_*` restarted).

## Files

| file | content |
|---|---|
| `run_radwave.py` | the runs_3x runner plus the key `rm_<par>=<v>` -> `<rad_m1>/<par>` |
| `time_table.py`, `med_table.py` | G1 order tables (`RW_ALL=1` also takes tol/amp tagged runs) |
| `stab_table.py`, `stiff_table.py` | stiff probes |
| `lists/` | g0h (G0), g1_*_L (G1 with levers), stiff_L / stiff_0, iso*, p1tau |
| `inp/` | box and slab inputs; `*_x` name every lever at its default (a command-line override of a parameter not in the file is fatal) |
| `scripts/` | GPU jobs (`sweep.sh` 1 GPU, `sweep2.sh` 2 GPUs, `prof.sh`, `g0_gpu.sh`, `gk_gpu.sh`), CPU runners, build scripts, `summarize.py`, `split.py`, `split2.py`, `hstcmp.py`, `bincmp.py` |
| `RESULTS_*.txt` | the tables above |
