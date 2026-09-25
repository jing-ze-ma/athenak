# runs_5r_fast4: making implicit M1 + VET cheaper (m1-fast4)

- **Date:** 2026-09-25, viper.
- **Branch:** `m1-fast4`, from rt-integration e91d6014, with rt-integration (m1-sp-order2b) merged in at 364e5899.
  - Worktree: `/viper/ptmp2/jinma/wt_fast4`.
- **Run tree:** `/viper/ptmp2/jinma/fast4_0925`. Build dirs are deleted.
  - `bin/`: binaries. The final one is `m8` = 76241f53 (+ lint-only 37cb379f).
    - md5 box GPU 1af9e1c8..., none GPU 431cadfe..., box CPU 0aebdd73..., none CPU 4025844e....
  - `runs/p*`: GPU jobs.
  - `cpu/`: CPU gates.
  - `rw/`: radwave.
  - `st/`: the runs_5q space-time set.
- **Cases** (as in prof_m1_0924 / runs_5p_coarse2):
  - **box** = He box `inp/box.athinput`:
    - vet_sc, full tensor, `implicit_precond = mg`, `implicit_mg_levels = 3`, hesdirk2;
    - 84x104x104 on 1 GPU; 84x208x104 on 2 GPUs (weak scaling).
  - **wedge** = sp He wedge `inp/wedge.athinput`:
    - 96x128x128, vet_col, `mg_gc` levels 1, with the m1-sp-order2b sp defaults;
    - 2 GPUs is strong scaling.
- **GPU timing:**
  - apudev MI300A, `HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1`.
  - Same binaries, interleaved, 3 reps, 60 cycles, ms/cycle over cycles 10-50 (`tsum.py`).
- **Fenced timers:** `<rad_m1>/implicit_timers = N` (diagnostic, read only when named), summarised by `tmr.py`. They split the cycle into:
  - radiation = the M1 stage chain (closure limits, opacity, the implicit solve including the tensor build, the rad-triggered hydro bvals + ConsToPrim, and the M1 bvals);
  - rest = hydro and everything else.
- **Profiles:** rocprofv3 `--kernel-trace --kokkos-trace`, per-label summary in `lsum.py`.

## 1. Where the time goes (job 11978969, t1 = e91d6014 + timers; wedge on the pre-merge sp defaults)

| ms/cycle | box 1 GPU | box 2 GPUs | wedge 1 GPU | wedge 2 GPUs |
|---|---|---|---|---|
| total | 47.5 | 53.9 | 31.4 | 19.1 |
| radiation | 33.0 | 38.9 | 24.5 | 14.8 |
| - Krylov (BiCGStab + mg) | 13.7 (17.2 it/cycle, 0.80 ms/it) | 16.3 | 5.3 (3.1 it, 1.7 ms/it) | 3.2 |
| - tensor (SC sweep / vet_col) | 4.2 | 6.0 | 7.7 | 4.0 |
| - pass setup (transverse terms, src, asm, vimp, stencil) | 4.0 | 4.4 | 3.2 | 2.0 |
| - solve pre + end + pass post | 6.6 | 6.7 | 6.2 | 3.8 |
| - rad-triggered hydro bvals + ConsToPrim | 2.5 | 2.8 | 0.9 | 0.8 |
| hydro and the rest | 14.5 | 15.1 | 7.2 | 4.6 |
| radiation / hydro | 2.28 | 2.57 | 3.41 | 3.21 |

- **Solves:**
  - hesdirk2 = 2 stage solves per step.
  - Picard passes per solve: 1.28 (stage 1) and 1.12 (stage 2) on the box; 1.12 on the wedge.
- **Summary.** The cost is a fixed per-solve price (setup kernels plus the tensor build), not Krylov iterations. lin_tol stays 1e-10.
- **Post-merge wedge** (sp defaults; job 11979229):
  - 33.4 ms/cycle, radiation/hydro 3.69;
  - `m1_vcol_team` 9.3 ms, 28 % of the cycle.

## 2. Levers

### Kept: accuracy-neutral kernel work (76241f53; all bitwise)

The cause: Kokkos' HIP RangePolicy is compiled for 1024 threads per block, which caps a kernel at 128 VGPRs. The heavy EOS/opacity kernels spilled.

- `par_for_lb` (`rad_m1_parfor.hpp`, `LaunchBounds<256,1>`) on:
  - `m1_impl_i1`, `m1_vimp_j`, `m1_opacity`, `m1_t2_vsf`, `m1_t2_vcp`.
- `M1_INL` (an always_inline lambda body) + `par_for_lb` on `m1_impl_asm` and `m1_impl_tsolve`. Their bodies were *called* from the Kokkos wrapper, with 176-700 B stack frames per lane.
- Flat 1-D reductions (max/min, order-independent) for `m1_impl_eck` and `m1_t2_adm`.
- One set of logarithms and bisections for both opacity tables.

Per kernel, us/call before -> after (profiles runs/p6/prof/*_m6 -> runs/p8/prof/*_m8):

| kernel | box | wedge |
|---|---|---|
| m1_impl_i1 | 371 -> 243 | 248 -> 123 |
| m1_opacity | 538 -> 348 | 249 -> 70 |
| m1_vimp_j | 527 -> 312 | - |
| m1_impl_eck (diagnostic) | 742 -> 116 | - |
| m1_t2_adm | 157 -> 41 | 366 -> 150 |
| m1_impl_asm | 186 -> 149 | 349 -> 304 |
| m1_impl_tsolve | 165 -> 226 | 808 -> 609 |
| m1_t2_vcp | - | 301 -> 232 |
| unchanged: m1_impl_src, m1_impl_stb, m1_t2_vsf, m1_vcol_team | 210, 519, 551, - | 308, 455, -, 9306 |

### Tried and not kept

| lever | result |
|---|---|
| forced inline of `m1_impl_src` | slower: 210 -> 304 us (box), 307 -> 408 us (wedge); reverted |
| forced inline of the EOS table evaluation (`EvalFromLogs`) | scratch of the EOS kernels unchanged; reverted |
| stencil-build loop unrolling / constant slots | the 176 B scratch array remains; reverted |
| vet_col team kernel: team size 64/128/256, chunk 8/32 (job 11979476), unrolled moment loops + coalesced transposed tables (job 11979711) | 9.3-9.6 ms either way (chunk 32: +13 ms); reverted |
| `implicit_gas_newton = true` on the wedge (round-off) | +0.1 ms; not recommended |
| `implicit_krylov_dev = 2` | refuses mg / mg_gc |

**`vet_sc_every = 2`** (tensor cadence with linear extrapolation in time):
- Committed, default off.
- CPU box accuracy vs a tight reference: identical to the default.
- `vet_sc_every = 4` is 45x worse in Frad in the seeded transient.
- On GPU it is **slower**: 44.7 -> 60.0 ms/cycle on the box. The extrapolated tensor spoils the predictor: pass-0 residual 2e-6 -> 1.6e-4, so one-pass acceptance never fires (2 passes per solve).
- Not recommended.

**Multi-rate radiation** (`implicit_mr_every = k`, `implicit_mr_theta`, `implicit_mr_peq`):
- Committed, default off. Per the user decision (no accuracy sacrifice) it is not recommended.
- What it is:
  - an operator-split SDIRK2 over k hydro steps, placed at the window centres (Strang);
  - window state and cadence state are in the restart file (header M1MRWIN1).
- Timing, ms/cycle:

  | case | default | k=2 | k=4 |
  |---|---|---|---|
  | box 1 GPU | 46.7 | 41.2 | 29.2 |
  | wedge 1 GPU | 33.4 | 32.5 | 18.1 |

- It is **first order**:
  - runs_5q time refinement (`st/RESULTS_mr.txt`): pulse_u_T E/F order 1.0-1.4; rsw_u_T F1 order 1.0-1.1. Causes: SDIRK2's first stage is at c = g, and the splitting in the coupled regime.
  - radwave P/tau = 1/10: first-order artificial damping (`rw/`, `scripts/rwtab.py`).

**`implicit_eos_cache_check_every = 10`** in box inputs:
- Bitwise state; the check is a diagnostic.
- Worth 1.5 ms before the eck fix, about 0.2 ms after it.

## 3. Final GPU timings (job 11979774, m4 = base code path vs m8, same job, interleaved, 3 reps)

| ms/cycle | box 1 GPU | box 2 GPUs | wedge 1 GPU | wedge 2 GPUs |
|---|---|---|---|---|
| before (m4 defaults = rt-integration bitwise) | 47.55, 47.46, 47.31 | 53.58, 53.66, 53.73 | 33.61, 33.87, 33.68 | 20.56, 20.38, 20.40 |
| after (m8) | 44.55, 44.89, 44.79 | 50.90, 51.24, 51.25 | 31.96, 31.96, 31.94 | 19.51, 19.47, 19.52 |
| radiation ms/cycle (timers), before -> after | 33.1 -> 30.4 | 38.9 -> 36.2 | 26.7 -> 24.9 | 16.1 -> 15.2 |
| radiation / hydro, before -> after | 2.26 -> 2.08 | 2.57 -> 2.40 | 3.74 -> 3.47 | 3.43 -> 3.24 |

Scaling from 1 to 2 GPUs (after):
- box, weak scaling: 0.875 (before 0.884);
- wedge, strong scaling: 0.82 (before 0.82).

## 4. Gates (m8)

- **CPU default path bitwise:** `gates.py` arms box_a1, box_m2, slab_a1, slab_m2, loose and tight, b0 (e91d6014) vs m8. Restart files and history are bitwise (`cpu/m8.out`).
- **Restart bitwise:** 14 cycles vs 7 + restart + 7 (`rst_mr.sh`).
  - m8: default.
  - m6: `vet_sc_every = 2`, `implicit_mr_every = 4` + theta.
- **runs_5q space-time set** (pulse_u_T, rsw_u_T, mvr_u_T, pulse_u, atm_u; CPU): m8 = m4, all 16 L1 rows identical (`st/RESULTS_f8.txt`). The orders stay 1.89-2.06.
- **GPU:** final restart files are bitwise at 1 and 2 GPUs (runs/p6, runs/p8):
  - box: t1 = m4 = m6 = m7b = m7c = m8;
  - wedge on the post-merge defaults: m4 = m6 = m7b = m7c = m8.
- **Non-converged:** 0 in every run.

## 5. 4-8 GPU scripts (not run)

    cd /viper/ptmp2/jinma/fast4_0925
    sbatch -p apudev -N 1 --ntasks-per-node=1 --gres=gpu:1 scripts/job_scale4.sh 1
    sbatch -p apudev -N 1 scripts/job_scale4.sh 2
    sbatch -p apu -N 2 scripts/job_scale4.sh 4
    sbatch -p apu -N 4 scripts/job_scale4.sh 8
    python3 scripts/tsum.py runs/s<NP>_m8 10 30

## 6. What is left (by size, per cycle on 1 GPU)

- **Wedge:** `m1_vcol_team` 9.3 ms. It is throughput-bound with an unknown limiter; access pattern, team size and chunk do not move it. The next step is hardware counters.
- **Box:** Krylov 13.7 ms, i.e. about 8.5 iterations per solve at 0.8 ms/iteration with mg 3.
- **Both:**
  - `m1_impl_stb`: 0.5 ms/call, 176 B scratch array.
  - `m1_impl_src`: call frame.
  - `m1_t2_vsf`: EOS inversion call.
  - general-EOS ConsToPrim: 0.7 ms/call, 8 calls per cycle on the box, 2 of them triggered by the radiation.
