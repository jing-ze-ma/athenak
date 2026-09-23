# Fewer Picard passes and inner iterations in the implicit M1 solve

Date 2026-09-23, viper. Base: HEAD d14f96de (rt-integration). NOT COMMITTED.
- Patch: `bench/m1_picard_0923/picard.patch` (md5 b10a9fa4...). It passes `git apply --check`
  on HEAD. Files: `src/rad_m1/rad_m1.{hpp,cpp}` and `src/rad_m1/rad_m1_implicit.cpp` only.
- Snapshots: `bench/m1_picard_0923/{base,new}`. Binaries: `bench/m1_picard_0923/bin/`.
  - `base` = HEAD.
  - `new4` = the patch. Its GPU md5 is 22e51fa2...
  - After new4 was built, two error-message strings in the patch were re-wrapped to fit
    90 columns. Nothing else changed.
- CPU runs are in `cpu/`, GPU runs in `runs/`. `stats.py` produces the CPU statistics and
  `summarize.py` the GPU timing.

## 1. How the loop converges (measured, `implicit_picard_log`)

New diagnostic: `implicit_picard_log = N` prints one line per pass for the first N
solves. Each line gives the E and T parts of the Picard residual, the x1 index of the worst
cell, the transverse change lresid, the inner iteration count and the inner starting
residual. Sources: `cpu/log_edd`, `cpu/log_vet`, and GPU job 11944947 (`runs/LOG_E`,
`runs/LOG_V`).

GPU, 3-D box, Eddington closure, a typical step (step 16):

| pass | Picard residual | lresid | inner its |
|---|---|---|---|
| 0 | 1.2e-4 | - | 68-79 |
| 1 | 3.0e-8 | 1.3e-3 | 30-34 |
| 2 | 1e-12 to 3e-12 | 7.5e-7 | 0-1 |
| 3 | 3e-16 | 1.5e-10 to 3e-10 | 0 |
| 4 (when needed) | 2e-16 | ~0 | 0 |

The 2-D CPU slab and vet_sc full show the same pattern.

**The convergence is quadratic.** The residual goes 1e-4, then 3e-8, then 1e-12, i.e.
C·r² with C of about 2. The gas Newton (Schur elimination with the consistent cached c_v)
is an exact Newton step for the coupled (E, T) system.

**Lagged terms, checked one by one:**
- Closure and tensor: frozen for the step by `closure_lag = step`. For vet_sc it is fixed
  per step by construction.
- Opacities: `opac_update = false`.
- EOS cache: frozen density, with the exact derivative of the cached e(T).
- The transverse and off-diagonal Eddington terms are already inside the BiCGStab operator.
- Still lagged: G0 in the v·G0 face term, and de0 through F. Both are O(v/c). Their
  estimated contraction is below about 1e-4 per pass.

No linear term that is lagged limits the rate. **The mode-3 lesson (put a lagged linear term
into the Jacobian) has nothing left to fix here.**

**Why the mean is 4.6 to 4.9 passes.** Only passes 0 and 1 do Krylov work. The rest come
from how convergence is tested:
- The step always spends one confirming pass: pass 1 lands at 3e-8, just above
  `implicit_tol` = 1e-8.
- `lresid` must also fall below `lin_tol` = 1e-10, relative to max E. It measures the change
  of U(E) = TRHS between passes, so it lags the Picard test by one pass. It is amplified by
  (c dt/dx2)² in the thin top cell (i = 86).
- Under bicgstab this test is redundant: the transverse coupling is in the operator and is
  solved to lin_tol on every pass.

**The inner iterations are where the cost is.** On the GPU, about 110 inner iterations per
step go into passes 0 and 1.

## 2. What the patch adds (every key is off by default; the default path is bitwise HEAD)

| key | effect |
|---|---|
| `implicit_picard_log = N` | the per-pass diagnostic above |
| `implicit_lres_test = false` | drop the lagged lresid test (only under bicgstab) |
| `implicit_conv_est = true` | also stop when the bound on what remains, r·q/(1-q) with q = r_k/r_{k-1} < 1/2, is below `implicit_tol`. This removes the confirming pass. |
| `implicit_predictor = step` | start the loop from E^n + (dt/dt_prev)·ΔE_prev and T^n + (dt/dt_prev)·ΔT_prev |
| `implicit_lin_ew_max = η` | Eisenstat-Walker (choice 2) inner tolerance, see below |
| `implicit_lin_cnorm = ε` | inner test on the per-cell norm max\|r_i\|/((1+SRCB_i)·E_i), see the note below |

Notes on each option:
- **Predictor.**
  - The increments are stored in a new `ipred(m,3,…)` array, allocated only when the
    predictor is on.
  - It is used only with a closure that does not read the iterate: Eddington, vet_sc or
    tau. Only the starting point moves: E^n, e^n, the EOS-cache window and all the scales
    stay those of the step.
  - It runs after `VetShortChar()`. The formal solution reads T^n from `M1_IW_TP` as its
    source; an earlier version that predicted T before it gave 5 % KE_2 differences in the
    vet_sc slab (`cpu/gv_D`).
  - The predicted E is exchanged into the ghost cells before the loop starts.
  - The first step after a start or a restart is cold.
- **Eisenstat-Walker.**
  - The tolerance is max(lin_tol·max|b|, η_k·max|r0|), with η_0 = η and
    η_k = min(η, 0.9·(|r0_k|/|r0_{k-1}|)²), plus the 0.9·η_{k-1}² safeguard.
  - It needs `implicit_bcg_sync ≥ 1`, and is wired into `ImplicitBiCGStabFused` only.
  - The inexact pass-0 error is removed in pass 1. That pass stays at lin_tol, because its
    r0 is already small.
- **`implicit_lin_cnorm`: the inner norm question.**
  - A per-cell norm that honestly bounds the error of E must weight by the row excess
    1 + SRCB, not by the diagonal: the transport rows sum to zero. For an M-matrix,
    |δE_i| ≤ max_j |r_j|/s_j.
  - With that weighting, cnorm = 1e-9 needs more inner iterations than the default
    (CPU: 67 358 against 49 076 over 1241 steps). cnorm = 1e-8 is about the same as the
    default (51 282).
  - A first version weighted by the diagonal A_ii and cut the inner iterations 5x, but it
    is not an error bound, so it was discarded.
  - **Conclusion: the max|b| test at 1e-10 is not over-solving.** It is equivalent to about
    2.5e-8 per cell in the thin top. The option is kept for experiments; it is not
    recommended.

Kernels are device-only, have no `this` capture and add no large functors. cpplint with the
repo CPPLINT.cfg reports 10 errors on the 3 files before the patch and 10 after. No line is
over 90 columns.

## 3. Gates

**Gate 1: the default path is byte-identical to HEAD.**
- CPU, seeded 2-D slab (`inp/slab2d.athinput` + `problem/vpert=1e-2`, `tlim = 200`,
  1241 steps):
  - Eddington (`cpu/base_edd` vs `cpu/g4_def`): hydro.hst and user.hst are identical, and
    the payloads of all 6 bin files are identical.
  - vet_sc full (`gv_base` vs `gv4_def`): hst identical.
  - 2 ranks (`m2*`, meshblock nx2 = 16), Eddington and vet_sc: hst identical.
- GPU: the `ref` arms of base and new4 have identical hst files (E1 and V1).
- The GPU runs are deterministic: repeat a equals repeat b bitwise in every arm.

**Gate 2: CPU, 2-D slab, 200 s statistics.**
- Values are relative differences from HEAD. `tol10` is HEAD run with `implicit_tol = 1e-10`;
  the ctl rows are HEAD with `rad_flux_inner` changed by ±1 ulp.
- NON-CONVERGED is 0 in every arm.

| run | F1top/Fin mean (t>100) | KE_1 end | KE_2 end | dt mean | totE end |
|---|---|---|---|---|---|
| **Eddington** (HEAD: 1.000002, 1.0852e26, 1.728e23, 0.16120 s) | | | | | |
| ctl +1 ulp / -1 ulp | 3e-12 / 2e-12 | 1.3e-9 / 4.7e-9 | 2.5e-8 / 3.6e-8 | 1e-13 | 1e-12 |
| tol10 | 1e-12 | 2.2e-9 | 3.7e-8 | 3e-12 | 1e-12 |
| B (lres off + conv_est) | 6e-11 | 2.3e-8 | 5.6e-8 | 4e-12 | 1e-11 |
| D (B + predictor) | 1.3e-10 | 4.1e-9 | 2.0e-8 | 4e-12 | 7e-12 |
| D2 (D + EW 1e-2) | 1.3e-10 | 1.2e-8 | 3.5e-8 | 4e-12 | 8e-12 |
| **vet_sc full** (HEAD: 1.000004, 1.7783e26, 1.735e23, 0.16053 s) | | | | | |
| ctl +1 ulp | 6e-12 | 1.3e-8 | 3e-10 | 3e-11 | 3e-12 |
| tol10 | 3e-14 | 1.3e-8 | 5.1e-8 | 4e-11 | 8e-14 |
| B | 2e-11 | 3.2e-8 | 3.4e-8 | 7e-11 | 3e-12 |
| D | 5e-11 | 1.2e-7 | 3.2e-8 | 6e-11 | 3e-11 |
| D2 | 7e-11 | 1.1e-7 | 1.3e-8 | 7e-11 | 2e-11 |

- The changes to the statistics are at the solver-tolerance level: at most 1.3e-10 on
  F1top/Fin and 1e-11 on the conserved energy.
- KE_1 in the vet_sc arms with the predictor moves by about 1e-7. That is about 10x the
  ulp control, and it is the per-step 1e-8 stopping tolerance amplified by the convecting
  slab. HEAD's own confirming pass carries the answer to about 1e-12, not to 1e-8.
- **Verdict: PASS.** There is one reservation on vet_sc: KE_1 moves 1e-7 against 1e-8 for
  the controls.

## 4. GPU timing

Setup:
- apudev, node vipa1327, jobs 11944974 (E1), 11944975 (E2), 11944976 (V1), 11944977 (V2).
- The 3-D box is 84x104x104 with 4 MeshBlocks, pcr and bcg_sync = 1, with
  HSA_NO_SCRATCH_RECLAIM=1.
- Same binary (new4) for every arm, interleaved; repeat b runs in reverse order.
- ms/cycle is measured from cycle 20 to cycle 120.
- Configurations:
  - B = `implicit_lres_test=false implicit_conv_est=true`;
  - D = B + `implicit_predictor=step`;
  - D3 / D2 = D + `implicit_lin_ew_max=1e-3` / `1e-2`;
  - DC8 = D + `implicit_lin_cnorm=1e-8`.

| arm | E, 1 GPU: ms (a, b) | Picard | inner/step | V, 1 GPU: ms | Picard | inner | E, 2 GPU: ms | V, 2 GPU: ms |
|---|---|---|---|---|---|---|---|---|
| HEAD default (base binary) | 185.2 | 4.575 | 111.9 | 187.0 | 4.942 | 112.9 | - | - |
| ref (new4 default) | 184.3, 185.0 | 4.575 | 111.9 | 185.0, 188.1 | 4.942 | 112.9 | 153.8, 152.7 | 166.9, 166.1 |
| B | 182.7, 177.8 | 2.025 | 110.9 | 171.9, 171.2 | 2.025 | 110.4 | 147.3, 145.5 | 157.6, 158.8 |
| D | 130.8, 132.4 | 2.042 | 79.2 | 137.5, 138.2 | 2.042 | 82.5 | 106.7, 106.2 | 127.8, 126.3 |
| D3 | 123.7, 122.2 | 2.050 | 70.9 | 117.1, 120.6 | 2.558 | 63.6 | 99.2, 99.3 | 109.6, 107.8 |
| **D2** | **117.2, 115.6** | 2.858 | **62.9** | **116.3, 117.5** | 2.983 | **59.4** | **94.7, 93.8** | **106.6, 105.8** |
| DC8 | 137.1, 137.9 | 2.042 | 83.5 | 144.9, 145.0 | 2.042 | 86.0 | 111.4, 111.7 | 129.8, 131.5 |

- **Speed-up of D2 over ref:**
  - Eddington: 1.59x on 1 GPU (184.7 → 116.4 ms/cycle) and 1.63x on 2 GPUs
    (153.3 → 94.3).
  - vet_sc full: 1.60x on 1 GPU (186.6 → 116.9) and 1.57x on 2 GPUs (166.5 → 106.2).
- **NON-CONVERGED** is 0 in every arm.
- **Accuracy on the GPU** (at t = 19.4, D2 against ref): KE_1 4e-12 to 1.7e-8, totE ≤ 7e-12,
  F1top ≤ 2e-11.
- **Where the gain comes from:**
  - Cutting the passes alone (B) saves only 4 to 15 ms/cycle. The zero-iteration passes
    are cheap.
  - The predictor cuts the inner iterations by 29 %: pass 0 starts closer to the answer.
  - EW 1e-2 cuts another 20 %, by stopping pass 0 early. It costs about 0.8 extra cheap
    passes.

**Recommendation for GPU M1 inputs:**
`implicit_lres_test = false`, `implicit_conv_est = true`, `implicit_predictor = step`,
`implicit_lin_ew_max = 1.0e-2`. Keep `implicit_lin_cnorm = 0`. The defaults were left
unchanged so that existing inputs and bitwise gates do not move.
