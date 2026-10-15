# NOTE: robust hesdirk2 predict (predfix-1011), 2026-10-11

Branch `predfix-1011` (local, not pushed), worktree /raven/u/jinma/ATHENAK/wt-predfix, based on prodcand-1011 e66ee669.
Fix commit: **eb75fa96**. CPU binaries: /raven/ptmp/jinma/predfix_1011/bin/athena_cpu_eb75fa96 and athena_he_cpu_eb75fa96.
Work dir: /raven/ptmp/jinma/predfix_1011 (run1.sh = suite_1011/run1.sh with the same arms, fixkeys and suite_1011/inp inputs; it only writes to this dir).

## Summary

The two literature failures have different causes. Only eq2 is caused by predict.

1. **eq2 (stiff uniform relaxation): predictor overshoot plus a one-sided periodic sweep.**
   - Mechanism (measured with the new `dbg_t2_pred`):
     - At the first predict build (cycle 1), P = U^n + dt K1 has E_P = 278.9 = 2 E^n.
     - eint_P <= 0 in all 1024 cells, so T_P = 0 and the SC source S = a T_P^4 = 0. The max relative change is r = 12.6.
     - The formal solution of P is then only the first-pass entry intensity.
     - With `vet_x1_periodic` and npass = 2, the upward rays enter with E + 3 F.Omega/c and the downward rays with 0. The SC therefore has H1/J = 0.77, uniform (0.74..0.80), in a uniform state.
   - Clip count: 0. The order-1 SC intensity is a positive combination of non-negative terms, so fmax(iv, 0) never acts. The symmetric-clipping hypothesis is not the mechanism.
   - Which term carries the flux:
     - `implicit_thin_corr = hcs` is OFF here: q *= esc, and esc = 1 - eps = 0 in a pure absorber.
     - The flux comes in through the thin switch (fse: g = max(|H|/J, eps C1/1.5); with S = 0, C1/J = 1, so g = 0.67 > h1 and w = 1). It is then carried by the half-range upwind face flux w c (h+_L E_L + h-_R E_R), with h+- = HP/J of the asymmetric SC.
     - A uniform E then gets the net flux (h+ - |h-|) c E. F1/cE goes 0 -> 0.77 in that one step.
     - After that, F^n feeds the upward entry intensity and keeps H1/J ~ 0.09. The gas is pushed to v = 1.71.
   - Lag and BE are also affected by the one-sided entry, but mildly: H1/J ~ 3 % and a non-uniform final E (143.35..149.45 around 145.79).
2. **Crossing beams ring 4: not caused by predict.**
   - On compute nodes `time2_vet_sc = lag` fails identically: same resid 2.067774e+25 at cycle 2, cell (0,0,13,146). This holds for 4c033dcf, e66ee669 and my build.
   - The suite's "lag passes" diagnostics (REC_diag_nopred/h2plain/c32) ran on the login node raven01. The outcome there depends on the machine and thread count: at OMP 36 the same case also passes.
   - The predicted state cannot matter in vacuum. S = a T^4 does not read E_P, and the SC of P is identical with and without the cap (mean J 3.125e-2, H1/J 4.34e-2).
   - Mechanism:
     - The run was killed by the **stage-1 solve itself**. Its Picard loop diverged (resid 1e28 -> 1e25 over 200 passes). It was correctly flagged NOT ADMISSIBLE, but the generic `implicit_resid_fatal` check ran afterwards and called exit, so the BE redo never ran.
     - A plain BE step from the same U^n converges in 4 passes (RECbe).
     - Dumps showed the redo's input state bitwise equal to RECbe's (u0, f0x, opac, vet_cell, hydro). The redo would have converged.

## Fix (eb75fa96)

- **`time2_pred_cap`** (rad_m1_time2.cpp, `Time2PredTheta`; both the vet_sc and the vet_col/vet_col_lat/vet_gd predict paths, which share `Time2VetColAt`):
  - P = U^n + theta dt K1 per cell, with theta = min(1, cap/r) and r = max(|dt K1_E|/E^n, |eint(P) - eint_b|/eint_b). eint_b is the stage-1 start gas without the coupling's K1.
  - One theta per cell is applied to E, the gas momentum and energy, and the vet_col inner face flux.
  - Default 0.5 on fresh runs; a restart whose file lacks the key gets 0 (off, bitwise). With cap 0.5, P stays positive.
  - Where r <= cap, theta = 1 and the arithmetic is bitwise the old predict. This was the case in every T-ladder and AG Car cell: 0 limited builds.
  - The report line counts limited cell-builds.
- **`vet_x1_periodic_inflow = source`** (default on fresh runs; `boundary` = old):
  - The first pass of the periodic x1 sweep enters both directions with the face source S. A uniform state then gives I = S on every ray and H = 0 exactly.
  - Measured: eq2 SC H1/J 0.77 -> 3e-16. eq3t final E spatial sd 2.6e-2 -> 4e-12.
- **No fatal on an already-failed stage**: `implicit_resid_fatal(_masked)` is skipped when `t2st && t2_fail`. The backward-Euler redo keeps the test.
- `dbg_t2_pred = N` (read only when named): prints the predicted state (min E_P, cells with E_P <= 0 or eint_P <= 0, max r, theta) and, for vet_sc, the SC-of-P mean/min/max H/J and the clip count.

## Validation (CPU, compute nodes)

| test | before (e66ee669 / suite 4c033dcf) | after (eb75fa96) |
| --- | --- | --- |
| eq2 REC (hesdirk2 + predict + tJe) | E 144.34654, T 3.46578, v 1.70575 | **E 145.78780, T 3.47480, v 0 (1e-21)**, uniform (min = max) |
| eq2 RECbe / lag | 145.78791 / 145.78790 (non-uniform 143.3..149.4) | lag 145.78781, uniform |
| eq2 cap only (inflow=boundary) | | 145.78790, v 1.5e-7 (= lag) |
| eq2 inflow=source only (cap 0) | | 145.78781, v 2e-21 |
| cb REC ring 4, c dt/dx 256 | FATAL cycle 2 | **completes, 20 cycles, bitwise = RECbe_r4** (16/20 stages still fall back to BE: positivity/offdiag fallback, as with lag) |
| cb REC r12, r4c0.5, cbg, cbh, shd r12/r36, shdg | pass | final outputs bitwise identical to the suite |
| cloud | xc 0.0276, max T 1.505, max rho 24.615 | identical to 4 digits (fields differ <= 5e-3; limiter active at fronts) |
| nltek1-8 | L1 / maxrel / E_top as suite | identical to 4 digits |
| sdifk | L1 0.0006/0.0028, peak 1.0001 | identical |
| ddiff, ddn, ddf | | bitwise identical |
| eq1 | E 96.79505 T 3.13663 | identical |
| eq3t (moving uniform gas) | E(200) 1.13329, sd 2.6e-2 | 1.13407, sd 4e-12 (uniform, as it must be) |
| G5 tJe_h2 k 12.8/128/1280/12800 | 0.999565 / 0.999874 / 1 / 1 | identical, PASS |
| marsh c1/c10 (vimp off, as suite) | | E max rel diff 5e-7 / 5e-6, F1 8e-6 / 8e-5 (68 limited cell-builds at the front, c1); F2 is round-off noise |

**T-order ladders** (od.py from thinfix_1011; inputs regenerated identical to thinfix/torder; levels 1-16, last orders):

| case | before | after |
| --- | --- | --- |
| pulse_k0.128 tJeP E | 1.94 1.97 2.00 | **2.01 2.01 2.01** (F1 2.01 2.01 2.01) |
| pulse_k12.8 tJeP / cenP E | 2.06 2.03 2.02 | identical (2.06 2.03 2.02) |
| pulse_k1280 tJeP E | 1.98 1.99 2.00 | identical |
| rw_t10a tJePvs E / F1 / dens | 2.67 2.00 1.91 / 1.96 1.92 1.87 / 1.98 1.96 2.08 | 2.67 2.00 1.91 / 1.97 1.84 2.01 / identical |
| rw_t1000a tJePvs (lev 2-16) E / F1 / dens | 2.05 2.09 / 2.00 1.99 / 1.99 1.99 | identical |

The limiter was never active in any ladder run (0 limited cell-builds), so the k12.8/k1280/rw_t1000a ladders are digit-for-digit identical. Only the periodic entry intensity changed: k0.128 (npass 16 at kappa L = 0.128, where the old one-sided entry was not negligible) improved to 2.01, and rw_t10a F1 moved slightly.

**vet_col stiff uniform case** (eq2c = eq2 with reflecting x1, zero-flux implicit BCs, closure = vet_col, hesdirk2 + time2_vet_col = predict, vimp):

| arm | E final | spatial range | v |
| --- | --- | --- | --- |
| be / lag | 145.787811 | uniform | 1e-16 |
| predict, cap 0 (old), relaxed or state source | 145.787796 | 145.7549..145.8251 | 3.5e-4 |
| predict, cap 0.5 (new) | 145.787811 | uniform | 1e-16 |

The vet_col path had the same overshoot (milder: its source does not go to 0). The cap limited all 1024 cells at the first predict build and gives the lag answer.

**Reduced AG Car A** (nx2 = nx3 = 32, 4 ranks, hesdirk2 + time2_vet_col = predict, 60 cycles):
- old e66ee669 vs new eb75fa96: hst bitwise identical, and the m1 / hydro_w / m1_thsw / m1_fs data are bitwise identical (the files differ only in the header's echoed input).
- 0 of 2.9e7 cell-builds limited, 0 stage fallbacks.

## Open / notes

- In vacuum, cb stage-1 solves still fail admissibility on 16/20 steps, so REC is BE there. The cause is the offdiag positivity fallback (od_now != impl_offdiag), the same with lag. It is now harmless but not fixed.
- The default of `vet_x1_periodic_inflow = source` changes every `vet_x1_periodic` result (lag/BE included) at the level of the old entry error, e^{-npass tau}. That is visible for npass = 2 (eq tests) and for k0.128. The vet_x1_periodic path is single-MeshBlock test-only and is not used in production.
- The suite agent's "lag passes cb" result is a login-node artifact. On compute nodes, prodcand REC lag fails cb ring 4 identically to predict.
