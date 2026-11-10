# NOTE 2026-11-10 ORION -> viper: raycheap15 -- bisect of the cs transport-only 1st order

Answer to `TASK-2026-11-10-orion-raycheap15.md`. No code changes. Binary /orion/ptmp/jinma/raycheap14_1109/athena_dea3e3732
(`git diff dea3e3732 930aad615 -- src CMakeLists.txt` is empty, so no gate). conv_loc.py from fork/raycheap3-1110.
Job 228167 (6 p.exclusive nodes, 24 ranks per cs run, 1 rank per 2-D run, all arms concurrent), 1 h 05 m. S0 levels
c 0.15 ... 0.01875 are the raycheap14 S0 runs (same binary, bitwise = raycheap13 S0); the 2 extra S0 levels are new.
Dir /orion/ptmp/jinma/raycheap15_1110 (gen15.py, ana15.py, ana_final.txt, convloc_S0_6.txt, part0_hesdirk2.txt).

**Summary**
- **Not pre-asymptotic: S0 over 6 levels plateaus at p 0.58.** p 0.16 / 0.46 / 0.58 / 0.58 (finest 4: 0.58 0.58);
  Richardson p_R (3 finest) 0.58. Picard mean 13.96 -> 3.10, never non-converged. The difference is spread uniformly:
  every conv_loc region has the same p as the whole (0.58 / 0.58 at the finest pairs), and every MeshBlock holds 1/24.
- **No arm reaches p >= 1.5 (none reaches 0.7), so no extra levels were run.** Best: B2 (time2_vet_col=rebuild) and C5
  (implicit_flux_x1min 0.3), both 0.59 at the finest pair (from 0.43 / 0.29), i.e. the same plateau.
- **The stiff transition and stage order are ruled out.** A1 (no clump) is p 0.16 / 0.46 = S0 (solution 5% different,
  error the same); A2 / A3 (trbdf2) 0.38 / 0.46. The Cartesian control K1 (kappa_s 128, c dt chi 31 -> 1.9, crossing the
  transition) with the D3 keys is **p 1.83 / 1.83 / 1.87**, K1T (trbdf2) 1.83 / 1.83 / 1.86: hesdirk2 stays 2nd order
  through c dt chi ~ 1 in Cartesian geometry. So the culprit is on the cs deck, not the hesdirk2 stage order.
- **Inert on this deck (bitwise = S0 in m1_e):** B3 (pred_cap 0; differs 7e-8 at c0.15 only, where 2424 of 3.0e7 builds
  were capped), C2 (cs_vertex_3cell false), C3 (k1_redist 0, posfix 0), C4 (pos_floor / floor_solve / gd_allow_clamp
  off). Part 0 agrees: every positivity / fix counter is 0 at every level.
- **Act on the solution, not on the order:** C1 (cs_cross_mono / cs_vertex_mono off) 1.1% different, p 0.16 / 0.45;
  C6 (thin switch none + picard + pass + trbdf2) 0.37 / 0.41 with d 3x S0; B1 (gd rebuild every pass) -0.13 / 0.47 at
  3.7-4.1x wall (finest 3885 s, all 4 levels run).
- **B4 (S5 + vet_source_noesrc) still does not converge** (p 0.18 / 0.10, d 1.07e-3 -> 8.8e-4) although S5's frozen
  tensor is BITWISE level-independent (m1_lat, all 12 variables, at t=0 and t_end). The frozen build is dt-free; the
  non-convergence comes from something per-step outside the tensor.
- **K2 (kappa_s 12.8) with the D3 keys FATALs at cycle 0 on all 5 levels**: `Picard NON-CONVERGED after 200 passes:
  resid=4.86 lin_resid=0.106`, then `Picard solve DIVERGED ... > implicit_resid_fatal=1e-2` (rad_m1_implicit.cpp:17211).
  Extra (not asked): the same opacities on the D0 deck (default step / step lags) run: K1d p 1.04 / 0.79 / 0.88, K1dT
  0.92 / 0.85 / 0.92, K2d 0.40 / 0.40 / 0.85 (d 0.49 -> 0.16, pre-asymptotic), K2dT 0.35 / 0.47 / 0.69.
- **conv_loc.py region split is mis-mapped for this layout:** it warns `logical x2 centres [-17, 6] x3 [-1, 23] (the
  region split assumes a panel spans [-1, 1])`, and "vertex" gets 221184 of 442368 cells, "seam" 156672. The region
  numbers below are as printed; the per-MeshBlock 1/24 shares and the identical region p are layout-independent.

## Part 0 (existing raycheap13/14 S0 and S5)
hesdirk2 end lines (`<rad_m1> time_scheme=hesdirk2: ...`; there is no line starting `<rad_m1> hesdirk2:`), per level
c 0.15 / 0.075 / 0.0375 / 0.01875:

| run | stage steps | BE | stage fb | vet clips | k1_redist cells / unpaid | predict fb | pred_cap limited builds | posfix steps / re-solves / cells | accepted-with-fb |
|---|---|---|---|---|---|---|---|---|---|
| S0 | 67 / 135 / 270 / 541 | 1 | 0 | 0 | 0 / 0 | 0 | 2424 of 2.96e7 / 0 of 5.97e7 / 0 of 1.19e8 / 0 of 2.39e8 | 0 / 0 / 0 | 0 |
| S5 | 67 / 135 / 270 / 541 | 1 | 0 | 0 | 0 / 0 | (lag: no predict line) | - | 0 / 0 / 0 | 0 |

`<rad_m1> m1-positivity:` identical in all 8 runs: clipped cell-passes 0, implicit_pos_gas cell-solves 0, E floor raises
0, energy created 0, |F| > cE scaled back cell-solves 0, sum(|F|/(cE)-1) 0.
Full lines: part0_hesdirk2.txt in the dir.

S5 tensor / face dumps vs c0.15: m1_lat (m1_tt0, m1_lat0-5, m1_rx/ry/rz, m1_jg, m1_hg) bitwise at t=0 and t=20 for all
levels; m1_face (fx1, fx1o, fx2, fx3) bitwise at t=0, at t=20 max rel 5.6e-2 / 1.4e-1 / 2.7e-1 (c 0.075 / 0.0375 / 0.01875).

## Part 1: S0, 6 levels (c dt chi clump 14.8 ... 0.46)

| level c | cycles | Picard mean/max/nc | wall (s) |
|---|---|---|---|
| 0.15 | 68 | 13.96/24/0 | 294 |
| 0.075 | 136 | 10.76/19/0 | 435 |
| 0.0375 | 271 | 7.93/16/0 | 666 |
| 0.01875 | 542 | 6.04/12/0 | 1037 |
| 0.009375 | 1083 | 4.14/8/0 | 1689 |
| 0.0046875 | 2165 | 3.10/6/0 | 2765 |

conv.py (6): d 8.985e-04 8.034e-04 5.852e-04 3.916e-04 2.625e-04, **p 0.16 0.46 0.58 0.58**; finest 4: p 0.58 0.58.
conv_loc (6), p per pair rest / seam / vertex / inner / outer / clumppatch:
- 0/1: 0.16 0.16 0.16 0.15 0.17 0.14
- 1/2: 0.46 0.46 0.46 0.44 0.46 0.46
- 2/3: 0.58 0.58 0.58 0.57 0.58 0.55
- 3/4: 0.58 0.58 0.58 0.56 0.58 0.56
Top MBs every pair: 0.04 each (= 1/24). Richardson p_R 0.58; err vs E_R 3.48e-3 2.58e-3 1.77e-3 1.19e-3 7.97e-4 5.34e-4.

## Part 2: arms (4 levels, c 0.15 ... 0.01875; all rc 0, NOT ADMISSIBLE 0, DT_RETRY 0, Picard never non-converged)

| arm | d (3 pairs) | p | vs S0 m1_e (c0.15 / c0.01875) | Picard mean per level | wall per level (s) |
|---|---|---|---|---|---|
| S0 | 8.979e-04 8.028e-04 5.848e-04 | 0.16 0.46 | - | 13.96 10.76 7.93 6.04 | 294 435 666 1037 |
| A1 no clump | 8.986e-04 8.034e-04 5.853e-04 | 0.16 0.46 | 5.1e-2 / 5.1e-2 | 13.96 10.68 7.90 6.03 | 302 456 685 1122 |
| A2 trbdf2 | 8.685e-04 6.675e-04 4.850e-04 | 0.38 0.46 | | 13.60 10.14 7.31 5.42 | 294 429 632 1006 |
| A3 A1+A2 | 8.686e-04 6.677e-04 4.852e-04 | 0.38 0.46 | | 13.56 10.05 7.26 5.41 | 298 419 628 1001 |
| B1 gd_rebuild_every 1 | 7.306e-04 7.972e-04 5.775e-04 | -0.13 0.47 | | 16.05 11.06 8.74 7.00 | 1198 1653 2529 3885 |
| B2 vet_col rebuild | 1.273e-03 9.462e-04 6.278e-04 | 0.43 0.59 | | 13.33 9.95 7.18 5.28 | 312 499 745 1245 |
| B3 pred_cap 0 | 8.979e-04 8.028e-04 5.848e-04 | 0.16 0.46 | 7.4e-8 / bitwise | = S0 | 280 435 664 1070 |
| B4 lag + every 1e6 + noesrc | 1.067e-03 9.435e-04 8.775e-04 | 0.18 0.10 | | 12.44 9.34 6.54 4.73 | 230 322 461 667 |
| C1 cross/vertex mono off | 8.979e-04 8.028e-04 5.866e-04 | 0.16 0.45 | 1.1e-2 / 1.1e-2 | 19.66 14.24 9.76 6.99 | 365 555 794 1218 |
| C2 cs_vertex_3cell false | = S0 | 0.16 0.46 | bitwise / bitwise | = S0 | 289 443 676 1081 |
| C3 k1_redist 0 posfix 0 | = S0 | 0.16 0.46 | bitwise / bitwise | = S0 | 290 456 668 1077 |
| C4 floors / clamp off | = S0 | 0.16 0.46 | bitwise / bitwise | = S0 | 288 437 690 1083 |
| C5 flux_x1min 0.3 | 1.732e-03 1.421e-03 9.411e-04 | 0.29 0.59 | | 13.49 10.19 7.26 5.39 | 273 407 618 990 |
| C6 thsw none + picard + pass + trbdf2 | 2.936e-03 2.276e-03 1.717e-03 | 0.37 0.41 | | 13.61 10.14 7.31 5.43 | 303 446 663 1055 |

Key checks in the logs: A2 `time2_tableau = trbdf2: TR-BDF2 implicit stages (stage order 2)`; B2 `time2_vet_col=rebuild`;
B3 the pred_cap counter is gone from the hesdirk2 line; C3 the k1_redist and posfix counters are gone; C2 the default
echo `<mesh>/cs_vertex_3cell = true` is gone; C4 `implicit_pos_floor = false`. The hesdirk2 / m1-positivity end lines of
every arm are as in Part 0 (all fix counters 0) except as noted.

## Part 3: Cartesian 2-D control (raycheap13 D3 deck, 5 levels, 1 rank)
kappa_s lives in `<rad_m1>` only (`kappa_s = 1280.0` in the D3 deck; the rad_m1_beam thick_pulse pgen reads no problem
opacity). Edited there.

| arm | d (4 pairs) | p | Picard mean/max per level | wall (s) |
|---|---|---|---|---|
| K1 kappa_s 128 + D3 keys | 2.760e-03 7.752e-04 2.173e-04 5.955e-05 | **1.83 1.83 1.87** | 31.5/57 29.3/41 27.2/33 24.0/37 21.3/38 | 29 39 63 98 160 |
| K1T K1 + trbdf2 | 2.774e-03 7.822e-04 2.207e-04 6.092e-05 | **1.83 1.83 1.86** | 31.4/57 29.3/41 27.1/33 24.0/37 21.3/38 | 30 41 62 98 160 |
| K2 kappa_s 12.8 + D3 keys | FATAL at cycle 0, all levels (Picard diverged, resid 4.86) | - | - | - |
| (extra) K1d kappa_s 128, D0 deck | 1.192e-02 5.795e-03 3.355e-03 1.829e-03 | 1.04 0.79 0.88 | 3.00/3 | 4 6 9 13 24 |
| (extra) K1dT + trbdf2 | 1.961e-02 1.040e-02 5.762e-03 3.051e-03 | 0.92 0.85 0.92 | 3.00/3 | 5 5 9 13 22 |
| (extra) K2d kappa_s 12.8, D0 deck | 4.868e-01 3.678e-01 2.793e-01 1.550e-01 | 0.40 0.40 0.85 | 3.00/3 | 6 7 11 16 28 |
| (extra) K2dT + trbdf2 | 4.984e-01 3.905e-01 2.815e-01 1.746e-01 | 0.35 0.47 0.69 | 3.00/3 | 7 8 11 16 27 |

## Not tested (suggestion)
Every cs arm sits on the same p ~ 0.5-0.6 plateau with the error spread uniformly in angle. One time-dependent input
not in the TASK's list is the modulated inner flux (`implicit_flux_x1min_mod_amp = 0.2`, `_mod_period = 40`): an arm
with `implicit_flux_x1min_mod_amp = 0` (steady inner face) would show whether the stage-time evaluation of that boundary
term is the 1st-order piece. Not run here (no TASK); cheap (4 levels ~20 min) if you want it.
