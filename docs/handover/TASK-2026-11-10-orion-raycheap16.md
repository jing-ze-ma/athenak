# TASK 2026-11-10 viper -> ORION: raycheap16 -- is the cs p 0.58 plateau a non-smooth START, or a per-step term? (CPU)

Answer with `docs/handover/NOTE-2026-11-xx-orion-raycheap16.md`. Answers raycheap15 (f70711079): the cs transport-only
S0 plateaus at p 0.58 over 6 levels, uniform in every region; no switch moves it; the Cartesian control at the same
stiffness is 2nd order; S5 with a bitwise level-independent frozen tensor does not converge at all.

## viper code reading (fork/raycheap2-1109 dea3e3732, read only)
Per-step quantities on the cs S0 path that depend on the cycle count or on dt in a non-smooth way:

| # | quantity | file:line | on S0? |
|---|---|---|---|
| 1 | vet_gd direction-set rotation: cycle clock `VetGdBlock = cyc/N`, golden / R_3 quasi-random 3-D axis (`vet_gd_rotate_axis` default 3d on cs), K-window / EMA averaging | rad_m1_vetgd.cpp:1311-1356 (angle), 1451-1458 (block), 1470-1475 (last step); rad_m1.hpp:3423 `VgdRotOn` | **OFF.** The deck NAMES `vet_gd_rotate_every = 0`, so the fresh-run default `vet_gd_rotate_dt = auto` (20 x first dt) is not taken (rad_m1_vetgd.cpp:1181-1215: `rev_named` blocks it); `vet_rotate_avg = 1` -> vra_on false (rad_m1.cpp:2284-2291). Part 0 confirms from the logs. |
| 2 | vet_col / vet_col_lat build cadence on the cycle clock (`cyc % vcol_every`, `cyc % vlat_every`) | rad_m1_vetcol.cpp:685-688, rad_m1_vetlat.cpp:499-500 | inert: both default 1 (every step) |
| 3 | inner face flux `F0 (1 + 0.2 sin(2 pi t/40))`, set at the stage time t_n + dt for both hesdirk2 stages and the predicted build | rad_m1_time2.cpp:666-675, 2107-2114; read at rad_m1_implicit.cpp:3352 | ON, smooth; but the first (backward-Euler) step uses the t = 0 value; Orion's raycheap15 suggestion |
| 4 | **start-up discontinuity**: the IC is E = e_out = 1, F = 0 everywhere (rad_m1_tests2.cpp:926-945), the inner face jumps to F_in = 1 (f ~ 1, a beam) at t = 0 and the outer Marshak face to bath E = 0 (`implicit_ebath_x1max = 0`) at t = 0: two fronts at c, a non-smooth-in-time forcing; the step-1 BE | deck; rad_m1_implicit.cpp:3361-3368 (ramp key exists, unused) | ON |
| 5 | thin switch fse + hc correction take the stage OLD-vector E (iw EN -> TB) and the cell thinness max(|H|/J, eps C1/1.5) at the first pass only, so the restage (S6) re-reads the same old vector: S6 bitwise does NOT prove state independence; the stage operator carries E_old = U^n + a21 dt K1, not the stage iterate | rad_m1_latface_body.hpp:54-59, 118-127, 222-226; rad_m1_implicit.cpp:14504-14505 | ON; but thin_switch none (S4, C6) is also p ~0.4, so not alone |
| 6 | implicit_thin_freeze mask c dt kappa_P < 1e-2 (DEFAULT 1e-2 on a fresh run, a dt-dependent binary mask) | rad_m1_implicit.cpp:1835-1846, 14226-14233 | inert: opacity = const, coupling = false; and no cell crosses 1e-2 on the 4 finest levels (ambient c dt chi 0.074 ... 0.0093 only at the 6th) |
| 7 | Picard stop: `resid < tol` with implicit_conv_est (q/(1-q) estimate, no confirming pass; default true for fixed closures) | rad_m1_implicit.cpp:1895, 16619-16632 | ON, tol 1e-13 abs: 1e-9 of d -- expected inert; tested (B1) |
| 8 | lateral HLL / blend `ch dt/dx` inputs (implicit_blend_xthin), implicit_hr_excess (c dt/dx gate) | rad_m1_latface_body.hpp:106, 246; rad_m1_implicit.cpp:2393-2401 | inert: xthin 0, hr_excess none |
| 9 | the face-flux elimination theta_f = 1/(1 + c dt_sub kt_f) | rad_m1_implicit.cpp ~4913, 5008, 5042 | consistent (stage dt), smooth |
| 10 | one-pass auto, eccheck_every, k1_redist, posfix, pred_cap | rad_m1_implicit.cpp:2997, 17316; rad_m1_time2.cpp:850-917 | inert: off or counters 0 (raycheap15 Part 0) |

So the user's rotation / cycle-clock hypotheses are NOT active on this deck (1, 2). What is left and common to every
raycheap13-15 arm is item 4: the run starts from a state that is not a smooth solution (two fronts, a switched-on beam,
a BE first step) and is measured at t = 20, while the radial light-crossing is 5 and the diffusion time of the
ambient (chi 1, depth 5) is ~25: the start-up transient is still in the field at t_end. A non-smooth start lowers the
observed order of an L-stable 2nd-order scheme (the raycheap12 Marshak front: p 1.18 from t = 0, 1.58 from a common
restart at t = 10). The 2-D Cartesian control is a smooth Gaussian (no front), and S5 freezes the tensor of the t = 0
state (isotropic: f ~ 1 beam under D = 1/3) -- both fit. Part A tests this first; Parts B/C close the other items.

## Binary
NO REBUILD, no new key: every key below exists in dea3e3732. Use /orion/ptmp/jinma/raycheap14_1109/athena_dea3e3732
and conv_loc.py from fork/raycheap3-1110 (as raycheap15). If a key is refused or a run stops, quote the FATAL line and
continue.

## Base deck and ladder
S0 of raycheap15 (Ar4 cs deck + coupling = false, time2_vet_col = predict, dbg_hydro_skip_t2 = true; tlim 20; 24
ranks). **Finest 4 levels: cfl 0.0375 / 0.01875 / 0.009375 / 0.0046875** (the raycheap15 S0 walls 666 / 1037 / 1689 /
2765 s). Overrides on top of S0 (rad_m1/ prefix unless noted). Order of priority: A3, A1, A2, A0, A4, B1, B2, C2, C1.

## Part 0 (no runs, minutes)
- `athena -n -i <S0 deck c0.0375> <S0 overrides>`: quote the effective values of vet_gd_rotate_every, vet_gd_rotate_dt,
  vet_rotate_avg, vet_gd_rotate_axis, vet_col_every, vet_col_lat_every, implicit_thin_freeze, implicit_conv_est,
  implicit_tol_rel, implicit_one_pass, implicit_hr_excess, implicit_blend_xthin (absent = say absent).
- grep the S0 run.logs (c0.15 and c0.0046875) for `rotate`, `one_pass`, `thin_freeze`: expected no
  `vet_gd_rotate_dt = auto` line (rotation off).

## Part A: smooth start (best guess)
| arm | overrides on S0 | tests |
|---|---|---|
| A0 | time/tlim=10 and time/tlim=40 (two 4-level ladders) | p vs t_end: an initial layer gives p rising with t_end |
| A1 | implicit_flux_x1min_ramp_t0=4.0 implicit_flux_x1min_ramp_w=1.0 | inner beam switched on smoothly (factor 0.5(1+tanh((t-4)/1)), 3e-4 at t = 0; the modulation stays) |
| A2 | A1 + implicit_ebath_x1max=1.0 | the outer Marshak face starts in balance with E = 1 (no outer front) |
| A3 | common restart: ONE run of S0 at cfl 0.01875 with `output4/file_type=rst output4/dt=20.0` to t = 20; then all 4 levels `-r <that rst at t=20> time/tlim=40 time/cfl_number=<c>` (+ the output overrides of S0) | the start-up transient removed (as raycheap12 Marshak); every level starts from the same bits (first step BE) |
| A4 | implicit_flux_x1min_mod_amp=0 | steady inner face (raycheap15 suggestion) |
For A3 say the restart file used and check that each level ends at t = 40 with the same cycle count rule as S0.

## Part B: the remaining per-step items
| arm | overrides on S0 | tests |
|---|---|---|
| B1 | implicit_tol=1.0e-15 implicit_lin_tol=1.0e-15 implicit_conv_est=false | Picard / Krylov stopping (item 7); report Picard mean/max and NON-CONV (maxit 200) |
| B2 | closure=m1 (drop vet_* / thin-switch keys only as far as the FATALs demand: try first `closure=m1 vet_gd=false implicit_thin_switch=none implicit_flux_beam=closure implicit_thin_corr=none`, then strip what is still refused, and list exactly what was removed) | transport alone vs the vet closure on the cs deck |
| B2s | B2 + A3 restart recipe IF A3 passes and B2 runs (its own rst from B2 at c0.01875) | |

## Part C: the rotation (user hypothesis; OFF in S0 by item 1)
| arm | overrides on S0 | tests |
|---|---|---|
| C1 | vet_gd_rotate_every=8 | cycle clock (3-D axis): more rotations per unit time on finer dt -- what the hypothesis would look like |
| C2 | vet_gd_rotate_dt=2.0 (keep vet_gd_rotate_every=0 in the deck; align and slope take their fresh defaults) | time clock; quote the `vet_gd_rotate_dt` / align / slope echo lines |

## Report (NOTE)
Per arm: conv.py d and p (3 pairs of the 4 levels), conv_loc region p for every arm with p < 1.5, Richardson p_R,
Picard mean / max / non-converged per level, the hesdirk2 and m1-positivity end lines (changed counters only), wall per
level (s), cycles per level. Pass: **p >= 1.8 at the finest pairs**. If an A arm passes, run it also on the 2 coarse
levels (0.15, 0.075) so the 6-level ladder shows where it enters the asymptotic range. CPU only, no code changes.
