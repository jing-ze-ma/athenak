---
name: he-presn-m1-wedge-0929
description: He presupernova model restarted on implicit M1 (user 09-29): sp wedge, 4.0 Msun, new pgen he_star_m1.cpp, 2nd order + fully coupled
metadata:
  type: project
---
User decisions 09-29: 4.0 Msun presn He star (Woosley 2019) on a spherical-polar WEDGE (theta pi/4..3pi/4 PERIODIC (symmetric about the equator), phi pi/2 periodic,
r 0.5 R .. tau 1e-2), NEW lean pgen src/pgen/he_star_m1.cpp (red_giant IC/gravity + box_convection M1 hookup; no
two-stream, no MLT closure), branch he-presn-m1. Must be 2nd order in time AND space with fully coupled gas-radiation:
the box 2nd-order key set (hesdirk2, force_reference_work split, implicit_opac_update, one_pass 0, tol 1e-8) + vet_col.
Plan: 1a column gate (force balance, shell L, energy ledger) -> 1b dt-order gate on the wedge (order >= 1.9 in rho, E,
v_r, v_h, T; cfl 0.075..0.9) -> 1c 3-D scout. Unverified on sp: split with point mass, vet_col tensor stage time
(predict vs extrapolate), implicit_opac_update, plm enthalpy on stretched r. Supersedes the two-stream
[[he4-presn-global-plan]]. Data dir /viper/ptmp2/jinma/hepresn_0929.

**Step 1a result 09-29 ~05:30** (branch he-presn-m1, worktree /viper/ptmp2/jinma/wt_hepresn, 6 commits, unpushed; data/runs
/viper/ptmp2/jinma/hepresn_0929): pgen + input + IC script built; nx1 640 uniform, r 0.5-1.0286 R, dfloor 1e-14 (table
min), dt 13.45 s set by the inner cells. Hopf (radiative-equilibrium) IC gate: force <= 2e-3 rho g below 0.9 R but
4e-3..1.3e-2 at tau 1-84 (1-2 cells per gas Hp there; dt-independent, survives damping); wall cell 2e-3 (1st order in dx);
L within 1.7e-3; ledger -2.1e-4 = one-time startup offset (verified); restart bitwise; gas-Newton fallbacks 9e-5
(location unknown); hopf unstable near the inner wall after t > 1000 s. Woosley/MLT IC (col) heats at ~1e-3/s in the FeCZ
(no convection yet), force up to 2.7e-2. OPEN: IC choice for 3-D (col vs hopf), stretched radial grid.
NAMING: the "Woosley column"/col IC is OUR 09-17 structure (bench/hestar_presn/column_sph.py, /viper/u2/jinma/ATHENAK/bench/
hestar_presn; Woosley gives only M, L, g): ideal mu 1.342 + aT^4/3 EOS, grey Eddington, KW MLT alpha 1.5. The new IC mode
mlt (b4290563) redoes it with the run's gas-only table EOS, M1 force form, Hopf atmosphere, discrete hydrostatics.

**Opus agent final report 09-29 ~08:30** (head d86c5695, binary builds/bin/athena_he_cpu_d86c5695, runs /viper/ptmp2/jinma/hepresn_0929/v2):
esrc core source (16ca48e9, bitwise off; hesdirk2 energy exact 1e-4/step; be loses 0.2 %/step systematically, unexplained);
he_ic_balance exact discrete WB; wall cell fixed (module face form, b33af5e6); frozen-MLT column balanced < 6e-4 rho g
except the photosphere (9e-3 at tau ~1 = vet_col radiates 0.2-1 % more than Hopf; relaxation 8.4 -> 5.4e-3 over 3000 s;
re-balancing per restart diverges). Grid P160 (nx1 160, grid_P160c.npy): dt 32.6 s vs 11.7 s on 640 (~11x cheaper);
FeCZ-edge force 1.4-2.7e-3 (> gate; refine edges). BLOCKER: inner-wall instability r/R 0.516-0.525 is NUMERICAL,
dt-dependent: none at cfl 0.15, onset 1485 s at 0.3, 595 s at 0.6, be 0.3 at 373 s; E drops first, T_rad < T_gas
although t_eq 3e-4 s (implicit exchange not holding LTE). Suspect rad_m1_implicit.cpp:6600-6717 stage old vector.
Blocks cfl >= 0.3 and the 3-D scout. Core fix awaits the user's go.
**WALL INSTABILITY FIXED 09-29 ~09:40 (5a908941 on he-presn-m1, + fb7a1270 merge of d26b7364):** root cause = vet_col
formal-solution source from the start-of-step GAS temperature (rad_m1_vetcol.cpp:748/1043), which carries the gas-only
compression transient that the stiff exchange removes (P_gas/P 0.01 -> 1/800 survives); per-column f_K differences
-> transverse checkerboard of v1 (column means hid it; "LTE broken" was a misreading). Fix: vet_col_source = relaxed
(default; local energy-conserving BE exchange T*,E*), gas = old bitwise. Gates: hesdirk2 cfl 0.3/0.6/0.9 + be 0.3 stable
to 3000 s, Newton on and off; column numbers unchanged; He box bitwise (vet_sc untouched). OPEN: wall cell -60 cm/s^2
(boundary WB, v -1e4 cm/s by 300 s); be esrc loss 1.9e-3 (separate); same fix possible for vet_sc later.

**09-29 ~18:45 USER GO: ADAPTIVE SHELL-MEAN DEFICIT CLOSURE** replaces the fixed ramp as the start design: F_sub -> max(0, F_req - <F_rad> - <F_conv,res>) relaxed over mlt_relax_time (start 0.2-0.5 turnover), via esrc, FeCZ only, restart state in PGENST01. Reason: frozen MLT w=1 run (M1/w1, job 12025960) double counted (L_top/L 1.17, mass loss ~2e20 g/s, dt 97 -> 3 s) once resolved convection grew. Lessons: [[red-giant-mlt-closure-lessons]], [[he4-porosity-luminosity-excess]].

**09-29 ~19:30 3-D RESULT: POROSITY, closure not at fault.** w1 (frozen) and ad3d (adaptive, relax 1500 s; 7b1e49b6) identical to ~3.5 turnovers; convective growth e-fold 2520 s (0.54 turnover), v_r/v_MLT 1.06 at 5 turnovers but L_conv/L <= 0.006 (very inefficient); from 3.5-4 turnovers rho rms/<rho> 0.26 -> 1.40, <F>/F_diff(mean) 1.02 -> 2.4 (3.3 at 0.75 R), corr(F,rho) -0.6..-0.8, local diffusion matches M1 flux 1-3 % (no M1 artefact); L_top/L 1.16, mass loss 3-4.5 %/turnover, dt 45 -> 3 s. Same as the 09-20 two-stream He4 porosity result. Open: steady porous state vs runaway (run longer, cheap ~10 min/turnover at dt 3 s), resolution convergence of rho, thin-column inner-wall instability at cfl 0.9 (1-D gate past 0.55 turnover). NHISTORY/NREDUCTION 20->22 global change to flag before merge. My earlier "w_mlt garbage" was my own column misread.

**09-29 ~21:30 DIAGNOSIS (diag_rst A/B/C):** restart faithful (B exact), ramp innocent (A no-ramp NaN at 5.14 tt), porous runaway cfl/BE-INDEPENDENT (C cfl 0.3 = orig to 0.5 % through 4.54 tt) -> physics up to 5 tt: runaway, not steady (P 1 -> 2.8, L_top 1.16, mass loss 4.5 %/tt). Implicit M1 then fails (NaN ~5.1 tt); failures concentrate on the MPI RANK faces (phi faces 52-96 % vs 6 % even) -> suspect cross-rank implicit exchange / lagged transverse terms. Proposed: 1 vs 2 ranks from rst 00016 (not run, awaiting user). r2 128x128 queued on apu (~05:17). Commits aef58241 06708947 603e9373 on he-presn-m1.

**RANK TEST 09-29 ~21:30 (diag_rank n1/n2): NO rank-boundary defect.** 1 rank vs 2 ranks round-off close before failure, first stage failure at the same cycle 355 (positivity), same counts and trajectory; rank face only pulls the worst-residual location. The earlier "phi-face 0.6-0.96" Newton-fallback location was a print-order artefact (first 16 cells per rank in m,k,j,i order). Stage failures come from the state (FeCZ top 0.93-0.95 R once drho ~0.7-0.8). Porous runaway to 4.5 tt robust (cfl, ranks, BE); NaN ~5.1 tt ends it. Next candidate (code change, not run): which positivity term fires first (rad_m1_implicit.cpp:8557-8568).

**POSITIVITY DIAGNOSIS 09-29 ~22:40 (dbg_t2_admiss, ec269634/a1c6833a; M1/diag_adm):** the stage SOLVE gives E < 0 (-1e3..-6e5 vs 3e6) in fast (v_r 3e7) underdense plume cells at the FeCZ top; E old vector fine; gas-eint old vector < 0 (risk 6.1) secondary; Picard reports converged with the iterate clipped at e_floor; BE redo SILENTLY floors E (adds energy). ALL 13 minimum cells at j = 3 = the theta-periodic seam (only minima printed). My proposal to the user: dump all failing cells + theta-reflecting control first, then aphll vs central; count BE floor clips.

**09-29 ~23:30 THETA-PERIODIC WEDGE IS AN ARTEFACT SOURCE (diag_seam P/R):** all 44 failing cells at global theta index 0/1/63 (the seam); phi-averaged FeCZ flow develops an equator-symmetric converging meridional circulation feeding the seam (v_r +24 km/s at the seam vs -1.5 at the equator, v_theta +-16 km/s toward the seam at t 18801, from ~2.5 tt), low-rho cells peak at the seam; seed phases random -> not the seed. cot(theta) flips sign across the seam: theta-periodicity is not a symmetry of the sphere. Reflect restart from the contaminated state breaks down worse. => the porous-runaway numbers of ad3d/w1 are CONTAMINATED. My 09-29 statement that periodic theta is fine (area match) was wrong. BE floor-clip counters 081bba66 (he-presn-m1). Next: t=0 run with reflecting theta walls (or other seam-free geometry); RUN 2 (same geometry) recommended cancel.

**RETRACTION 09-29 ~23:10 (wr3d, theta reflecting walls from t=0):** the converging meridional circulation is NOT a seam artefact: with walls the same equator-symmetric box-scale convective mode grows (upflow near the theta edges, peak just inside the wall; sinking at the equator; v_r ~90 km/s at 4.5 tt), macroscopic L_top, M, drho, P, dt, hk4 identical to ad3d within ~1-5 %, first stage failure 19338 vs 19300, worst cells again along the plume at the theta edges. The cot(theta) mechanism I told the user was WRONG. ad3d physics is not seam-contaminated; porous runaway + solver breakdown belong to the box-scale upflow plume at the FeCZ top. The 128x128 run was cancelled on a wrong premise. Next candidate: positivity-preserving implicit radiation advection in fast plumes (aphll/blend where Peclet > 2).

**ROOT CAUSE 09-29 ~23:45 (diag_seam D/N, verified run logs):** the opacity-NEWTON face term (M1OpnCell, implicit_opac_newton, DEFAULT ON since the m1-perf merge d26b7364) makes the sp-row diagonal TB negative in fast thin plume cells (term -1e4 vs diffusion +800): s = q kt/bk, q = -0.5 chat dt th g_f, sign follows face flux and d kappa/dT (Fe bump), magnitude unbounded (rad_m1_implicit.hpp:435-441, .cpp sp rows ~7911/~8010). implicit_opac_newton = false: 0 floor clips (463), 0 NC (5), Picard mean 10.3 (19.4), pre-failure solution identical to tolerance (1e-8). Proposed fix: per-face guard keeping the diagonal >= 0.5 x its no-Newton value (skip diag + rr part together), counted; fixed point unchanged. Remaining milder trigger: vimp positivity fallback (16 stages). AFFECTS every M1 run with the Newton default (He box too?).

USER 09-29 ~23:55: GO guard fix (agent, branch m1-opn-guard, merge after combined gate, NOTE for DeltaAI); He presn wedge input implicit_opac_newton = false until then (committed on he-presn-m1); 128x128 run resubmitted WITH implicit_opac_newton=false (12028935-37); vimp positivity fallback diagnosis LATER.

USER 09-30 ~00:55: He presn wedge input back to GUARDED Newton (ae9e5c77); 128x128 run resubmitted with binary athena_he_gpu72_0b8b6c0d (guard) and implicit_opac_newton=true: jobs 12030056 12030057 12030058.

USER 09-30 ~09:15: He presn wedge input rsolver -> hllc (abf547b5 on he-presn-m1) and M1/he3d_M1.athinput (backup .pre_hllc); 128x128 run RESUBMITTED as 12037023 -> 24 -> 25 (old 12030056-58 cancelled). Not directly comparable with the lhllc 64x64 ad3d/wr3d.

USER 09-30 ~09:20: 64x64 hllc run from t=0 (M1/ad3d_hllc, r3.sh, jobs 12037044 + 12037045; same keys/binary as the 128 run: guarded Newton, adaptive closure, ramp 23500/4700, tlim 37600) for a clean hllc resolution pair with the 128 run.
