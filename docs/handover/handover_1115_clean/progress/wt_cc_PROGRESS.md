# STATE AT HANDOVER (11-15, viper -> Orion)

Branch clean-cc, head dc67dba7b.  Worktree /viper/ptmp2/jinma/clean_1114/wt_cc.
Last binary: /viper/ptmp2/jinma/builds/bin/athena_none_cpu_dc67dba7_cc19 (none_cpu).
Debug-only branch dbg-bstress (worktree /viper/ptmp2/jinma/clean_1114/wt_dbgB, local, NOT for merge): CLHIST/CLMAG
(realizability maps), CLBUD (stage budgets), CLMOM (total momentum), clean_dbg_itstage/closn/pform/bstress.

## Keys (<rad_m1>, scheme = clean)
- clean_layout = staggered (default) | cell.
- clean_closure = predict (staggered default; rays once per step on the predicted state) | iterate (rays of the
  iterate every Picard pass) | extrap (CELL DEFAULT: rays once per step on the converged U^{n+1}, moments extrapolated
  linearly in time; stage 1 rebuilds the closure once from the rays of its first iterate; K1 re-evaluated with the
  end-of-step closure, cell) | extrap2 (also the stage-1 moments) | eddington (D = I/3, h+- = +-1/4, k+ = I/6, C1 = 0,
  no formal solution).
- clean_alpha = tau (default) | tau_sensor | tau_ray | one | zero; clean_recon = plm (default) | none; clean_psi_eps 1e-12.
- clean_maxpic 12 (Picard passes for lagged slopes / closure), clean_maxnl 50, clean_tol 1e-10, clean_tol_nl 1e-11
  (cell: effective max(tol_nl, 10 tol)), clean_maxit 4000, clean_pc_sweeps 1 (multigrid V-cycle preconditioner).
- clean_pos_limit true (a-posteriori low-order faces, clean_maxpos 3; BE ends with every interior face low-order;
  in that all-low-order solve the boundary stress is the upwind flux against the ghost, else the physical half-range
  stress c k_out E + incoming beam/bath stress), clean_pos_boundary false, clean_diag_cells false,
  clean_ray_retard false, clean_ray_min 1e-3, clean_fill_passes 8, clean_w_c0/c1 0.1/0.5, clean_vfreeze 0 (diag).
- A BE solve that does not converge fails the cycle (dt_retry) with a CLEAN WARNING.

## Passed / failed (numbers below in this file)
- rw time order (P 3.3 / 1/3, tau 10, N256, amp 1e-2): iterate 2.00/2.00; extrap + 1 stage-1 closure pass 2.01-2.00 /
  2.00-1.99 at 2.05 formal solutions per step (debug build 40df0f05 = same logic as dc67dba7b; cc19 itself NOT yet
  run).  Plain extrap 0.90 / 1.65; production-form predictor 1.00 / 0.89.  prod-1118 and unify-1114 (vet_sc, predict):
  2.03 / 1.99.
- rw space (amp 1e-6): cell 1.9-2.1 (ray closure, cc4..cc11); fixed D (eddington) cell 1.87-2.37 every case;
  staggered eddington rows INVALID (Newton stagnation at tol_nl 1e-13 -> failed stages -> BE / dt_retry).
- nltek1 x1 vs exact 9.4e-5 (staggered 9.5e-5); hcyl x0.5 0.0468 (staggered 0.025); shd/shdg x0.25 0.264/0.256
  (staggered ~0.02); cb/cbg/ws66 x0.5 score 1.68/1.62/1.35 (staggered L23 1.47/1.40/1.26).
- Momentum (rw periodic, all cells, CLMOM): conserved to round-off, 5.5e-17 (tol 1e-14) / 1.3e-16 (tol 1e-10) of
  sum|rho v|.  The reported 2e-14..7e-11 drift is an artefact of rw/conserve.py (tab outputs are slice_x2, one of the
  4 rows; inter-row momentum exchange through x2 faces is missed).

## Open problems (priority)
1. Momentum conservation: the code conserves to round-off (above); no non-conservative term was found.  FIX THE
   MEASUREMENT (conserve.py must sum all rows: drop slice_x2 or use bin outputs) and re-verify on Orion.
2. 2-D beam realizability: every high-order stage solve on cb/cbg/ws66 is non-realizable (~3100-3300 cells; shd/shdg
   ~6500), so every step falls back to BE; after the all-low-order BE 238/538/320 |F|>cE cells remain (cb/cbg/ws66
   x0.5, first step): the first-order D = I implicit scheme itself is not positive in these 2-D transients.
3. Accuracy shd 0.26 / hcyl 0.047 vs staggered 0.02 / 0.025.
4. Cheap closure: stage 1 needs the closure of its own state Y1 (one rebuild from the first iterate suffices, 2 solves
   per step).  A 1-solve predictor (production's predict) is order 2 in prod but first order in clean; not understood.
5. Fixed-D order matrix (both layouts, space 1e-6, time 1e-4/1e-5, 1e-6 noise and fine-grid checks) not run;
   staggered needs clean_tol_nl >= 1e-12.  cc19 run lists: /viper/ptmp2/jinma/clean_1114/farm_req/.

## Where things live
- rw driver /viper/ptmp2/jinma/clean_1114/rw/run_radwave.py (caselists rw/*.txt; refs radwave_disp.py: vetqs_d for the
  ray closure, edd for D = I/3), ladcmp.py (successive differences), vsit.py (arm vs reference arm), edtab.py (space
  orders), conserve.py (FLAWED for nx2 = 4); farm front end /viper/ptmp2/jinma/clean_1114/farm/rwfarm.py.
- IDORT runner suite/run_ladder.sh, scorers cc/score_any.py, cc/nlte_score.py (exact transport ana/nlte_exact.py),
  decks /viper/ptmp2/jinma/blend_1022/p2/decks; beam1d/45/45n decks cc/ret/I_beam*/in.athinput.
- Run dirs cc/rw_* (ladders), cc/ret/* (beam diagnostics), cc/momdbg (CLMOM), suite/<deck>/<ARM>_x<f>.

# clean-cc progress (cell-centred clean layout), 11-15 PST

Worktree /viper/ptmp2/jinma/clean_1114/wt_cc, branch clean-cc off clean-1114 c20520ad5 .

## Code
- `<rad_m1>/clean_layout = staggered | cell` (default staggered = old code path).
- cell layout keys: `clean_alpha = tau (default) | tau_sensor | tau_ray | one | zero`,
  `clean_recon = plm (default) | none`, `clean_psi_eps` (sensor floor, 1e-12).
- Files: src/rad_m1/rad_m1_clean_cell.hpp (device functions: CleanDissipationWeight, jump
  coefficients, face Jacobian/flux, block solve), rad_m1_clean_cell.cpp (host steps);
  CleanSolve / CleanBiCGStab generalised in rad_m1_clean.cpp / rad_m1_clean_solver.cpp.
- Unknowns E, y = F/c; PLM (van Leer) jumps with lagged weights, Picard passes folded into
  the Newton loop (tol = max(clean_tol_nl, 10 clean_tol)); x1-line block-tridiagonal
  (4x4) preconditioner; boundary: E flux = the staggered boundary forms, F flux = zero
  gradient stress (reflect: normal stress only).

## Commits (clean-cc, local)
- 5b14afdcc cell layout (cc1 = fd5509021 with 90-col fixes)
- f81101c5b first-order jumps in the preconditioner, clean_maxpic, iteration log (cc2)
- 06ccd4145 multiplicative x1/x2/x3 line preconditioner (cc3)
- dcaddd553 upwind boundary stress flux against the beam / bath ghost (cc4)
- 84cae4ad1 + 1e5ba13a4 a-posteriori low-order face limiter (clean_maxpos, default 3), clean_pc_sweeps (cc5)
Binaries: /viper/ptmp2/jinma/builds/bin/athena_none_cpu_<sha>_cc{1..5}; baseline athena_none_cpu_c20520ad_cbase.

## Results
(pending)

### Bitwise (staggered default) - PASS
- rw P1/3 tau10 N64 nt128 (gridA6 line 4): cbase (c20520ad5 = clean-1114) vs cc1 (fd550902), 34/34 tab files identical
  (/viper/ptmp2/jinma/clean_1114/cc/rw_bit_{base,cc1}).
- nltek1_FIN3_d x1, 60 cycles (stopped at the run's own wall limit), tab outputs + every CLEAN log line identical
  (/viper/ptmp2/jinma/clean_1114/cc/bitid).

### First cell-layout smokes (cc1/cc2)
- nltek1 (1-D): runs, every stage converges.
- cbg x0.25 (64x16, c dt/dx = 256, crossing beams in vacuum): BiCGStab with the x1-line preconditioner does not converge
  (4000 its) with PLM; first order converges at 130-800 its then also fails; 56 negative-E cells and 416 |F| > cE
  cells after the first BE step.  -> cc3: multiplicative x1/x2 line preconditioner; cc4: beam/bath ghost for the F flux.

### rw grid, amplitude 1e-6, cell layout (cc4 dcaddd55: plm + tau + maxpic 12 default), 1 period, nt = 2N
Runs /viper/ptmp2/jinma/clean_1114/cc/rw_grid (log rw_grid.log); staggered reference = rw/gridA6_27.out (cln27).
L1 error (rho+E) vs the exact eigenmode, order between N and 2N; [staggered]:

| P_rad/P_gas | tau | N64 | N128 | N256 | orders | staggered N64/128/256 |
|---|---|---|---|---|---|---|
| 1/3 | 0.1 | 5.65e-3 | 1.34e-3 | 2.97e-4 | 2.08, 2.17 | 5.65e-3 / 1.34e-3 / 2.98e-4 |
| 1/3 | 10 | 1.90e-3 | 4.60e-4 | 1.07e-4 | 2.04, 2.10 | 1.47e-3 / 3.40e-4 / 8.53e-5 |
| 1/3 | 1000 | 2.37e-2 | 1.19e-2 | 1.69e-4 | 1.00, 6.1 | 2.37e-2 / 1.18e-2 / 1.52e-4 (same gas-damping transition) |
| 3.3 | 0.1 | 3.01e-3 | 7.60e-4 | 1.88e-4 | 1.99, 2.02 | 3.00e-3 / 7.62e-4 / 1.90e-4 |
| 3.3 | 10 | 1.01e-3 | 2.67e-4 | 6.25e-5 | 1.92, 2.09 | 6.29e-4 / 1.64e-4 / 4.25e-5 |
| 3.3 | 1000 | 6.59e-2 | 3.37e-2 | 3.04e-4 | 0.97, 6.8 | 6.58e-2 / 3.37e-2 / 2.53e-4 |
tau 1000: the cell and staggered layouts give the same numbers; the N64/128 error is the damping rate (gas side), common to both.

Conservation (rw/conserve.py, periodic rw runs, max drift over 17 outputs, normalised by sum|.|):
cell layout energy 1.5e-16..1.7e-15, momentum 2e-14..7e-11; staggered (runs_gridA6_27, same cases) energy 8e-15..3e-12,
momentum 1e-11..5e-10.

### IDORT, first attempt (CC1 = cc4, clean_maxpic 4) - aborted
Step 12796548.183 on vipa1235 was SIGKILLed at 01:51 PST (all 25 runs).  Before that (first 1-13 stages):
cb x0.5 1579 negative-E / 2795 |F|>cE cells in the first stage, cbg x0.5 398-438 / 2380-2762, ws66 x1 2946 / 9616;
shd, shdg, ws66, cb: BiCGStab not converged (4001 its) with the x1-line preconditioner (shd: c dt/dx = 1e4).
nltek1/4 and hcyl: every stage converged, 0 negative, 0 |F|>cE.  nltek1 x0.5 (finished): error 2.24e-2 (staggered
L26 1.13e-2).
Small cbg/ws66 smokes with cc5 (low-order limiter, 2 sweeps): linear solves converge (1500-4500 its per stage summed
over 4 slope passes), but after 3 low-order passes cbg x0.5 still has 678-706 negative-E and ~1900 |F|>cE cells per
stage, ws66 x0.5 232 / 1118: the first-order a = c flux is not positive either in these 2-D beam transients.
Single normal beam: 0 negatives; single oblique (45 deg) narrow beam: 1 negative / 3 |F|>cE.

### rw time ladder, P 3.3, tau 10, N256, amplitude 1e-2 (cc4) - FAIL
| arm | successive differences nt 512/1024/2048/4096 | orders |
|---|---|---|
| cell plm (rw_lad) | 3.49e-7, 3.78e-7, 3.45e-7 | -0.11, 0.13 |
| cell first-order jumps (rw_lad_r0) | 3.51e-7, 3.57e-7, 3.20e-7 | -0.03, 0.16 |
| staggered (rw/runs_lad29) | 2.77e-8, 8.91e-9, 2.88e-9 (amp 1e-2 rows: 2.77e-7, 8.88e-8, 2.84e-8) | 1.64 |
Not the limiter (first-order jumps give the same floor).  Frequency converges (13.346464, 13.346067, 13.345977,
13.345966) but the damping rate does not (-3.65969, -3.65962, -3.65956, -3.65950: steps ~6e-5 per halving).
Picard converges (5-7 passes), no failed stages, no BE steps beyond the first.

### nltek (CC5 = cc5 + clean_maxpic 4 + clean_pc_sweeps 2, vipa1235; scored with cc/score_any.py)
| deck | x0.5 | x1 | x2 | x4 | staggered L26 x0.5..x4 |
|---|---|---|---|---|---|
| nltek1 | 2.237e-2 | 2.239e-2 | 2.245e-2 | 2.255e-2 | 1.127e-2 ... 1.126e-2 |
| nltek4 | 3.092e-3 | 3.094e-3 | 3.104e-3 | (running) | 1.549e-3 ... 1.550e-3 |
Resolution-independent in both layouts (closure / boundary-layer error vs the Eddington analytic profile); the cell
layout overshoots E in tau 1..0.005 by up to 18 % (staggered 11 %), top E/ex 0.973 (staggered 1.000).
0 negative / 0 |F|>cE stages.  hcyl x0.5: 1.150e-1 (staggered A16 x1: 2.14e-2), 0 negatives.

### Beam / shadow decks, CC5 (cc5, 2 sweeps, maxpic 4, low-order LF limiter 3 rounds) - stopped (02:1x CET)
| run | stages done | max neg-E cells/stage | max |F|>cE cells/stage | unconverged solves |
|---|---|---|---|---|
| cb x0.5 | 13 | 1288 | 2501 | 0 |
| cbg x0.5 | 13 | 706 | 1920 | 0 |
| cbg x1 | 1 | 4770 | 9084 | 0 |
| ws66 x0.5 | 13 | 246 | 1118 | 0 |
| shd x0.25 | 4 | 38 | 1414 | 0 |
| shdg x0.25 | 4 | 54 | 1748 | 0 |
Linear solves converge with 2 sweeps (1100-4500 BiCGStab its per stage summed over passes); every stage is
non-realizable, so every step goes to BE.

### cc6 (085045f2): final limiter chain (char flux on low faces, BE ends with all faces low) - smoke only
cbg x0.5: stages mark 10-13k faces over 3 rounds, BE with every face low-order still leaves 0-2 negative-E and
2500-3000 |F|>cE cells; ws66 x0.5: 20-30 negative, 2250-3000 |F|>cE.  Normal beam (1-D): 0 marks.
The remaining violation is consistent with the transverse stress rows (c D_bd E central with the rays' D) on
x1 faces, which the 1-D char flux does not cover.

## Diagnosis A: 2-D realizability (cbg x0.5, cc7 6c00ed0e, first step = BE with every face low-order)
Run dirs /viper/ptmp2/jinma/clean_1114/cc/smoke/diag_{false,true} (clean_char_all), log lines CLDIAG.
- char flux on the normal stress row only: 2582 cells |F|>cE after the all-low-order BE (2 negative-E).
  D - f f^T < 0 for the OLD state: 0 of 2582; for the rays' own flux h = H/J: 0 of 2582 (rays self-consistent).
  The old state there is vacuum (E_old = 1e-12, f_old = 0); the solution reaches |f| = 1.3-3.3 (e.g. i 1..12, j 0:
  E 6e-3..4e-2, f = (1.6..2.9, -0.5..0.1)), in cells where the rays carry no radiation (J/E ~ 0, filled closure
  D = (0.58, 0.41, +-0.25..0.49), C1 = 1).
  => fraction explained by an incompatible closure (D - f_old f_old^T < 0): 0 %.  The face flux produces it.
- char flux (D = 1 on the diagonal, 0 off-diagonal) on ALL stress rows of low-order faces: 610 cells (6 negative-E),
  and these are marginal: |f| = 1.001-1.003, f along x1 (f2 ~ 0), row j = 0, i 17..125.
  => the transverse stress rows c D_bd E (rays' D, central) cause ~76 % (1972 of 2582) of the violations; the rest
  is a 0.1-0.3 % excess along x1 (not diagnosed further; candidates: the boundary faces, which are never low-order
  and use the D-based ghost flux, and the tolerance of the 3 marking rounds).

## Diagnosis B: rw time floor (P 3.3, tau 10, N256, amp 1e-2, cc7 6c00ed0e), runs cc/rw_B_*
| arm | successive differences nt 512/1024/2048/4096 | orders | damping rate nt 512..4096 |
|---|---|---|---|
| hesdirk2 (baseline, = rw_lad) | 3.49e-7, 3.78e-7, 3.45e-7 | -0.11, 0.13 | -3.65969 -3.65962 -3.65956 -3.65950 |
| BE (ts=be) | 2.37e-4, 1.20e-4, 6.01e-5 | 0.99, 0.99 | -3.73574 -3.69767 -3.67857 -3.66900 |
| vfreeze=1 (advective a at v_old) | 6.01e-5, 2.94e-5, 1.44e-5 | 1.03, 1.03 | -3.63251 -3.64606 -3.65278 -3.65612 |
| vfreeze=2 (work term at v_old) | identical to baseline (bitwise same differences) | | |
BE is first order as it must be (its differences are far above the 3.5e-7 floor, so BE cannot show it).
Lagging the velocity in h chi a E makes the scheme first order (expected; that term must stay implicit), the work
term's dm^2/2 part has no effect.  Neither carries the floor.
| clean_alpha = one (RC off, full dissipation) | 3.48e-7, 3.77e-7, 3.45e-7 | -0.12, 0.13 | -3.65969 -3.65962 -3.65956 -3.65951 |
| clean_alpha = one, BE | 2.37e-4, 1.20e-4, 6.00e-5 | 0.99, 0.99 | same as BE |
At tau 10, N256 the cell tau is 0.04, so A(tau) ~ 0.98: the dissipation weight and the RC term are nearly inactive,
and switching them fully on/off changes the 7th digit only.  With first-order jumps (rw_lad_r0) the floor is the
same.  So the floor is NOT in the face flux (RC, dissipation, limiter), NOT in the work term, and BE is clean.
Remaining suspects (common to all arms): the once-per-step closure from the predicted state (test running:
clean_predict = false, cc/rw_B_np), and the stage coupling of the cell F with the hydro/IMEX machinery (the F slots of
t2inc/t2k and the Heun average of u0 at stage 2, which in the staggered layout acted on a derived F).
| clean_predict = false (closure from U^n, no prediction) | 3.17e-5, 1.59e-5, 7.93e-6 | 1.00, 1.00 | -3.6492 -3.6543 -3.65686 -3.65813 |
The closure time level carries the time error: lagged at U^n it is first order, 100x above the floor; with the
predicted closure (default) the differences sit at 3.5e-7 without converging.  Narrowed suspect for the floor: the
once-per-step closure from CleanPredictState (one D for both stages, built from a state predicted with the K1 slopes;
note the rays get iw F = u0 F of U^n, not a predicted F).  Next test (not run): closure rebuilt at each stage
time, or the predicted state's F also advanced with K1_F, to see whether the floor turns into order 2.

## Diagnosis A, step 2 (cc8 393aa19b): low-order faces carry the stress with D = I (now the definition, no key)
First step (BE, every interior face low-order at the end), |F|>cE cells (CLDIAG threshold 1e-6) / negative E:
| deck x0.5 | interior low-order only (clean_pos_boundary=false) | + physical boundary faces low-order (true) |
|---|---|---|
| cb_FIN3_r4 | 298 / 12 | 1372 / 6 |
| cbg_FIN3_g8x32 | 610 / 6 | 1374 / 6 |
| ws66 | 338 / 16 | 1352 / 6 |
Low order on the boundary faces (upwind characteristic flux against the beam/bath ghost) makes it worse: the new
violations sit at the beam injection (cbg i 1..7, j 6..7, |f| up to 1.27) where the imposed beam flux is replaced by
the ghost exchange.  Boundary faces stay out of the limiter (clean_pos_boundary default should be false; not yet
changed).  Marking tolerance (5 rounds) running.
Marking tolerance: clean_maxpos = 5 gives the same counts (cb 298, cbg 610, ws66 338; the final all-low-order solve
decides), and cbg with clean_tol 1e-14 gives the same 610: the remainder (|f| up to ~1.13 among the printed cells) is
a property of the fully first-order implicit scheme itself (dimension-by-dimension upwind with D = I), not of the
marking or the solver.

## Diagnosis B, step 2 (cc8): closure variants, P 3.3 tau 10 N256 amp 1e-2 (cc/rw_C_*)
| arm | differences nt 512..4096 | orders |
|---|---|---|
| clean_predict_flux (rays see F + dt K1_F) | 3.49e-7, 3.78e-7, 3.45e-7 | identical to baseline (periodic rays do not use F) |
| clean_closure_stage (stage 2 closure from Y1, no rates) | 1.30e-5, 6.45e-6, 3.18e-6 | 1.01, 1.02 |
Neither gives order 2.  Next (user, physical fix): the closure rebuilt from the rays of the current iterate inside the
Picard loop.

## Diagnosis B, step 3 (cc9 a1de69bf): clean_closure = iterate (closure from the rays of the iterate, every Picard pass)
Staggered layout, same binary, default closure (cc/rw_D_s33b): differences 3.47e-7, 3.77e-7, 3.45e-7, orders
-0.12, 0.13 -> the staggered layout has the SAME floor with the current code (the 1.64 of runs_lad29 was an older
build), so the floor is common to both layouts.
| arm (cc9, clean_closure = iterate) | differences nt 512..4096 | orders | formal solutions (nt 1024 run) |
|---|---|---|---|
| cell, P 3.3 tau 10 | 5.81e-7, 1.46e-7, 3.65e-8 | 2.00, 2.00 | 11385 for 1024 steps (~11 per step, ~5.6 per stage) |
| cell, P 1/3 tau 10 | 2.26e-7, 5.65e-8, 1.42e-8 | 2.00, 2.00 | 11195 |
| staggered, P 3.3 tau 10 | 5.80e-7, 1.45e-7, 3.62e-8 | 2.00, 2.00 | 10491 |
Default closure (predicted, 1 formal solution per step): floor 3.5e-7 in both layouts.
=> The time-order floor was the once-per-step predicted closure.  Closure from the rays of the iterate in every
Picard pass gives order 2.00 in time in both layouts, at ~11 formal solutions per step instead of 1 (each Picard
pass of each stage; clean_maxpic caps it).  The cell layout now defaults to clean_closure = iterate (commit
"the cell layout takes its closure from the rays of the iterate by default"); the K1 predictor stays only for
clean_closure = predict (staggered default, unchanged bitwise).

## A (physical fix: retarded formal solution) - NOT implemented, plan and code locations
The quasi-static rays are VetShortChar -> VetShortCharSweep (src/rad_m1/rad_m1_vet.cpp, ~3727-4310).  Per cell the
sweep reads a direction-independent opacity vc_(M1_VET_CHX) and source vc_(M1_VET_SRC) (set at ~3842/3884), and per
ray r integrates segments with dtau = 0.5 (c_up + c0) dx/|mu| and VetSegW(dtau, ex, w0, wu) (~4000, ~4077, and the
2-D/3-D upwind-interpolation branches ~4054-4115 that read CHX/SRC at the four upwind points).
Retarded version: CHX -> chi + 1/(c h) (direction independent, one line at ~3842/3884); SRC becomes per ray:
S'_r = (chi S + I_old,r/(c h))/(chi + 1/(c h)).  Needs (i) storage of the cell intensity of every ray from the last
accepted step (nray x cells; for the IDORT decks cbg x0.5 128x32 cells x nray), written after the step's final
sweep, (ii) every SRC read in the sweep (19 places, including the four-point upwind interpolations and the boundary
extrapolation 1.5 s0 - 0.5 s1) replaced by the per-ray S'_r, (iii) h = the stage step (dt_sub) passed in.
With clean_closure = iterate the rays then run inside the Picard loop with I_old of the step start.

## clean_pos_boundary default false (ba476c29e)

## clean_closure = extrap | extrap2 (52c5bb35b, binary athena_none_cpu_52c5bb35_cc11)
- One formal solution per step on the CONVERGED state after the last solve of the step (CleanClosureAfterSolve,
  state = u0 + gas after write-back, no prediction); its unnormalised moments M = (J, K = D J, H+, H-, <|I-S|>) are
  kept at t^n and t^{n-1} (cln_mn, cln_mp).  Stage closure from M(t_s) = M^n + (t_s - t^n)/(t^n - t^{n-1}) (M^n - M^{n-1}),
  then D = K/J, h+- = H+-/J, C1; guard per cell: J > 0, |H| <= J, K/J - h h^T PSD (all principal minors), else M^n
  (counted, "guard" on the CLEAN line).  extrap2: also the moments of the converged stage-1 solution; stage 2 on the
  line (t^n, M^n) - (t_1, M_1).  First step / restart: CleanClosure (state at t^n, no K1 rates); BE solves and the step
  after a BE step use M^n only.
- NOTE: for hesdirk2 both implicit stages end at t^{n+1} (A rows [1-g, g] and [1/2, b2, g] sum to 1), so c_1 = 1, not
  gamma: extrap uses M(t^{n+1}) for both stages, and extrap2's stage 2 takes the stage-1 moments themselves.
- nltek1 x1 (cc11, keys clean_maxpic 4, clean_pc_sweeps 2): error 2.239e-2 for extrap, extrap2 and iterate (same
  digits); formal solutions 61 / 120 / 135 (60 steps); guard fired 0 times; wall (login node, loaded) 117 / 155 / 190 s.

## clean_ray_retard (f9922d94c, binary athena_none_cpu_f9922d94_cc12)
- VetShortCharSweep reads chi and S through two per-ray helpers: chi' = chi + 1/(c h), S'_r = (chi S + I_old,r/(c h))/chi'
  (key off: + 0.0, same arithmetic).  I_old,r = the ray intensities of the last formal solution of the previous step
  (CleanRetardSetup swaps cln_iold/cln_inew at the first sweep of a new mesh time; first sweep of a run / restart:
  I = max(E + 3 F.Omega/c, 0) of U^n).  C1 keeps the physical S.  h = time from the step start to the stage end
  (CleanStageOffset: dt for both hesdirk2 stages and BE, since both stages end at t^{n+1}; dt_sub = g dt would
  overstate dI/dt by 1/g = 3.4).  FATAL with the banded multi-block sweep, a rotating set or vet_milne.
- Storage 2 x nray x (cells incl. ghosts) x 8 B: cbg/cb/ws66 x0.5 (128x32, 512 rays) 39 MB, x1 145 MB; a 3-D
  128^3 block with 512 rays 17 GB, 256^3 137 GB (not feasible for production without the moment reconstruction).
### First step (BE, all faces low order at the end), |F|>cE / negative-E cells, x0.5 decks, cc12 (cc/ret/G_*)
| deck | iterate, retard off | iterate, retard on | extrap (first step = state closure), retard on |
|---|---|---|---|
| cb_FIN3_r4 | 238 / 6 | 238 / 6 | 434 / 12 |
| cbg_FIN3_g8x32 | 538 / 6 | 528 / 6 | 760 / 10 |
| ws66 | 320 / 14 | 306 / 12 | 216 / 18 |
Retardation changes nothing material: in these decks c dt = 8 = 256 dx (box 4 x 1), so 1/(c h) adds an optical depth
of only 0.5 over the whole box; the light crosses the box twice per step and there is no radiation "in transit" at
the step scale.  The quasi-static rays are the right physics here; the violations are not a closure-time problem.
(iterate without retard vs the cc8 predict closure: 238/538/320 vs 298/610/338.)

## nltek1 vs EXACT transport (ana/nlte_exact.py; scorer cc/nlte_score.py)
| run | L1 vs exact | E/J_exact top cell |
|---|---|---|
| staggered L26 x1 | 9.52e-5 | 0.9999 |
| cell CC5 x1 (predict) / I11 x1 (iterate) / E11x x1 (extrap) / E11x2 x1 | 1.097e-2 / 1.098e-2 / 1.097e-2 / 1.098e-2 | 0.9725 |
| cell R0I x0.5 / R0E x0.5 (cc12) | 1.094e-2 / 1.094e-2 | 0.9727 / 0.9726 |
| cell + outgoing boundary stress c D E (branch clean-cc-bstress 2a68179c7, cc13x), extrap, x1 | 9.48e-5 | 0.9999 |
Profile (cell, iterate x1): E/J_exact = 1.000 for tau > 10, 1.023 at tau 1, 1.075 at 0.1, 1.086-1.089 for tau 1e-2..1e-4,
then 0.973 in the top cell.  Top cells: F/E = 0.48-0.49, top cell 0.46, escaping F 0.1265 vs staggered F/E 0.543,
F 0.1305 (escape 3 % low -> radiation piles up below).  The rays' boundary E flux is right (F_bdry/E_N = 0.541 = h+),
the closure is not the cause (same for every closure option).
Cause: the top face's STRESS flux.  Since dcaddd553 the open boundary face takes the upwind (a = c) stress between
the inner cell and the bath ghost: fY = c/2 (D E_N + D_g E_g) + c/2 (y_N - y_g).  With no inflow this is
c/2 (D_11 + f) E_N, exact only for D = f = 1; in the thin top (D_11 ~ 0.7, f ~ 0.54) it understates the outgoing
stress and pins the top cell's F to a value inconsistent with the half-range E flux (E dip in the top cell, F/E 0.46).
Fix tested: on an open face without inflow (no beam, zero bath) the stress is the outgoing half-range stress of the
inner cell, c D E (exact there: K = K+ when nothing comes in); beam / bath faces keep the upwind flux.  nltek1 x1 error
drops 1.10e-2 -> 9.48e-5 (= staggered), top E/J 0.9999, F/E 0.543.  Not on clean-cc yet (side branch; beam decks
and hcyl gates running).
Beam decks with the bstress experiment (cc13x, iterate, first step, cc/ret/H_*): |F|>cE / neg = cb 3770/14, cbg 3780/8,
ws66 3914/22 (vs 238/6, 538/6, 320/14 with the upwind boundary stress) -> 7-16x worse; the c D E stress on the
outflow faces of the beam decks (closure there often filled, J_ray/E ~ 0) is not acceptable as is.  A consistent
alternative would be the rays' outgoing half-range second moment K+ (needs a new ray sum); not done.

## rw time ladders, extrap / extrap2 vs iterate (cc11 52c5bb35, P 3.3 and 1/3, tau 10, N256, amp 1e-2; cc/rw_E)
| arm | P 3.3 diffs nt 512..4096 | orders | P 1/3 diffs | orders | formal solutions (nt 4096) | CPU s nt 4096 (login, loaded) |
|---|---|---|---|---|---|---|
| iterate cell | 5.81e-7 1.46e-7 3.65e-8 | 2.00 2.00 | 2.26e-7 5.65e-8 1.42e-8 | 2.00 2.00 | 35372 / 33673 | 3427 / 3082 |
| iterate staggered | 5.80e-7 1.45e-7 3.62e-8 | 2.00 2.00 | 2.26e-7 5.65e-8 1.41e-8 | 2.00 2.00 | 33583 / 33580 | 3196 / 2947 |
| extrap cell | 5.11e-7 2.74e-7 1.46e-7 | 0.90 0.91 | 2.73e-7 7.90e-8 2.52e-8 | 1.79 1.65 | 4098 / 4097 | 893 / 1315 |
| extrap staggered | 5.12e-7 2.75e-7 1.46e-7 | 0.90 0.91 | 2.74e-7 7.90e-8 2.52e-8 | 1.79 1.65 | 4098 / 4097 | 450 / 758 |
| extrap2 cell | 3.66e-7 2.46e-7 1.40e-7 | 0.57 0.82 | 2.65e-7 7.69e-8 2.47e-8 | 1.79 1.64 | 8398 / 8396 | 1436 / 1901 |
| extrap2 staggered | 3.67e-7 2.47e-7 1.40e-7 | 0.57 0.82 | 2.66e-7 7.69e-8 2.47e-8 | 1.79 1.64 | 8398 / 8396 | 777 / 1281 |
(ladcmp.py groups sorted by arm; formal solutions and CPU are P 3.3 / P 1/3.)  Guard fired 0 times in every run.
extrap: 1 formal solution per step (vs ~8 for iterate), 3.8-7x less CPU, but FIRST ORDER in time at P 3.3 (0.9) and
falling at P 1/3 (1.79 -> 1.65).  Diagnosis running (dbg-extrap branch, EXTDBG: extrapolated stage closure vs the rays
of the converged U^{n+1}).

## Physical boundary stress (c1a85e83a on clean-cc, cc14): outgoing half-range K+ of the rays + prescribed inflow
New ray sums K+^d_bd = sum_{mu_d>0} w I Omega_b mu_d (vet_cell slots 28..36, clean only), closure k+ = K+/J (CLN_KP,
filled/extrapolated like D); open boundary face stress = c k_out E_N (high face k+, low face D - k+) + the incoming
stress (beam column c sum w Omega_b mu1 I from ImplBeamInit, ibm_s; bath c E/6 normal).  E flux unchanged (h+- E +
beam/bath).  Touches rad_m1_clean_solver.cpp in ONE line (CleanMakeArgs: x.bstress = ibm_s).
| test | upwind ghost stress (before) | bstress c D E (cc13x) | physical K+ (cc14) |
|---|---|---|---|
| nltek1 x1 L1 vs exact / top E/J | 1.098e-2 / 0.9725 | 9.48e-5 / 0.9999 | 9.59e-5 / 0.9999 |
| nltek1 x2 | 1.104e-2 / 0.9724 | 9.49e-5 / 0.9999 | - |
| hcyl x0.5 (score_any) | 0.1154 (iterate and extrap) | 0.0468 | 0.0468 |
| cb / cbg / ws66 x0.5 first step |F|>cE / neg | 238/6, 538/6, 320/14 | 3770/14, 3780/8, 3914/22 | 3864/8, 3764/14, 3944/0 |
| beam1d (normal incidence) | 0 | - | 0 (3 solves) |

### extrap order loss: diagnosis (dbg-extrap branch c586665a, debug only, not for merge; cc/rw_dbgE, cc/rw_W)
- Extrapolated stage-2 closure vs the rays of the converged U^{n+1} (EXTDBG, max over the row, mean over steps):
  dD11 2.85e-8 / 7.3e-9 / 1.94e-9 and dH+ 9.2e-8 / 2.3e-8 / 6.1e-9 for nt 512/1024/2048 -> O(dt^2) as designed;
  dC1 8.6e-6 / 3.7e-6 / 1.4e-6 -> O(dt) (<|I-S|> is not linear in I), but C1 is inert here: clean_w_c0/c1 = 10/20
  gives the same differences to all digits.
- |E_extrap - E_iterate| (cell, P 3.3): 1.80e-6 / 7.40e-7 / 3.38e-7 / 1.62e-7 for nt 512..4096 -> O(dt), extrap2 the same;
  at the first dump (t = P/16) 2.3e-6 / 6.3e-7 / 2.0e-7 / 7.7e-8, growing to the final value: an accumulating O(dt^2)
  per-step error, not a start-up effect.
- Likely mechanism (not yet proven): the FSAL slope.  The flux is stiff (F relaxes to -(c/chi) grad(D E) within the
  step), so an O(dt^2) closure error moves the stage solution's F by O(dt^2), and the stored slope
  K = (Y2 - old vector)/(g dt) is then O(dt) off f(U^{n+1}; D(U^{n+1})); K1 enters the next step with a21 dt ->
  O(dt^2) per step -> first order.  iterate is immune because its stage-2 closure is D(Y2) itself (to tolerance).
  Candidate fix: after the end-of-step formal solution, re-evaluate K1 = f(U^{n+1}) with the new closure instead of
  the FSAL difference quotient (one operator application, no solve).

## beam45n (narrow 45 deg beam, 1 step, iterate): old upwind boundary (cc12) 3222 |F|>cE / 5 neg; physical boundary
(cc14) 15564 / 94, and its cycle-1 BE solve broke down (it -4001) and was accepted -> fixed below.

## Merge of clean-cc-pc (7e4c9bacb) + BE failure fix (9cf0c0aad), binary athena_none_cpu_9cf0c0aa_cc15
- Conflict only in rad_m1_clean_cell.hpp (bstress vs physical boundary): kept clean-cc's physical half-range form
  (pc did not touch that file).  CLEAN line: "... guard G solves S".
- A backward-Euler solve with nit < 0 now sets dtr_fail (dtr_why 1) when <time>/dt_retry is on and prints
  "CLEAN WARNING cycle N: backward-Euler solve not converged (...)"; stages keep failing through t2_fail.
  Check on beam45n with cc15 still running (cc/ret/K_beam45n.log).
- rw smoke with cc15 (cc/rw_M15): iterate time ladder P 3.3 cell 5.813e-7, 1.456e-7 (2.00), staggered 5.800e-7,
  1.448e-7 (2.00) = cc11 to 4 digits; space grid P 3.3 tau 10 N64/128/256 9.885e-4 / 2.707e-4 / 7.095e-5 = cc11
  iterate (rw_gridI 9.885e-4 / 2.707e-4 / 7.096e-5).
- rw space grid extrap vs iterate (cc11, amp 1e-6, nt = 2N): P 3.3 tau 10 1.037e-3/2.846e-4/7.48e-5 vs
  9.885e-4/2.707e-4/7.10e-5; P 1/3 tau 10 1.954e-3/4.97e-4/1.283e-4 vs 1.894e-3/4.774e-4/1.225e-4; tau 0.1 and 1000 equal.

## 11-15 (worker): task A, extrap order
### K1 re-evaluated with the end-of-step closure (b7efa9f62 on clean-cc, binary athena_none_cpu_b7efa9f6_cc16)
CleanCellClosureSlope(-1) before / (+1) after CleanClosureAfterSolve: K1 -= closure-dependent part of f (h div of E, F
rows; chi a E with a = v + D v), closure update, face speeds, K1 += new part.  rw ladder (cc/rw_F):
| arm | P 3.3 diffs 512..4096 | orders | P 1/3 diffs | orders |
|---|---|---|---|---|
| extrap cell + K1 re-eval | 3.57e-7 2.46e-7 1.40e-7 | 0.54 0.81 | 2.71e-7 7.85e-8 2.51e-8 | 1.79 1.65 |
|E_extrap - E_iterate| at nt 512/1024 (P 3.3): 1.54e-6 / 6.89e-7 (was 1.80e-6 / 7.40e-7).  => FSAL K1 is NOT the mechanism.
### Stage bisection (dbg-bstress branch, key clean_dbg_itstage, binary athena_none_cpu_740321de_dbgB3, cc/rw_G), P 3.3 cell
| arm | diffs nt 512/1024/2048 | order | |E - E_iterate| nt 512/1024/2048 |
|---|---|---|---|
| stage 1 closure from rays of the iterate, stage 2 extrap | 3.03e-7 7.53e-8 | 2.01 | 4.03e-7 1.01e-7 2.52e-8 |
| stage 1 extrap, stage 2 iterate | 4.05e-7 2.50e-7 | 0.70 | 1.48e-6 6.80e-7 3.26e-7 |
=> the order loss is in the STAGE-1 closure: stage 1 (theta-like, c1 = 1) needs the closure of its own state Y1, which
differs from U(t^{n+1}) by O(dt^2); the extrapolated D(t^{n+1}) is O(dt^2) from D(Y1) and that difference costs
O(dt^2) per step.  Stage 2 may use the extrapolated closure.
### Production-form predictor arm (coordinator): rays once per step on P = (stage-1 start) + dt K1 (E, F, gas momentum,
gas energy; no exchange relaxation), clean_closure = predict, key clean_dbg_pform (dbgB5 36f0b09c, cc/rw_H)
| arm | P 3.3 diffs 512..4096 | orders | P 1/3 diffs | orders |
|---|---|---|---|---|
| pform cell | 3.21e-6 1.60e-6 8.01e-7 | 1.00 1.00 | 2.16e-7 9.32e-8 5.04e-8 | 1.21 0.89 |
| pform staggered | 3.21e-6 1.60e-6 8.00e-7 | 1.00 1.00 | 2.16e-7 9.32e-8 5.04e-8 | 1.21 0.89 |
=> first order; same class as extrap (any closure O(dt^2) off D(Y1) in stage 1 loses the order).  K1 there is not
re-evaluated (that needs end-of-step rays, i.e. 2 solves per step).
### Cost of the 2.00 arm "stage 1 iterate + stage 2 extrap" (Gs1): formal solutions per step 7.7 / 6.7 / 6.2 for nt
512/1024/2048 vs iterate 13.2 / 11.1 / 10.2 vs extrap 1.0 (login node, loaded: CPU not comparable, not reported).
Next (running, cc/rw_J): stage 1 starts from the extrapolated closure and updates it from the iterate's rays only in
passes 1..n (n = 1, 2): 1 + n formal solutions per step.

## Task B: beam decks with the physical boundary stress - DIAGNOSIS (dbg-bstress branch, CLHIST/CLMAG, cc/ret/L_*, M_*)
- First step (BE, all interior faces low-order at the end), physical K+ (b0): bad cells ~3800 of 4096 in every deck,
  almost all gross (|f|-1 >= 0.1), spread over the whole box (most >= 16 cells from an x1 face).
  Ray-valid vs filled: cbg 2958+104 valid / 602+76 filled; ws66 3456 valid / 432 filled; cb 244 valid / 3494 filled
  (cb has only 272 ray-valid cells).  The x1max (outflow) inner cells are ray-VALID in cbg and ws66 (32/32), filled in
  cb (30/32).
- Mode 2 (upwind ghost stress only where the face's inner cell is filled): cb 284 (fixed), cbg 3762, ws66 3990 (not
  fixed).  => the hypothesis "violations come from filled-closure outflow faces" is REFUTED: in cbg/ws66 the outflow
  faces are ray-valid and the physical k_out there still drives ~3800 violations.
- Mode 3 (physical K+ in every high-order solve, upwind ghost stress in the solve with every interior face low-order):
  cb 238/6, cbg 538/6, ws66 320/14 = the old upwind numbers.  Mechanism: the all-low-order solve is the first-order
  (D = I, a = c) scheme; a boundary stress c k_out E with no y term is not its upwind closure (no dissipation), and with
  c dt = 256 dx its error spreads over the whole box in one step.  Stage solves (high order) are non-realizable with
  every boundary form (~3000 cells, t2s 2), as before.
- Implemented on clean-cc (b4dc0e75 = 18b70413 + comment fixes, binary athena_none_cpu_b4dc0e75_cc17); full runs CC17
  running (suite/*/CC17_x*, cc/ret/N_beam*).
- BE-breakdown fix verified (cc15, cc/ret/K_beam45n.log): "CLEAN WARNING cycle 1: backward-Euler solve not converged
  (dt_retry, 1 so far)" then "DT_RETRY cycle=1 ... dt 4e-3 -> 2e-3 reason=NON-CONVERGED".

## Momentum conservation (rw periodic), worst case P 3.3 tau 0.1 N256 nt512 amp 1e-6
- Debug CLMOM (dbg-bstress 09b97a42, key clean_dbg_mom): sum over ALL cells of gas m1 and radiation F1/c^2 at the start
  and end of every CleanSolve, long double sums (cc/momdbg/t14.log, clean_tol 1e-14, 30 cycles, cell layout, iterate):
  max |M_tot(t) - M_tot(0)| = 9.1e-21 = 5.5e-17 of sum|rho v| (1.66e-4 for the 4 rows): ROUND-OFF.  Stage increments
  of t2inc: gas + radiation sum to 1e-25..3e-25 (same weights on both; inc gas -1.0e-13 vs rad +1.0e-13).
- The 2e-14..7e-11 "drift" came from rw/conserve.py, which reads the tab outputs: these are slice_x2 = 0.125, ONE of the
  nx2 = 4 rows of the vet mesh, so the sum misses the momentum exchanged between rows through the x2 faces (transverse
  stress at solver-tolerance level: K12 5.3e-12, K14 2.0e-13, amp 1e-2 3.3e-16; cc/rw_K, cc17).  Not a
  non-conservative term.  (clean_tol 1e-10, cc/momdbg/t10.log, 20 cycles: 2.2e-20 = 1.3e-16 relative -> tolerance independent, round-off.)

## clean_closure = eddington (87e045577, binary athena_none_cpu_87e04557_cc18): D = I/3, h+- = +-1/4, k+ = I/6, C1 = 0,
no formal solution, same code path (CleanClosure).  Space ladder amp 1e-6, nt = 2N, ref=edd, clean_tol 1e-14,
clean_tol_nl 1e-13 (cc/rw_EDS; table cc/rw_EDS.table; L1 rho/E, orders rho/E):
```
EDSc   P 0.333 tau 0.1 | N64 5.626e-03/4.157e-03  N128 1.334e-03/1.179e-03  N256 2.958e-04/3.098e-04  N512 6.642e-05/8.048e-05 | 2.08/1.82  2.17/1.93  2.16/1.94
EDSc   P 0.333 tau 10 | N64 1.903e-03/1.142e-03  N128 4.711e-04/3.014e-04  N256 1.194e-04/7.784e-05  N512 3.003e-05/1.973e-05 | 2.01/1.92  1.98/1.95  1.99/1.98
EDSc   P 0.333 tau 1000 | N64 4.145e-03/4.035e-03  N128 8.056e-04/7.818e-04  N256 1.697e-04/1.663e-04 | 2.36/2.37  2.25/2.23
EDSc   P 3.33 tau 0.1 | N64 3.024e-03/7.016e-03  N128 7.663e-04/1.863e-03  N256 1.910e-04/4.741e-04  N512 4.808e-05/1.157e-04 | 1.98/1.91  2.00/1.97  1.99/2.03
EDSc   P 3.33 tau 10 | N64 9.328e-04/8.832e-04  N128 2.552e-04/2.318e-04  N256 6.694e-05/5.938e-05 | 1.87/1.93  1.93/1.96
EDSc   P 3.33 tau 1000 | N64 4.811e-03/4.753e-03  N128 1.195e-03/1.186e-03  N256 2.897e-04/2.888e-04 | 2.01/2.00  2.04/2.04
EDSs   P 0.333 tau 0.1 | N64 5.621e-03/3.788e-03  N128 1.332e-03/1.019e-03  N256 3.102e-04/2.954e-04  N512 6.286e-05/2.615e-05 | 2.08/1.89  2.10/1.79  2.30/3.50
EDSs   P 0.333 tau 10 | N64 2.418e-03/2.636e-03  N128 2.801e-03/2.501e-03  N256 3.785e-05/1.568e-04  N512 3.501e-05/4.307e-05 | -0.21/0.08  6.21/4.00  0.11/1.86
EDSs   P 0.333 tau 1000 | N64 3.826e-03/3.466e-03  N128 5.962e-04/5.626e-04 | 2.68/2.62
EDSs   P 3.33 tau 0.1 | N64 3.218e-03/6.278e-03  N128 8.804e-04/1.546e-03  N256 4.212e-04/1.090e-03  N512 1.010e-04/2.959e-04 | 1.87/2.02  1.06/0.50  2.06/1.88
EDSs   P 3.33 tau 10 | N64 4.922e-03/1.438e-02  N128 5.901e-03/1.348e-02 | -0.26/0.09
```
CAVEAT staggered: with D = I/3 the staggered Newton loop stagnates at clean_tol_nl 1e-13 ("it -1": BiCGStab 0 its,
Newton not converged after maxnl): N128 P 1/3 tau 10 has 38 failed solves of 606 -> failed stages are redone by
backward Euler and, since 9cf0c0aad, a failed BE solve triggers dt_retry (dt halved): the staggered rows above are
contaminated (erratic orders).  The cell layout has 0 failures (its tol_nl is max(tol_nl, 10 tol)).  Staggered order
runs need clean_tol_nl >= 1e-12 (or the same max(tol_nl, 10 tol) floor), and the time ladders at clean_tol_nl 1e-15
(cc/rw_EDT2) will hit it for staggered.  The farm matrix should use those tolerances.

## Task A RESULT: stage-1 closure pass (cc/rw_J, dbgB6 40df0f05, key clean_dbg_closn on top of extrap + K1 re-eval)
Stage 1 starts from the extrapolated closure and rebuilds it ONCE from the rays of its first iterate (pass 1); stage 2
uses the extrapolated closure; end-of-step rays as before.
| arm | P 3.3 diffs 512..4096 | orders | P 1/3 diffs | orders | formal solutions / step |
|---|---|---|---|---|---|
| n = 1 pass (cell) | 3.03e-7 7.53e-8 1.89e-8 | 2.01 2.00 | 2.27e-7 5.68e-8 1.43e-8 | 2.00 1.99 | 2.05 |
| n = 2 passes (cell) | same to 4 digits | 2.01 2.00 | same | 2.00 1.99 | 3.1 |
| iterate (cc11, rw_E) | 5.81e-7 1.46e-7 3.65e-8 | 2.00 2.00 | 2.26e-7 5.65e-8 1.42e-8 | 2.00 2.00 | 8.6-13 |
| extrap (cc11) | 5.11e-7 2.74e-7 1.46e-7 | 0.90 0.91 | 2.73e-7 7.90e-8 2.52e-8 | 1.79 1.65 | 1.0 |
| production prod-1118 b2d852bb (vet_sc, time2_vet predict) | 6.06e-7 1.43e-7 3.49e-8 | 2.09 2.03 | 2.26e-7 5.72e-8 1.44e-8 | 1.99 1.99 | 1.05 |
| unify-1114 7a3f4e23 | identical to prod-1118 (for vet_sc both default time2_vet = predict; extrap2 applies to vet_col only) | | | | 1.05 |
Implemented on clean-cc as the extrap default: dc67dba7b (binary athena_none_cpu_dc67dba7_cc19); extrap is now the
cell-layout default.  Not yet rerun with cc19 itself (tests moved to Orion): run lists in
/viper/ptmp2/jinma/clean_1114/farm_req/ (rw_cc19_time.txt, rw_cc19_space.txt caselists for run_radwave/rwfarm;
idort_cc19.run).  CPU timings were login-node, loaded: not reported.

## Task B accuracy, CC17 (b4dc0e75: physical K+ in high-order solves, upwind ghost stress in the all-low-order solve)
| deck | score_any CC17 | R0I (iterate, upwind boundary, cc12) | L23 staggered | CC5 |
|---|---|---|---|---|
| cb x0.5 | 1.677 | 1.677 | 1.471 | 1.095 |
| cbg x0.5 | 1.624 | 1.624 | 1.398 | 1.475 |
| ws66 x0.5 | 1.345 | 1.345 | 1.256 | 1.453 |
| hcyl x0.5 | 0.0468 | 0.1154 | | |
| shd x0.25 | 0.264 | | | |
| shdg x0.25 | 0.256 | | | |
| nltek1 x1 (vs exact) | 9.40e-5, top E/J 0.9999 | 1.10e-2 | | |
Beam decks: stages always non-realizable (max nreal 3100-3300 per stage, shd/shdg 6400-6600), every step goes to BE;
first-step BE counts = 238/538/320 (old upwind values).  No failed solves, no warnings.  beam1d/45/45n full runs were
killed by the coordinator before finishing.
Stage budgets (CLBUD, feedback off, dbgB7, clean_tol default 1e-10): relative residual of (dE_cells + h div boundary
flux - h(exchange)) per stage: cbg 6e-11..1.7e-10 (E, y1), hcyl 5e-12..1.8e-11, nltek1 1-2e-5 (E; stiff cells, scale
1e6): at the solver tolerance (feedback off = the gas is a fixed bath; with feedback the write-back makes energy and
momentum exact by construction, see CLMOM).
