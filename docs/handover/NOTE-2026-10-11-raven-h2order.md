# NOTE h2order: hesdirk2 time order below 2 (rlat, rswv, P = 10 waves)

Worker H, 2026-10-11. Work dir /raven/ptmp/jinma/h2order_1011.
Branch `h2order-1011` (worktree /raven/u/jinma/ATHENAK/wt-h2order), based on prodfinal-1011 8279c0a6. Not pushed.

- Commits:
  - c6f48bb7: opt-in `time2_start = slope`, and the test hook `problem/gas_hold_v`.
  - cad4a6ce: opt-in `implicit_vimp_bc`; `time2_start_frac` default 1.
  - c38eb350: the probe keeps `implicit_f0_init`.
- Binaries (bin/): athena_cpu_s1 ... s7 and athena_he_cpu_s6 / s7. **s7 = c38eb350 is the final build.**
- Harness:
  - ana/h.py: run, rundeck and tord, on the same metrics as sptestsC cmpmk.py and ladder od.py.
  - ana/ev.sh, plus the diagnostics ana/{modes,ph,hist,hist2,prof,bdiff}.py.
  - pool.sbatch / pool2.sbatch / rst.sbatch. rst.sbatch does a pre-run at dt/8, then restarts every level from the
    same state; it isolates the start-up error from the per-step error.
- Bases: base/*.athinput. These are copies of the sptestsC rlat_T / rsw_T / rswv_T HFfM_1 decks, the suiteC ladder
  rw_p10t{1,100} tJePvs decks, and rsw_cart with gas_vr.

## 1. Verdict

There are three independent defects. None is in the conservative source (worker G).

1. **The first step is Lie split.** A hesdirk2 run with no stored slope does Heun hydro, then one BE solve. Its O(dt)
   start-up error survives as a free wave mode.
   - The error is O(dt) and does not accumulate: it rotates between v and rho with the wave phase (ana/hist.py).
   - The final-dump T ladder picks v at its maximum, so this is the P = 10 "v and F1 order 1.0-1.2, E 0.6" finding.
   - It is also half of rlat's 1.5 / 1.3.
   - Fix: `time2_start = slope`, opt-in.
2. **implicit_bc_advect and implicit_vimp do not match at a Marshak end face.**
   - With implicit_bc_advect (default on for fresh runs), the Marshak x1 end face carries the enthalpy flux a E.
   - Its a is the end cell's a at the old-vector velocity.
   - implicit_vimp gives every interior face the stage's radiative kick dv, but skips physical faces ("physical
     faces carry no enthalpy flux", a stale comment).
   - The missing kick is an O(dt) lag in the E row of the top cell. This is a per-step error, so it is first order.
   - It is the dominant error of rlat (lateral wave) and of rswv / Marshak-top + exchange.
   - Fix: `implicit_vimp_bc = true`, opt-in. Reduced AG Car A is **bitwise** with it (see section 4).
3. **A residual floor remains in rswv (radial background flow).** It is about 1e-5 of the pulse amplitude, at
   order 0.2-0.7 at the finest pair.
   - The pre-relaxed restart ladder (RV) shows the floor directly.
   - It is not the enthalpy limiter (central is identical), the tolerance, the predict cap, the closure or the thin
     switch (section 3.2).
   - **Open.**

Not defects:
- The P = 10, tau = 100 residual first-order term at amplitude 1e-5 is the linear-solver stop.
  - The stop is lin_tol relative to max|b| = the background E, with Eisenstat-Walker.
  - With implicit_tol 1e-14 / lin_tol 1e-15 every variable is 2.00.
- The face coefficients frozen for the solve (recon_lag step / thin switch) are not a lag in the production
  half-range path.
  - Rebuilding them at Picard pass 1 from the iterate was bitwise identical on all cases (code reverted).
- Rotation (vet_gd_rotate_every 10) is irrelevant for rlat and rswv.
- Hydro reconstruction (dc vs plm) is irrelevant.

## 2. Orders before / after

These are self-convergence T ladders, read from the last dump. The norms follow sptestsC (L1/mean|q| for bin, and
|q - mean| for tab) and od.py (P = 10).

| case | arm | key variables, orders (1-2 2-4 4-8 [8-16]) |
|---|---|---|
| rlat_T (lateral, vet_gd) | prodfinal default (H0) | dens 1.54 1.34; vely 1.51 1.30; E 1.54 1.33 |
| | slope only | 1.27 1.14 (all) |
| | **slope + vimp_bc** | **dens 1.94 2.16; vely 1.94 2.05; velz 1.89 2.08; E 1.94 2.26; F3 1.89 2.07**. L1(4-8) 1.7e-6 vs 8.3e-5 before |
| rswv_T (v_r 0.01) | default | dens 1.11 1.02; E 1.08 0.97; v 1.07 0.96 |
| | slope + vimp_bc (rotate 0 or 10) | dens 1.20 0.98; E 1.53 0.72; v 1.54 0.60. L1 4-8x lower (E 9.2e-6 vs 6.6e-5); floor, open |
| rswv_cart (vet_sc twin) | default | 1.1-1.2 / 0.95 |
| | slope + vimp_bc | E 1.63 0.74; dens 1.27 1.00 (floor as above) |
| Marshak top + exchange (rswv_cart, T_gas out of eq., no flow) | default start | 1.09 1.04 (all) |
| | pre-relaxed restart + vimp_bc (RU) | **E 1.95 1.83; F1 2.04 1.82; dens 2.43 1.73; v 2.04 1.81** |
| rw_p10t100 (P 10, tau 100) | default (B) | v 1.08 1.17 1.02; F1 1.10 1.18 1.03; E 2.32 2.27 1.04 |
| | slope | v 1.95 1.93 1.69; F1 1.95 1.92 1.68; E 1.97 1.89 1.62; dens 1.06 1.24 1.00 (tolerance floor) |
| | **slope, tol 1e-14 / lin 1e-15** | **E 2.01 2.01 2.00; F1 1.99 1.99 1.98; dens 2.02 2.05 1.99; v 1.99 1.99 1.99** |
| rw_p10t1 (tau 1) | default (B) | E 2.00 0.63 0.66; F1 1.35 1.11 1.08; dens 1.97 1.92 1.55; v 1.74 1.63 1.41 |
| | slope | **dens 1.99 2.06 2.00; v 1.96 2.02 1.98**; E 1.79 1.73 1.29; F1 1.58 0.92 1.13 (open, small: E L1 1.9e-5) |

Regression with `time2_start = slope` (s7):

| test | default | slope |
|---|---|---|
| rw_t10a T | velx 1.91 1.87 1.89 | 1.99 1.97 2.03 (others equal) |
| rw_t10a C 32-256 | as suiteC | identical to 3 digits |
| rw_t1000a T | 1.9-2.1 | 1.9-2.0, errors equal |
| rw_t1000a C 32-256 | as suiteC | 1.97-2.00, errors equal |
| eq3t (hesdirk2, moving equilibrium) | E 1.136724, P_tot 3.002079 | E 1.136750, P_tot 3.002068 |
| rsw_T (static pulse) | dens 2.02 1.90; v 1.73 1.23 | dens 2.03 1.76; v 1.78 1.29 (+ vimp_bc: equal) |

- rsw's v / F1 decay of order at 4-8 is the same with and without the fixes. It is not one of the three cases.
- **Default path bitwise** (s7, no new keys): rw_t10a T_1, rw_p10t100 l1 and rlat_T l1 are byte-identical to the
  prodfinal / suiteC binaries.

## 3. Details

### 3.1 time2_start = slope (rad_m1_time2.cpp Time2StartProbe, driver.cpp)

- On a hesdirk2 step with no K1, the probe runs first.
  - It is a BE solve of dt_p = time2_start_frac x dt from U^n (radiation, faces and gas; no hydro), with
    t2_solve = BESTORE.
  - This stores K1 = (Y_p - U^n)/dt_p.
  - The probe then puts U^n back (FluidRefresh for the gas) and the step runs the two stages.
- **The frac default is 1.** A tiny dt_p is wrong wherever U^n is off the stiff manifold.
  - Examples are a pulse with F = 0, or a moving gas in a radiation field at rest.
  - There K1 ~ jump/dt_p, and the stage-1 old vector U + (1-g) dt K1 overshoots by dt/dt_p.
  - Measured with frac 1e-3: floors in rsw and rlat. eq3t became wrong (P_tot 3.043) through a separate bug, fixed
    in c38eb350.
- For a start that is on the manifold and non-stiff, frac 1e-3 is better.
  - Marshak + exchange: 2.2 with frac 1e-3, 1.3 with frac 1 (FINXu).
  - The per-step order is 2 either way (RU).
  - This is a start-up trade-off, not a scheme defect.
- The probe's implicit_f0_init (the cold-start face flux F0 = F - a E of a moving gas) is copied into the face backup.
  The init therefore survives the put-back.
  - c6f48bb7 lost it: the faces started at 0, and eq3t gave P_tot 3.043.
  - A redo in the stage-1 solve also failed. It takes a at the old-vector velocity, which gives an O(dt)-wrong face
    start: P = 10 v went to first order with error 1.3e-3.
- The probe is reported on the Time2Report line (`time2_start=slope probes=N`).
- MHD: FATAL (the field restore is not wired).

### 3.2 implicit_vimp_bc (rad_m1_implicit.cpp ImplicitVimpBuild, kernel m1_vimp_j)

- With the key, the x1 Marshak end face with implicit_bc_advect takes the interior faces' vimp terms:
  - with the end cell's own dv^k, D_dd and P row;
  - ef = E^k of the cell on outflow, or the bath E on inflow;
  - on Cartesian (dt/dx) and on sp (A_f dt/V).
- What showed it:
  - rswv_cart with gas T out of equilibrium and no flow (`Xu`) is first order.
  - The error sits in the top 2-3 cells.
  - It is 2nd order with a reflecting top.
  - It is 2nd order with dbg_gas_force = false.
  - marshak_q is irrelevant: q is h+ of the formal solution there.
  - bmom_half changes its size.
  - A restart from t = 0.5 is still first order, so the error is per-step.
- rlat has the same top-cell physics: the lateral wave's gas is pushed against the reflecting top.
- AG Car A, reduced (480x32x32, 60 cycles, FRESH deck agc/inp/A0 = prodfinal agcar_A_fresh0 + output dcycle):
  - vimp_bc is **bitwise** (all bins).
  - The top cells either flow in with E_bath = 0, so ef = 0, or do not move in 60 cycles.
  - time2_start = slope changes the run by velx 6e-5 and F1 9e-4 (max rel; the first step only).
  - s7 default is bitwise vs athena_he_cpu_18ff91c8.
  - Runs: agc/runs/j.31082955.

### 3.3 Open: the rswv floor

- RV = rswv_cart (v 0.01): pre-run to t = 0.5 at dt0/8, then levels 1-8 to t = 2, slope + vimp_bc.
  - L1 9.5e-6 8.0e-6 7.2e-6 (E): no convergence.
  - RU (the same without flow, gas out of equilibrium) converges at 1.8-2.0.
- Identical to 3 digits with any one of:
  - implicit_enthalpy = central (no limiter);
  - tol 1e-14 / lin 1e-15;
  - time2_pred_cap 0;
  - closure eddington;
  - thin switch none;
  - hydro dc;
  - hydro hlle.
- implicit_enthalpy = upwind is cleanly FIRST order (2.9e-4 1.4e-4 7.0e-5). That is a separate first-order term in
  the donor-cell enthalpy path under hesdirk2 + vimp, not production (plm is the default).
- implicit_bc_advect = false (RVnb) and vimp_bc off (RVnv) both bring back a clean first-order term of 2-4e-5.
- The floor sits in the pulse region (r 1.9-2.3), not at the top.
- Candidates not yet tested:
  - the vimp dv^k Picard lag (the plm/central deferred enthalpy correction is taken at the previous iterate, and
    Picard stops at 2-4 passes on |dE|/E only);
  - the t2 "DA" enthalpy split.

### 3.4 Not used / reverted

`time2_face_refresh` (rebuild the frozen face coefficients at pass 1) had no effect, bitwise, and was not committed.

## 4. Recommendation

- `implicit_vimp_bc = true` is a clear bug fix (vimp skipped a face that carries enthalpy flux) and is bitwise on
  reduced AG Car A.
  - Recommend making it the default on fresh runs (restart without the key keeps false).
  - It stays opt-in until the user decides.
- `time2_start = slope` (frac 1) is a start-up option for convergence studies and test problems.
  - Production runs restart, so it only affects a FRESH run's first step.
  - Keep it opt-in.
- Test ladders on the IDORT / Jiang linear waves at amplitude 1e-5 need implicit_lin_tol <= 1e-15. Otherwise the
  Eisenstat-Walker stop makes a first-order-looking floor.

## 5. Files

- ana/: h.py, ev.sh, modes.py, ph.py, hist.py, hist2.py, prof.py, bdiff.py.
- runs/<case>/<arm>_<lev>:
  - B = suiteC baseline (symlinks);
  - H0 / Hr0 = prodfinal default;
  - S* = slope;
  - FIN = s7 final (slope + vimp_bc);
  - X* / RU / RV* = the Marshak + exchange and rswv bisections;
  - N* = gas_feedback false. **Invalid as a diagnostic:** the E-T Schur solve moves T' although the gas energy is
    not written back.
- jobs_*.txt (job lists), logs/.
- agc/ (the reduced AG Car script and decks).
