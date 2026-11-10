# NOTE 2026-11-10 ORION -> viper: raycheap17 FAST SUBSET (pulse_u, fs_u, lat_th_mid, mvr_u) + bisection

Answer (fast subset) to `TASK-2026-11-10-orion-raycheap17.md`. Binary /orion/ptmp/jinma/raycheap14_1109/athena_dea3e3732
(no rebuild). Work dir /orion/ptmp/jinma/raycheap17_1110 (scripts/ = your st17.py etc.; all_eval.txt = every eval line),
job 228210 (1 p.exclusive node, interactive steps). Python: st0925.py needs >= 3.9 (`dict | dict`); used the
python-waterboa/2024.06 module (3.12, numpy 1.26). Part 1b (other 5 sp cases) and Part 2 (cs) not run yet.

**Summary**
- **S09 C reproduces 09-25 on all 4 cases** (pulse_u E 1.93 2.06 2.08 / F1 1.87 2.05 2.07; fs_u 1.95 2.09 2.08 /
  2.00 2.07 2.07; lat_th_mid 2.01 2.00 / 2.00 1.88; mvr_u E 1.98 1.99 2.05, F1 2.39 2.28 1.98, T 1.99 1.98 2.01).
- **lat_th_mid and mvr_u: SP passes in C, X and T** (all variables p ~2 at the finest pair; lat_th_mid F1 1.83/1.82 as
  09-25's 1.88). mvr_u SP C/T E show a low first pair (-0.52, then 4.75, 2.06): the coarse level, same in every arm.
- **Two production losses on sp, both located:**
  1. **TIME: `implicit_flux_beam = halfrange` with vet_gd -> no time convergence (pulse_u).** SP T: E d 3.0e-2 3.6e-2
     3.8e-2 4.0e-2 (p -0.25 -0.08 -0.07): a dt-INDEPENDENT ~4 % difference between levels. Bisection on T:
     SP-col (vet_gd off) p 2.02 2.03 2.02; S09+gd p 0.99 1.37 1.59 (and 2.03 1.88 with 2 more levels, l = 32, 64);
     S09+gd+{lag,pos,frz} = S09+gd; S09+gd+x2 0.92 1.38 1.65; **S09+gd+fse = SP (flat 4e-2)**. Inside fse: dropping
     implicit_hr_*, implicit_thin_* or implicit_blend* alone: unchanged (flat); **dropping implicit_flux_beam
     (halfrange -> closure; blend and thin keys must go with it, FATALs) restores p 0.99 1.37 1.56** (= S09+gd).
     Not moved by: SP-norot (rotation off, log "every 0 cycle(s) (0 = fixed)") 0.24 0.04 0.16, SP-pred, SP-tolrel,
     SP-tight, SP-nc64, SP-nofse, SP-lag (bitwise SP), and the half-range options implicit_hr_face=cell, hr_q=col,
     hr_recon_pred=true, hr_recon_fresh=1, hr_recon=dc (all bitwise SP or 1e-6 off). fs_u (vet_gd dropped, see
     Part 0) passes T 2.03 2.03 2.02.
  2. **SPACE: `vet_col_ncore = 8` caps the finest pair (n 256 -> 512)**, as your 09-25 note said. fs_u SP C E 1.94
     2.08 **1.48**, X 1.96 2.09 **1.45**; pulse_u SP-col C 1.92 2.05 **1.04**, X 1.94 2.07 **0.94**. With
     vet_col_ncore = 64: fs_u SP-nc64 C **1.94 2.08 2.08** / F1 1.99 2.07 2.07; pulse_u SP-colnc64 C **1.92 2.05 2.08**
     / F1 1.87 2.04 2.07.
- So on these 4 cases: **SP with implicit_flux_beam = closure (no blend/thin) + vet_col_ncore 64 would be 2nd order**;
  vet_gd's time order is pre-asymptotic on the coarse T levels (needs l >= 16 to show 2).

## Part 0 (2-cycle SP runs) -- FATALs and what was dropped (scripts/drops.json, applied to every SP-based arm)
| case | FATAL | dropped / set |
|---|---|---|
| pulse_u | `implicit_vimp needs hydro with coupling, gas_feedback and dbg_gas_force on` | implicit_vimp = false |
| fs_u, lat_th_mid | `vet_gd: reflecting outer x1 not supported; use implicit_bc_x1max = marshak`, then the vimp one | vet_gd = vet_gd_replace = false, implicit_vimp = false (so SP = SP-col there) |
| mvr_u | the vet_gd reflecting-x1 one | vet_gd = vet_gd_replace = false |
SP-central also FATALs as listed: `implicit_flux_faces = all needs implicit_flux = berthon | blend` and `implicit_blend
= idort_f needs implicit_flux_beam = halfrange` -> the arm also sets implicit_flux_faces = x1, implicit_blend = tau_f.
Echoes: implicit_precond = mg, time2_vet_col = extrap2, implicit_tol_rel = 0.1 on all 4.

## st17.py changes (Orion, all marked "Orion" in scripts/st17.py; st17.orig.py = yours)
- a file lock around s9.mk + parse_deck (parallel st17.py calls rewrote the shared 09-25 decks: `KeyError: 'rad_m1'`
  and corrupt decks in the first launch, all rerun);
- DROPS (Part 0), SP-central faces/blend keys, an extra arm SP-colnc64, `S09+a+b` multi-feature arms, and two
  bisection hooks (env GDX = keys excluded from the S09+ set, EXTRA = ad-hoc key=val) with their own run dirs.

## Key eval lines (st17.py eval; d per pair | p)
| case arm mode var | d | p |
|---|---|---|
| pulse_u SP C E | 1.18e-1 9.17e-2 8.77e-2 7.77e-2 | 0.36 0.07 0.18 |
| pulse_u SP X E | 1.37e-1 9.40e-2 8.45e-2 6.27e-2 | 0.54 0.15 0.43 |
| pulse_u SP T E | 3.01e-2 3.58e-2 3.79e-2 3.97e-2 | -0.25 -0.08 -0.07 |
| pulse_u SP-col C E | 3.28e-2 8.64e-3 2.08e-3 1.01e-3 | 1.92 2.05 1.04 |
| pulse_u SP-col T E (l to 64) | 2.14e-2 5.29e-3 1.30e-3 3.19e-4 7.91e-5 1.97e-5 | 2.02 2.03 2.02 2.01 2.01 |
| pulse_u SP-colnc64 C E | 3.28e-2 8.64e-3 2.08e-3 4.93e-4 | 1.92 2.05 2.08 |
| pulse_u S09+gd T E (l to 64) | 3.87e-3 1.95e-3 7.54e-4 2.51e-4 6.13e-5 1.67e-5 | 0.99 1.37 1.59 2.03 1.88 |
| pulse_u S09+gd+fse T E | 3.33e-2 3.86e-2 3.99e-2 4.12e-2 | -0.21 -0.05 -0.05 |
| pulse_u S09+gd+fse minus beam/blend/thin T E | 3.87e-3 1.96e-3 7.57e-4 2.56e-4 | 0.99 1.37 1.56 |
| pulse_u SP-central T E | 7.80e-3 4.12e-3 1.59e-3 5.05e-4 | 0.92 1.38 1.65 |
| pulse_u SP-central C E | 3.47e-2 4.19e-2 4.92e-2 4.54e-2 | -0.27 -0.23 0.12 |
| fs_u SP C E / X E | 2.16e-2 5.65e-3 1.33e-3 4.78e-4 / 2.02e-2 5.18e-3 1.22e-3 4.47e-4 | 1.94 2.08 1.48 / 1.96 2.09 1.45 |
| fs_u SP-nc64 C E | 2.16e-2 5.65e-3 1.33e-3 3.16e-4 | 1.94 2.08 2.08 |
| fs_u SP T E | 2.13e-2 5.21e-3 1.27e-3 3.14e-4 | 2.03 2.03 2.02 |
| lat_th_mid SP C / X / T E | | 2.02 2.00 / 2.01 2.01 / 2.02 2.01 2.01 (F2 2.04 2.04 / 2.03 2.04 / 1.99) |
| mvr_u SP X E, T T | | 2.03 2.06 2.10; T 2.02 2.01 2.00 |
SP-col / -nofse / -lag on fs_u, lat_th_mid, mvr_u: = SP (vet_gd already off there; thin switch / lag inert).
mvr_u SP-central X: the finest level (n 512) jumps to d 1.6 (run ends normally at t = 10, 1928 steps; central flux at
the X-mode fixed dt, i.e. the largest c dt/dx): an SP-central artefact, not production.

## Next (Orion, continuing; user: sp 2nd order in space and time first, then cs)
Locate the halfrange-beam term that does not converge in dt (code: impl_beam_hr paths, rad_m1_implicit.cpp
~14529-14800 hrx branch of m1_impl_aphll, the hr sf build per stage), then Part 1b on the other 5 sp cases with
SP + beam closure + ncore 64, then Part 2 (cs). Also: an Orion side finding on the berthon x1 face coefficients under
implicit_recon_lag = step (frozen at the stage start: 1st order in time on sp with eddington; picard restores 2) --
on these 4 cases SP-lag = SP bitwise, so it is not what pulse_u shows.
