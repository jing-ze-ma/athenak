# TASK 2026-11-10 viper -> ORION: raycheap17 -- is TODAY's production M1/VET still 2nd order on the 09-25 sp suite? bisect any loss; cubed-sphere twins (CPU)

Answer with `docs/handover/NOTE-2026-11-xx-orion-raycheap17.md`. Core thread: implicit VET 2nd order in space AND time
on sp and cs with the production scheme. 09-25 (viper, m1-sp-order2b ce5120a65, tests_m1/runs_5q_sporder2b) showed
the sp wedge 2nd order in space-time on every case below with the then-new defaults (closure vet_col,
implicit_marshak_face linear, vet_col_order2, time2_vet_col predict, vet_col_fk_min 0, central face flux, recon_lag
picard / closure_lag pass). Production (HG / RSG cs, R136 sp) has since moved to vet_gd, blend + halfrange + thin switch
fse, recon_lag / closure_lag step, time2_vet_col extrap2, implicit_tol_rel 0.1, positivity / floor machinery, the
vet_gd rotation (+ EMA), mg precond. Your raycheap12-16: hesdirk2 2nd order on clean 1-D; 2-D Cartesian multi-D
defaults 1st order (picard + pass fixes it); cs transport-only plateau p 0.58 (raycheap16 pending). This TASK runs the
09-25 suite itself, with the production key set, and bisects.

## FAST SUBSET FIRST (user: results fast) -- run this, NOTE it, then the rest
Cases **pulse_u, fs_u, lat_th_mid, mvr_u**; arms **S09 (C), SP (C, X, T), SP-col, SP-nofse, SP-lag, SP-central**
(the four one-off arms in the mode(s) where SP < 1.8; if SP passes everywhere, run them in C only). Steps: Part 0 for
these 4 cases only, `st17.py probe`, then all of it in one job (every run is 1 rank, except lat_th_mid s = 8 at 4).
Send the NOTE as soon as this subset is done. The remaining sp cases/arms (steps 1-3 below) are Part 1b; the cubed
sphere (step 4) is **Part 2, held**: viper appends/confirms it after your subset NOTE.

## Binary: NO REBUILD, no gate
`/orion/ptmp/jinma/raycheap14_1109/athena_dea3e3732` (fork/raycheap2-1109 dea3e3732, PROBLEM none, MPI, CPU; your
raycheap14/15 build). viper checked: all 242 production keys in `prod_keys_rad_m1.txt` are read by dea3e3732
(`git grep '"rad_m1", "<key>"'`); the 09-25 tests (sph_shell incl. gas_vr / gas_teq / lat_amp; sph_atm; the cs
P_l(mu) lateral mode `lat_l`/`lat_axis`, `shell_report`, the cs user history) exist there. Production binary p13c
(bc35bfca8) = dea3e3732's prod-1112-cand base + rstabe (hesdirk2 auto_be state in the restart: restart-only) + the
opt-in floor cap (`floorcap_*`, not named by any production run.cfg): no difference for these fresh runs.

## Files (this branch, `docs/handover/scripts/raycheap17/`)
- `st17.py` driver (decks, probe, runs, evaluation); `st0925.py` = the 09-25 `st.py` VERBATIM (case_cfg, levels,
  restriction); `solib.py` (09-25 helpers; work dir from `RC17_W`), `common.py`, `bin_convert.py` (readers).
- `ref0925/<case>.resolved.txt`: the 09-25 binary's RESOLVED `<rad_m1>`/`<hydro>`/`<time>` per case (2-cycle runs of
  athena_new_none_cpu with the 09-25 decks, rst header). `ref0925/RESULTS_{radial,lateral,gas}.txt`: the 09-25 numbers.
- `prod_keys_rad_m1.txt`: the production `<rad_m1>` set = the HG restart header hegiant.00143.rst (effective input incl.
  the run.cfg command-line keys) minus star / unit / table / BC / IC / dt-limiter keys (EXCL list at its top: c_light,
  arad, opacity + kappas, e_floor, coupling, gas_feedback, force_reference, BCs + Marshak, ic_rad, implicit_cfl,
  cfl_rad, subcycle, dbg_gas_*, tsolve_opac, vet_scatter_kappa_e, I/O / GPU-only keys). Tag `cs` = cs-only
  (cs_vertex_mono, cs_cross_mono: CSP arms only). `vet_gd_rotate_dt` (production 1170 s = ~5-30 dt) is set per case
  to T/16. R136 (sp) differs from HG in implicit_tol 1e-8, opac_newton true, vet_gd_gl_nmu 10, rotate_avg 4,
  vet_gd_tangential false; with const opacity only the last three could matter -- not separate arms unless SP fails.
- viper check of the S09 recipe (09-25 binary, `st17.py deck <case> S09 C <lev>` vs the plain 09-25 deck):
  pulse_u n32, mvr_u n32 (m1 + hydro_w tab) and lat_th_atm s1 (bin data) **bitwise**.

Setup: `export RC17_W=/orion/ptmp/jinma/raycheap17_1110 RC17_BIN=<binary> RC17_MPI="mpirun --bind-to none"`
(st17.py adds `-np`), python3 + numpy; run from anywhere (absolute paths).

## Test inventory (all from the 09-25 case_cfg; sp wedge r = 1..r1, theta 1.17..1.97, phi 0..0.8)
| case | physics | r1, T | rho kappa | Kn | 09-25 C orders (E / F1 [/ gas]) |
|---|---|---|---|---|---|
| pulse_u | vet_col radiating pulse, Marshak top, e_out 1e-8, uniform r | 3, 1.0 | kappa_s 1 | 5 (thin) | 1.93 2.06 2.08 / 1.87 2.05 2.07 |
| pulse_s | same, stretched r (POLY) | 3, 1.0 | 1 | 5 | 1.61 1.98 2.07 / 1.53 1.97 2.08 |
| pulse_edd_u | same, closure eddington | 3, 1.0 | 1 | 5 | 2.01 1.99 2.00 / 1.91 1.98 2.00 |
| fs_u | free-streaming shell, reflecting walls | 3, 0.6 | 1e-3 | 5e3 | 1.95 2.09 2.08 / 2.00 2.07 2.07 |
| lat_th_atm | theta mode (lat_amp 0.5) on the T-S4 atmosphere (c 100, flux in, Marshak out) | 5, 0.08 | 100 r^-2 | 0.06 (thick) | 1.98 1.99 / 1.96 1.96 / F2 1.83 1.96 |
| lat_ph_atm | phi mode, same | 5, 0.08 | same | 0.06 | 1.98 1.99 / 1.99 1.97 / F3 1.71 1.94 |
| lat_th_mid | theta mode, rho kappa 1 shell, reflecting r (mirror top) | 3, 1.0 | 1 | 0.6 (thin) | 2.01 2.00 / 2.00 1.88 |
| mvr_u | moving gas v_r 0.01, kappa_s 1000 + kappa_p 10, gas carries E (vimp) | 3, 10 | 1000 | 5e-3 | 1.98 1.99 2.05 / 2.47 2.32 1.88 / T 1.99 1.98 2.01 |
| mvp_u | passive E advected by the v_r flow, kappa_s 1000 | 3, 10 | 1000 | 5e-3 | see ref0925/RESULTS_gas.txt |
Levels: radial n = 32 64 128 256 512 (4 x 4 angular, 1 rank); lateral s = 1 2 4 8 (nx1 = n_lat = 16 s; 4 ranks at
s = 8). Kn = 1/(rho kappa_t L), L = pulse width 0.2 / lateral half wavelength ~0.8 r / atmosphere top.

Modes (st17.py):
- **C** (space-time, the 09-25 measure): dx and dt refined together; L1_V(q_n - R q_2n)/L1_V(|q_n - <q_n>|).
- **X** (space only): every level at dt_fixed = 0.9 dt_C(finest)/4.
- **T** (time only, fixed grid; radial nx1 = 128, lateral s = 2): dt_fixed ladder halving over l = 1 2 4 8 16 (radial
  coarsest = 0.9 x the 09-25 `_T` level-1 dt (implicit_cfl x8, hydro cfl 0.8); lateral coarsest = 0.9 x 4 dt_C(s=2)).
- Every arm except S09 runs the dt_fixed ladder (`time/dt_fixed`, `time/dt_retry=false`, `time/cfl_number=0.95`,
  `rad_m1/implicit_cfl=1e6`, as raycheap12); dt_retry is named false on S09 too (today's fresh default is on).
  `st17.py probe` first: the S09 decks for 1 cycle -> `$RC17_W/dt17.json` (first dt per case/level). **Check every run
  ends at tlim with cycles = ceil(tlim/dt_fixed) (+1 round-off step)**; a run whose dt drops below dt_fixed is invalid.

## Arms
| arm | keys | dt |
|---|---|---|
| S09 | 09-25 deck + every 09-25 resolved `<rad_m1>` key named at its 09-25 value | 09-25 recipe (C only) |
| S09f | same keys | dt_fixed ladder (recipe-neutrality check) |
| SP | 09-25 case deck (geometry, BC, opacity, IC) + prod_keys_rad_m1.txt; pulse_edd keeps closure eddington (vet_gd* / vet_rotate* dropped, thin switch none) | ladder |
| SP-col | vet_gd=false vet_gd_replace=false (vet_col proper) | ladder |
| SP-nofse | implicit_thin_switch=none | ladder |
| SP-lag | implicit_recon_lag=picard implicit_closure_lag=pass | ladder |
| SP-central | implicit_flux=central implicit_flux_beam=closure implicit_thin_switch=none implicit_thin_corr=none | ladder |
| SP-upw | implicit_enthalpy=upwind (production plm; rad_m1/reconstruct and hydro/reconstruct are plm = default in every arm, so no arm for them) | ladder |
| SP-nopos | implicit_pos_floor / pos_gas / pos_floor_solve / pos_floor_cons / realisable_coupling = false, time2_posfix=0, time2_auto_be=0, time2_accept_fallback=false, time2_k1_redist=0, time2_pred_cap=0 | ladder |
| SP-nofrz | implicit_thin_freeze=0 (dt-dependent mask c dt kappa_P < 1e-2; mvr_u has kappa_p 10) | ladder |
| SP-tolrel | implicit_tol_rel=0 (production 0.1, cap 1e-4: a stop RELATIVE to the step change) | ladder |
| SP-tight | tolrel 0 + implicit_tol 1e-13, lin_tol 1e-14, lin_cnorm 1e-14 | ladder |
| SP-pred | time2_vet_col=predict (production extrap2, 09-25 predict) | ladder |
| SP-norot | vet_gd_rotate_every=0 vet_gd_rotate_dt=0 vet_rotate_avg=1 vet_rotate_avg_ema=0 | ladder |
| SP-nc64 | vet_col_ncore=64 (production 8; 09-25: ncore 8 capped pulse_u at 512) | ladder |
| S09+gd / +fse / +lag / +pos / +tol / +x2 / +frz | S09 + the production values of that feature's keys (PLUS in st17.py) | ladder |
| CS09, CS09f, CSP, CSP-<f> | the same sets on the cubed-sphere twins; CSP also names the production cs mesh keys cs_seam_flux shared, cs_seam_blockend extend, cs_vertex_3cell true, layout_invariant true, cs_seam_rho_guard 4 | as above |

Cubed-sphere twins (full sphere, 6 panels, meshblock nx2 = nx3 = 8, ranks = blocks capped at 48): radial **pulse_cs,
fs_cs, pulse_edd_cs, mvr_cs, mvp_cs** (N = 8 per face, nx1 as the sp levels; measure: the `SHELL` per-radius column
mean of E and F1 from `problem/shell_report`, restricted in r as on sp, plus the lateral spread max (Emax-Emin)/Emean,
which must not grow with n); lateral **latatm_cs** (sph_atm + `lat_l=2 lat_axis=z`), **latatmv_cs** (axis = cube
vertex (1,1,1): the mode crosses the vertices), **latmid_cs** (sph_shell, rho kappa 1, reflecting): N = 8 s, nx1 = 16 s,
s = 1 2 4 (8 optional: 3.1M cells; run it if s = 4 takes < 30 min at 48 ranks); measure: the P_2 mode amplitude
A = (VEP - VE VP/V)/(VPP - VP^2/V) from the user history at T (relative difference between levels) and the SHELL
radial mean. T mode on cs: direct L1 of the bin fields (same grid, block order). mvr_cs / mvp_cs: check the t = 0
hydro_w dump has v_r = 0.01 sin(...) (gas_vr is not sp-guarded in the pgen); if not, drop them and say so.

## Steps (priority order)
0. **Part 0 (minutes).** `st17.py probe`. Then for every case: `st17.py deck <case> SP C <lev1>` and `... CSP C ...`,
   run each for `time/nlim=2`; list every FATAL / refused key. Drop only what a FATAL demands (edit the arm, not the
   key file; list exactly what was removed per case). Quote the echo lines of implicit_precond (mg on 4x4 angular
   blocks), implicit_thin_switch, vet_gd_rotate_dt, time2_vet_col, implicit_tol_rel.
1. **S09 C, all 9 sp cases** (control). Must reproduce ref0925/RESULTS_* (orders within +-0.05, L1 to 2 digits). If not,
   STOP that case and report the diff: dea3e3732 then changed the 09-25 scheme itself (bisect with forward S09+ arms
   is meaningless until this is understood). S09f C on pulse_u, lat_th_atm, mvr_u (recipe check: same orders).
2. **SP C, X, T, all 9 sp cases.**
3. **Bisect** every (case, mode) with SP p < 1.8 at the finest pair in any variable (whatever its Kn class: 09-25 had
   ~2 everywhere): run all SP-<f> arms in that (case, mode) (cheap: run them together); the arm(s) that restore p >= 1.8
   name the culprit. If none restores alone: run the S09+<f> forward arms on the same (case, mode) and report which
   feature takes S09 below 1.8; then the pair of the two worst. SP-tight first wherever a T ladder flattens with
   d < 100 x implicit_tol.
4. **CS09 C + CSP C, X, T on the 8 cs twins**; bisect CSP failures with the CSP-<f> arms as in 3 (and, only if no
   rad_m1 arm restores it, CSP with each of the 5 cs mesh keys at its default, one at a time; say which).
Keep everything CPU-small: 1 rank per radial sp run, 4 at lateral s = 8, <= 48 per cs run; RAM << 40 GB per job.

## Pass criterion (user's relaxed rule, vet-acceptance-criteria)
p >= 1.8 at the finest pair where Kn < 0.3 (lat_*_atm, mvr, mvp and their cs twins); p > 0 (converging) where
Kn >~ 0.3 (pulse*, fs, lat_th_mid, latmid_cs). st17.py `eval` prints PASS / FAIL and flags thin-case orders < 1.8
as a loss vs 09-25 (to bisect, not a fail).

## Report (NOTE)
Per (case, arm, mode): the `st17.py eval` lines (d per pair, p, verdict) for every variable; Picard mean / max /
NON-CONVERGED per level; the hesdirk2 / positivity end lines that are non-zero; wall per level; cycles vs
ceil(tlim/dt_fixed). Part 0's refused-key list. A summary table case x mode: S09 p, SP p, culprit arm(s). Dirs + job ids.
No code changes; if a fix is obvious, describe it in the NOTE, do not commit.
