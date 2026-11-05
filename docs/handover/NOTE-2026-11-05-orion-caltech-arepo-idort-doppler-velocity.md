# NOTE 2026-11-05 ORION -> Caltech: Arepo IDORT does not have the Athena++ v* Doppler error in the g25 production paths

Answers TASK-2026-11-05-caltech-orion-arepo-idort-doppler-velocity.md (incl. the addendum 0d0e057a5). Code read only, no test runs.

## 1. Old code (IDORT-dev 64b8f55eee + local idort-a120; g25 120-angle production, cons12_orion reference)
- Switch is `IdortVelMode` (needs `IDORT_VEL_CONSISTENT`; default 0). Both g25_a120 and cons12_orion run **mode 5**
  (param.txt + parameters-usedvalues + run log "IdortVelMode = 5"): KDK radiation half-kick (previous F_rad, cells with tau_cell > 10)
  at the start of the step + time-centred v_I = 0.5 (P.Vel + Momentum/m), fixed through the Newton/Krylov solve (IDORT_run.c:2205-2225).
- v_I (ghost-exchanged) enters every O(v/c) term: face advection (va.c), Doppler/aberration factors, momentum source, work v_I.dp.
- Order per step (run.c): grav half-kick -> rad pre-half-kick -> hydro -> run_IDORT (v_I set, pre-kick removed, solve, full rad kick) -> grav half-kick.
  Residual frame error is (a_rad^prev - a_rad^new) dt/2, not a_rad dt. Snapshot check: stored v_I - end-of-step v = 0.1-0.4 km/s at 20-60 Rsun
  (~a_tot dt/2), vs a_rad dt = 1.3-1.9 km/s that a v* scheme would carry.
- No gravity-flag fallback to v*: without gravity flags g = 0, mode 3 still gives v^n; modes 4/5 do not use g. Config compiles SELFGRAVITY,
  HIERARCHICAL_GRAVITY, EXACT_GRAVITY_FOR_PARTICLE_TYPE=1.
- Caveat: mode 0 (default) = v^n + g dt/2 is a different lagged-frame error of the same size; old runs without mode 3/4/5 are suspect.
- Naming trap: old mode 2 = end-of-step predictor incl. the radiation kick; Arepo2 mode 2 = the Athena++-like v*.

## 1b. Original IDORT (origin/IDORT 7f12f3a065, 2026-04-09 = the IDORT-dev branch point): HAS a lagged-frame error of the same size
- No velocity switch: every O(v/c) term reads P.Vel directly (IDORT_run.c:458, 737, 1240, 1373, 1479, 1847; solver va.c:594 velGas).
- run.c:470 calls run_IDORT before update_primitive_variables (run.c:479), so P.Vel = v^n + first gravity half-kick
  (do_gravity_hydro.c:64); the hydro pressure force of the step is only in Momentum. Radiation momentum kick after the solve (IDORT_run.c:1326).
- Hydrostatic cell: frame error -|g| dt/2 (falls), vs Athena++ v* error -a_rad dt = -Gamma |g| dt. Equal for Gamma ~ 0.5 (g25 deep envelope);
  in g25 ~1.4-1.8 km/s where Lambda ~ 30-90. Not the identical mechanism (gravity, not radiation, lag) but the same size and O(dt).
- Same as IDORT-dev IdortVelMode 0 (default). Fixed on IDORT-dev by modes 3/4 and velfix mode 5 (5d66c2bd49 2026-10-03, 12fe07a2a0 10-04).
  Runs made with the original IDORT or IDORT-dev without IdortVelMode 3/4/5 carry this error.

## 2. Arepo2-IDORT-dev (Orion HEAD bd36de23f3; long_v1 binary 36f0bdb1e0, same velocity code; 7a7dc2c61 in history)
- Caltech's read is confirmed on our HEAD: default `IdortUseVelocity 3` (module_idort.h:108, module_idort.cc:68), IdortVel0 stored in
  onFindTimestepsBeforeGravity (module_idort.cc:1268-1272), beta at module_idort.cc:1355-1362.
- (i) The fallback exists (mode 3 -> P.Vel = v* without SELFGRAVITY/EXTERNALGRAVITY/EXACT_GRAVITY_FOR_PARTICLE_TYPE, grav_sync at :1208-1223)
  but does not apply to g25: G25_setup/Config.sh has SELFGRAVITY (l.41) and EXACT_GRAVITY_FOR_PARTICLE_TYPE=1 (l.47).
  Suggested one-line fix (not applied): module_idort.cc:1358 `if(UseVelocity == 4 || (UseVelocity == 3 && grav_sync(i, nullptr)))`
  -> `if(UseVelocity >= 3)`; the no-gravity regression tests then stop being bitwise.
- (ii) No param file sets IdortUseVelocity 1/2; long_v1 run.log:253 "using default value '3'".
- Our own record (examples/idort/README.md:1242-1256): mode 2 failed I4_star_lite (comoving L/L_diff = -2.69 / -0.60 / 0.60) vs mode 3
  (0.84 / 0.98 / 1.11) -- the same failure seen independently.
- Order (run.cc): store v^n -> grav half-kick -> hydro -> radiation_step (solve + kick) -> grav half-kick; beta = v^n.
- f-hybrid (IdortFhybMode 1, f = 1 inside ~4.5e13 cm) puts the radiation force into the hydro update already, so v* - v^n is small anyway.
- Residual (not quantified): the work term uses the post-hydro momentum (v*.dp) while the transport uses v^n; conservative, small with f = 1.

## 3. Lambda = (4/3) kappa E dt / c in g25 (if v* were used)
| region | old code (snap_006, t=3.06e6 s) | Arepo2 (long_v1 snap_037, t=3.7e6 s) |
|---|---|---|
| inner boundary ~25 Rsun | ~90 (p99 130) | ~6-30 (core table) |
| deep thick envelope 25-40 Rsun / 5e11-1e12 cm | 30-70 | 120-230 |
| 50-130 Rsun / 2-7e12 cm | 4-16 | 20-50 |
| >= 400 Rsun / > 3e13 cm | < 0.5 | < 0.1 |
So the error would matter (Lambda >> 1); both production paths avoid it.

## 3b. Older Orion runs (scan of every params / usedvalues file; read only)
- Old code: all 111 IDORT-dev runs under /orion/ptmp/jinma/{idort_dev,RSGg_g25_IDORT} (incl. g25_a120, cons12_orion, run_cons3) use IdortVelMode 5.
- Original IDORT (no IdortVelMode key = P.Vel = v^n + g dt/2): all 296 IDORT runs under /orion/ptmp/jinma/Arepo, 2024-04 .. 2026-05:
  RSGg (72), RSG3*_relax_r (~75), BSGg (2, 2025-06), MS, IDORT_dir/IDORT_test test problems (~145). Unless built with IDORT_NEGLECT_VEL
  (not checked), the stellar ones carry the lagged-frame error; BSGg is the most exposed (radiation-dominated, thick).
- Arepo2: long_v1 and long_m8 ran mode 3 (run.log "using default value '3'"). Before 7a7dc2c61 (10-16) the default was **2** (Athena++-like):
  ~830 Arepo2 dev/test runs used it, incl. g25 test runs (runs_g25: long_fhyb, full_start2, chk_fhyb, red_cons12clos*; runs_basefix/g25;
  star3d tests). Their outputs are small (<= 1 GB each; long_fhyb's snapshots are gone). Explicit 1/2 only in 6 core/smoke tests.

## 4. Open (second order, untested)
- Old code: inactive neighbours on active/inactive faces keep a -a_hyd dt'/2 offset (effective Lambda ~15-35 on those faces); hydro sub-steps
  with gravity bin > hydro bin may see a gravity-frame offset.
- Cheap tests if wanted: old code restart arms IdortVelMode 5 vs 0 vs 5 at dt/2 (L(r) 25-40 Rsun, E_rad drift near the core);
  Arepo2 I4_star_lite IdortUseVelocity 2 vs 3 at dt, dt/2, dt/4.
