# NOTE vcpred-1010: predicted vet_col / vet_gd tensor under hesdirk2, and hesdirk2 on AG Car

Branch `vcpred-1010` (local, NOT pushed), worktree /raven/u/jinma/ATHENAK/wt-vcpred. It is based on
thinsw-1010 61cab6c0, which already contains torder 8ee1e217 and the thin switch.

- **Commit:** b2377eea `time2_vet_col_source = relaxed | state`, plus a vet_gd_async warning.
- **Binaries** (bin/):
  - athena_cpu_9d8fd8c9 and athena_he_a100_9d8fd8c9 are the base. Their code equals 61cab6c0, copied from thinsw_1010.
  - athena_cpu_b2377eea (md5 88a54e72...) and athena_he_a100_b2377eea (md5 ef3cf09a...) are the new build.

## Verdict

1. **Nothing new was needed to predict the vet_col tensor: the key exists and works for vet_col and vet_gd.**
   - The vet_col analogue of `time2_vet_sc` already exists: `<rad_m1>/time2_vet_col = lag | predict | rebuild | extrap` (m1-sp-order2b).
   - On the sph wedge its default is already `predict` for a fresh hesdirk2 run.
   - It builds ONE formal solution per step at P = U^n + dt K1 (Time2VetColAt). time2_vet_sc was in fact cloned from it.
   - vet_gd rides in the same build: VetColBuild -> VetLatBuild -> VetGdBuild. Its source (VetLatSweep(0)) reads the swapped-in iw EN/TP/KT and opac at P.
   - So predict covers vet_col, vet_col_lat and vet_gd, with no extra formal solution per step.
   - **Key choice: keep `time2_vet_col`. I did not extend time2_vet_sc or add time2_vet_col_pred.**
     - The key is per-closure, as time2_vet_sc is.
     - It is already in production inputs: agcar_A_xtf names `time2_vet_col = predict`, which is ignored under be.
     - It is already the sp default, restart-recorded, and validated in runs_5q_sporder2b.
     - A second key would only duplicate it.
2. **Measured T orders (CPU): hesdirk2 + predict = 2.0 for vet_col and vet_gd; lag ~1.0; be ~0.9.** Two interactions break this:
   - **(a) The relaxed VET source (default `vet_col_source = relaxed`) is first order in time with predict wherever gas-radiation exchange is resolved and out of equilibrium.**
     - Cause: the relaxed source relaxes P by a further dt (it lands near t^{n+2}).
     - Fix: new key `time2_vet_col_source = state` (b2377eea), read only when named. It gives 1.9-2.0.
     - **But on AG Car `state` is unusable.** Picard rises from 4.4 to 14.3 per solve (max 47), the cost goes to 2.3x, and the photosphere changes by O(1). This is the vet_scatter `ma` + current-J amplification by kappa_P/kappa_R. So keep `relaxed` for production. The order loss only matters in resolved non-equilibrium exchange.
   - **(b) `vet_gd_rotate_every = N` rotates the ray set every N CYCLES.**
     - The ray-noise pattern therefore depends on dt, and the T ladder becomes erratic (orders 0.1-1.5).
     - It does not grow with dt: the error is bounded by the ray-noise amplitude. Its size is about 7x predict's own time error but still about 10x below be's (pulse_u_T, nmu 4).
     - This is by design (temporal ray-effect averaging) and was not changed.
3. **hesdirk2 runs on AG Car A (4 A100, gpudev) with the production keys and no fix.** rc 0, no FATAL, NON-CONVERGED 0, no NaN, min E 1.43 as in be. Per arm:
   - **predict (+ vimp + vstage):** 0.618 s/cycle vs be 0.420, i.e. 1.47x at the same dt (dt is hydro-limited, cfl 0.3). Picard 4.37 per solve vs be 4.33, 2 solves per step.
   - **lag:** 1.33-1.40x.
   - **Inside 0.95 R_ph:** a dt/2 self-difference estimate puts predict's time error 5-10x below be's (E, dens, F1 in 0-0.8 R_ph). That gain comes from vimp + vstage + the second-order stages. The tensor predict changes the interior by only ~1e-7.
   - **Photosphere (0.95-1.05 R_ph):** the hesdirk2 dt-self-difference is NOT smaller than be's (dens 5.3e-3 vs 2.3e-3). The thin-switch keys do not change this.
4. **Default bitwise:**
   - 150 tab files from the new vs base CPU binaries (lag / predict / be / vet_gd arms, two cases) are identical.
   - AG Car `be` from the new GPU binary is bitwise equal to the base binary's (bin + hst).
5. **Restart-exact (CPU, restart at cycle 15 of 32, mid rotation block):** bitwise for vet_gd + hesdirk2 + predict + rotate_every 10 + twin/twin_fuse/twin_det + vet_scatter + vimp/vstage, with both source = relaxed and source = state, and for vet_col predict.

## 1. Code (b2377eea)

- **`time2_vet_col_source = relaxed | state`** (rad_m1_time2.cpp Time2Init and Time2VetColAt; rad_m1.hpp `t2_vcsrc`).
  - Read only when named; default `relaxed` (bitwise).
  - With `state`, the predict / rebuild build sets `vcol_srelax = false` for that build only. The source is then eps B + (1 - eps) E of P, or VetScatterMix with J = E_P.
  - It applies to vet_col, vet_col_lat and vet_gd (all read vcol_srelax). lag / extrap / backward Euler are untouched.
  - Echoed in the Time2Report line as `source=state`.
- **vet_gd_async + hesdirk2 + time2_vet_col = predict | rebuild:** a startup WARNING (rad_m1_vetgd.cpp).
  - With async, the tensor of cycle n is build n-1, i.e. D ~ t^n, so the run is first order in time.
  - It is a warning, not a FATAL, because predict is the sp default and a FATAL would break existing inputs.
- **Checked, with no change needed:**
  - `vet_source_noesrc` only acts under be.
  - VetGdRealDiag is a diagnostic after each solve.
  - twin / twin_fuse / twin_det (the fused twin reads only this build's source; restart-exact).
  - `vet_col_lat_offdiag = operator` (D_r,lat comes from the build, through VetLatOdMax).
  - The per-shell exact halo: vet_gd's lateral band is always exact, and `implicit_halo_direct` is unaffected.
  - The rotation angle comes from ncycle, so a rebuild's second build uses the same angle.
  - `vet_col_lat_every > 1` with predict is already FATAL.
  - **Not tested:** MPI-split vet_gd under predict (the CPU tests use one rank; the GPU runs use 4 ranks / 4 blocks and are deterministic).

## 2. Order tables (CPU, T ladders: dt/1..16 on a fixed grid; self-convergence L1; last three orders)

Battery pulses / rw with closure = vet_col are **not allowed**: vet_col FATALs on a periodic x1 ("a periodic x1 boundary has no top"), and vet_gd needs a 3-D wedge with a Marshak top. I used two substitutes:

- **(A)** The battery pulse on the Cartesian slab with outflow + Marshak x1 ends (`opulse_k*`, cpu/ana/od.py, nx1 = 256).
- **(B)** Uniform-r sp-wedge cases from runs_5q_sporder2b, on the 3-D wedge phi = pi/3 (cpu/ana/vst.py: NT 128 x 4 x 4, `_T` ladders). Arms:
  - `gd*` add the AG Car vet_gd keys: replace, GL 4x12 staggered, twin + fuse + det, noesrc, vet_scatter ma/relaxed, offdiag = operator, op_stencil.
  - `R` = rotate_every 10.
  - `vimp` + `vstage` are on wherever the gas moves.
  - `rswm_u_T` is the rw analogue: radiation-dominated moving slab, P_rad/P_gas ~ 10, gas-coupled, vimp, Marshak top with a bath.
  - `neq1_u_T` is a new case: static gas out of radiative equilibrium (aT^4 = 0.41 vs E = 0.1), rho kappa_P = 1, Marshak top.

### Cartesian vet_col (A): E / F1

| case | hesdirk2 lag | **hesdirk2 predict** | be |
| --- | --- | --- | --- |
| opulse_k0.128 | 1.41 1.97 1.85 / 1.61 1.89 1.70 | **1.39 2.00 2.01 / 1.61 1.98 1.91** | 0.78 0.71 0.55 / 0.77 0.73 0.63 |
| opulse_k12.8 | 1.18 1.58 1.42 / 1.36 1.49 1.37 | **1.32 1.97 1.99 / 1.53 1.98 2.01** | 0.76 0.73 0.57 / 0.74 0.72 0.64 |
| opulse_k1280 (levels to 64) | 0.57 0.37 0.81 0.96 (E) | erratic 4.45 -2.62 6.70 -0.15 1.50; L1 at levels 8-64: 4.9e-6 5.5e-6 1.9e-6 | 1.01 1.01 1.00 |

- **k1280 is thick.** D ~ 1/3 except at the walls, and the dt-differences sit in the top boundary layer (x > 0.98).
- Levels 1-4 of BOTH hesdirk2 arms have negative solve E at the Marshak walls and one stage fallback to BE.
- Predict's L1 at fine dt is 5-10x below lag's (2e-6 vs 1e-5 at 32-64), but the ladder is pre-asymptotic.

### sp wedge (B): E (F1, T)

| case | lag | **predict** | predict, source = state | be |
| --- | --- | --- | --- | --- |
| pulse_u_T vet_col | 1.09 1.06 1.03 | **2.06 2.05 2.03** | 2.06 2.05 2.03 | 0.76 0.88 0.94 |
| pulse_u_T vet_gd | 1.13 1.08 1.05 | **2.06 2.05 2.03** (F1 2.04 2.04 2.02) | same | 0.77 0.88 0.94 |
| pulse_u_T vet_gd, rotate 10 | 0.86 1.99 0.33 | 1.27 0.71 1.49 (L1 1.98e-3 vs 2.84e-4 unrotated; be 2.45e-2) | same | - |
| fsm_u_T vet_col (thin shell, Marshak top) | 1.07 1.04 1.03 | **2.00 2.03 2.02** | same | 0.75 0.81 0.88 |
| fsm_u_T vet_gd | 1.18 1.11 1.06 | **2.04 2.03 2.02** | same | 0.83 0.89 0.93 |
| fsm_u_T vet_gd, rotate 10 | 1.32 0.43 1.09 | 1.14 0.11 1.19 | same | - |
| rswm_u_T vet_col (vimp + vstage) | 1.99 1.96 1.92 | **2.00 1.98 1.96** (T 1.94 1.95 1.95, M1 2.0) | identical | 0.96 0.98 0.99 |
| rswm_u_T vet_gd (rotate 0 or 10) | 1.99 1.96 1.92 | **1.99 1.96 1.92** | identical | 0.96 0.98 0.99 |
| **neq1_u_T vet_col** | 1.93 1.90 1.92 (T 1.94 1.96 1.95) | **1.71 1.57 1.40 (T 1.67 1.47 1.23)** | **1.92 1.90 1.93 (T 1.95 1.98 1.99)** | 0.73 0.83 0.91 |
| **neq1_u_T vet_gd** | 1.93 1.94 1.89 | **1.72 1.60 1.38 (T 1.68 1.49 1.24)** | **1.91 1.94 1.89 (T 1.95 1.98 1.99)** | 0.73 0.83 0.91 |

Notes on the table:
- **rswm is thick** (rho kappa_s 1000): D ~ 1/3, so lag is already 2. It shows the gas-coupled vimp + vstage hesdirk2 = 2 vs be = 1 (this is the rw role).
- **neq1:** lag's relaxed source at U^n happens to land at t^{n+1} for the exchange, so lag looks second order there. predict + relaxed overshoots by a further dt; predict + state fixes it.
  - `predG` (vet_col_source = gas, all builds) gives the same as `state`: 1.92 1.87 1.91.
- **Rotation:** the ray set changes every 10 cycles, a dt-dependent pattern. The ladders are erratic at the ray-noise level (pulse_u_T vet_gd L1 at dt/16: 2.0e-3 rotated vs 2.8e-4 fixed, be 2.5e-2).
  - With nmu 10 (AG Car) the noise is smaller.
  - A rotation keyed to simulated time would still jump, which is first order with an amplitude of ~delta D.
  - The only way to restore order 2 with rotation would be a smooth rotation in t; not done.

## 3. GPU: AG Car A (agcar_A_vc.athinput = agcar_A_xtf + named keys)

Setup: 1 node x 4 A100, gpudev. Runs are 30 cycles at production dt (1252.6 s; tlim 3.75e4), and the H arms are 60 cycles at cfl 0.15 (dt/2). Arm keys on top of the production input:
- be = production;
- hesdirk2 adds `time_scheme = hesdirk2`;
- P = `time2_vet_col = predict`;
- L = `lag`;
- v = `implicit_vimp = true` + `time2_vstage = true`;
- ts = `implicit_blend_xthin = 0`, `implicit_thin_switch = fse`, `implicit_thin_corr = hcs`;
- S = `time2_vet_col_source = state`.

All runs: rc 0, no FATAL, NON-CONVERGED 0, no NaN / inf in log, hst or dumps, min E 1.43 (hps 0.85). The cost column is s per 1e4 s of simulated time.

| arm | s/cycle | cost | Picard mean / max | BiCGStab inner |
| --- | --- | --- | --- | --- |
| be (production; base / new binary bitwise) | 0.399-0.420 | 3.2-3.35 | 4.33 / 9 | 3.37 |
| be + ts | 0.407 | 3.25 | 4.43 / 10 | - |
| hesdirk2 L (vimp off) | 0.524 | 4.18 | 3.56 / 9 | 2.93 |
| hesdirk2 L v | 0.557 | 4.45 | 3.54 / 9 | 2.93 |
| **hesdirk2 P v** | **0.618-0.621** | **4.93** | 4.37 / 9 | 2.84 |
| hesdirk2 P v ts | 0.608 | 4.85 | 4.44 / 10 | 2.86 |
| hesdirk2 P v S | 1.415 | 11.3 | **14.3 / 47** | 2.89 |
| be dt/2 (beH) | 0.367 | 5.86 | 3.65 / 9 | 3.08 |
| hesdirk2 P v dt/2 (hpH) | 0.478 | 7.63 | 3.03 / 8 | 2.64 |
| hesdirk2 L v dt/2 (hlH) | 0.471 | 7.52 | 2.74 / 7 | 3.10 |

- **Cost:** there are 2 stage solves per step.
  - Predict costs ~19 % more Picard per solve than lag at production dt (11 % at dt/2), i.e. +11 % wall.
  - Run-to-run spread of s/cycle is ~5 % (be 0.399 vs 0.420).

### Field differences

Metric: max over (theta, phi) of |a - b| per shell, divided by the shell max; then the max over r/R_ph bands 0-0.8 | 0.8-0.95 | 0.95-1.05 | > 1.05. The dt/2 rows are time-error estimates: the first-order error is about 2x that row, the second-order error about 4/3x.

| pair | dens | m1_e | m1_f1 | velx |
| --- | --- | --- | --- | --- |
| hesdirk2 L v - be | 2.9e-4 3.1e-5 4.0e-3 1.1e-2 | 3.1e-4 7.1e-5 1.0e-3 9.5e-5 | 3.2e-2 2.1e-3 1.1e-4 5.3e-5 | 0.19 0.010 1.3 0.17 |
| hesdirk2 P v - be | 2.9e-4 3.1e-5 8.1e-3 1.1e-2 | 3.1e-4 7.4e-5 1.6e-3 8.5e-5 | 3.2e-2 2.1e-3 1.2e-4 8.8e-5 | 0.19 0.011 2.6 0.17 |
| hesdirk2 P v - hesdirk2 L v (dt) | max 6.8e-8 overall | 6.9e-8 | 1.3e-6 | 2e-2 (outer) |
| hesdirk2 P v - hesdirk2 L v (dt/2) | 1.2e-7 2.1e-7 1.6e-3 8.7e-4 | 1.1e-7 2.3e-6 3.1e-4 5.9e-5 | 1.1e-7 2.6e-6 1.0e-5 2.1e-5 | 1e-7 2.7e-4 0.09 0.011 |
| **be: dt - dt/2** | 1.3e-4 1.0e-5 2.3e-3 1.6e-2 | 1.8e-4 2.6e-5 6.2e-4 2.5e-3 | 2.3e-2 9.0e-4 1.0e-4 2.6e-3 | 0.05 0.006 0.49 0.11 |
| **hesdirk2 P v: dt - dt/2** | **1.1e-5** 1.3e-5 5.3e-3 **9.0e-4** | **2.9e-5** 4.1e-5 1.2e-3 2.6e-3 | **2.5e-3** **1.3e-4** 1.3e-4 2.6e-3 | 0.03 0.009 0.60 0.14 |
| hesdirk2 P v ts: dt - dt/2 | 1.1e-5 1.3e-5 5.9e-3 8.5e-4 | 2.9e-5 5.7e-5 1.8e-3 2.6e-3 | 2.5e-3 1.3e-4 1.7e-4 2.6e-3 | 0.03 0.011 0.74 0.14 |
| be ts: dt - dt/2 | 1.3e-4 1.0e-5 2.3e-3 1.6e-2 | 1.8e-4 2.6e-5 5.8e-4 2.5e-3 | 2.3e-2 9.0e-4 1.0e-4 2.6e-3 | 0.05 0.006 0.33 0.10 |
| hesdirk2 L, vimp off: dt - dt/2 L v | 1.4e-4 1.4e-5 3.7e-3 7.1e-4 | 2.8e-4 4.1e-5 1.3e-3 2.6e-3 | 4.6e-2 4.6e-4 1.3e-4 2.6e-3 | 0.33 0.008 0.77 0.16 |
| hesdirk2 P v S - hesdirk2 P v (dt/2) | 9.7e-5 1.2e-5 **1.3e-1** 2.3e-2 | 1.9e-4 2.6e-3 **0.77** 0.54 | 1.8e-2 2.1e-3 0.47 0.56 | 0.18 0.07 **83** 0.9 |

- **Interior:** at its production dt, hesdirk2 + P + v has a smaller time-error estimate than be at dt/2 in the interior (0-0.8 R_ph):
  - dens 1.5e-5 vs 1.3e-4;
  - E 3.9e-5 vs 1.8e-4;
  - F1 3.3e-3 vs 2.3e-2.
  It does this at lower cost per simulated time (4.93 vs 5.86).
- **At 0.8-0.95 R_ph** the two schemes are comparable.
- **At the photosphere (0.95-1.05)** hesdirk2's self-difference is 2x be's (dens 5.3e-3 vs 2.3e-3; velx 0.60 vs 0.49). So the first-order be error estimate (2x) and hesdirk2's (4/3x) are similar there (~5e-3 vs 7e-3 dens). The thin switch does not change this.
- **Source of the interior gain:** vimp + vstage + the stages (`hesdirk2 L, vimp off` is like be in the interior: F1 4.6e-2).
- **Caveat:** these are 30-cycle runs from the IC (seed transient); they are an estimate, not a converged error.

## 4. Recommendation

- **Use as the hesdirk2 production set for vet_col / vet_gd:**
  - `time_scheme = hesdirk2`, `time2_vet_col = predict` (the sp default), `implicit_vimp = true`, `time2_vstage = true` (the sp default);
  - keep `vet_col_source = relaxed` and `time2_vet_col_source = relaxed` (default);
  - no vet_gd_async.
- **It is robust on AG Car** (0 NON-CONVERGED, no fallbacks, restart-exact), at 1.47x the be cost per step at the same dt.
- **Do NOT make it the AG Car production scheme yet.** At hydro-limited dt it buys accuracy in the interior (5-10x lower time error, cheaper than halving be's dt) but not at the photosphere, where the AG Car science is. Before switching:
  1. Run a longer A run (>= 200 cycles, past the seed transient) with dt / dt/2 for both schemes, and compare the photospheric L_top, Mdot_top and v_r errors.
  2. Find why the photosphere time error does not drop. Suspects:
     - the dt-dependent `implicit_blend_xthin = 30` weight of the production hr flux (torder: hr is first order in T). ts with xthin = 0 did not help in 30 cycles, but on vet_gd fse reduces to hj;
     - the opacities and EOS at the stage-start state (torder "remaining issues");
     - rotate_every 10 ray noise at the photosphere.
  3. Only then flip `time_scheme` in the AG Car input.
- **`time2_vet_col_source = state`:** keep it opt-in for problems with resolved non-equilibrium exchange and no strong kappa_P/kappa_R. With vet_scatter `ma` on AG Car it is unstable / costly (Picard 14, photosphere O(1) change).
- **Rotation:** with `vet_gd_rotate_every`, hesdirk2 has a time-error floor at the ray-noise level. If second order matters more than ray-effect averaging, use rotate_every 0, or a larger nmu.

## Runs / scripts

- **cpu/:**
  - ana/vst.py: sp T ladders (copy of runs_5q st.py; arms lag / pred / be / predG / predS / gd*, cases rswm, neq; phi pi/3).
  - ana/od.py: Cartesian opulse_k* with the vcL / vcP arms.
  - sp.sbatch, cart.sbatch.
  - sp/st/<case>/<arm>_<l>, cart/runs/.
  - rst/rstest.sh (restart tests gdpR_base, gdpRS_new, pS_new: 0 diffs).
  - Jobs: 31060057 (e1), 31060148 (c1), 31060212 (e2), 31060244 (c2), 31060463 (e3, new binary + bitwise).
- **gpu/:**
  - gpu_arms.sh (ARMS / TL / BH env); inp/agcar_A_vc.athinput; ana.py, ana_gpu.py.
  - Runs: runs/j.31060000 (be hl hp, 30 cycles), j.31060253 (be beH hp hpH, tlim 3.75e4), j.31060474 (be[new] hps hpsH), j.31060475 (hl0 hpt hlH), j.31060854 (hptH bet betH).
- **build.sh:** job 31060211; CPU + A100 he_star_m1 builds of b2377eea.
