# icgen-1011: a transport-consistent radiation IC for AG Car A (he_star_m1)

Branch `icgen-1011` in `/raven/u/jinma/ATHENAK/wt-icgen`, based on icrad-1011 7045bc14.
Work dir `/raven/ptmp/jinma/icgen_1011`.

## 1. Problem

`ic_agcar_A_ge.txt` (make_ic_agcar_ge.py, Viper, not reachable) is an inward MLT++ column.
Below R_ph its E = a T^4 exactly. Above R_ph it is an isothermal (T_ph = 9050.8 K),
gravity-only atmosphere, with E "diluted" and flux factor F/(cE) = 0.25 -> 0.74. The
builder also used the TOPS kappa_R. The run uses the ext2 Ferguson/TOPS pair, which
differs below log T = 4.2, i.e. for r > 0.98 R_ph.

Prototype (proto/rt1d.py, proto/run1.py) with the RUN's kappa_R:

* E_sph1d / E_col is 0.95–0.97 at 0.98–0.995 R_ph and 0.38–0.48 at 1.005–1.5 R_ph.
  Deep down, E_sph1d - E_col tends to a constant (-0.12 % at 0.9 R_ph, tau 3000): the
  photospheric kappa mismatch shifts E by a constant at depth.
* E_sph1d / E(cycle 50, relaxed default run j.31064135/be) is 1.002 at 1.02, 0.999 at
  1.05, 0.997 at 1.1, 1.000 at 1.2, 1.007 at 1.5, 1.014 at 2.0 and 1.020 at 2.9 R_ph.
* The 0.98–1.01 R_ph layer is still moving hydrodynamically at cycles 50–200
  (rho changes by -20 % at R_ph, and E200/E50 is 0.9–1.18 there). It is not a target.

## 2. Design (A): opt-in pgen solve `problem/he_ic_rad = file | sph1d`

`file` is the default and is bitwise the old IC. Keys are read only when named.
The solve runs inside `HeStarM1` after the file is resampled to the fine grid and the
opacity tables are loaded, and BEFORE the column goes to the device. So everything
downstream sees the new column (hd, he, hE, hTcol): kt and a_ref = kt F/c, tau, Phi_eff,
the WB pair, he_ic_balance, the BCs and the cell fill. F_r and the MLT fraction are not
changed. The solve is deterministic host code. Opacities come from the module's own
M1TableOpacities, evaluated by a device kernel per point, total kappa_t = of + os
(= kappa_R, the vet_col/vet_gd extinction). The EOS is the code's own `ThermoAt`. The
solve runs on every start, restarts included, because the BCs and WB need the same column.

Solve (grey, spherical, radiative equilibrium S = J = E, the vet_gd source in RE incl.
the vet_scatter split, which mixes to S = E when B = E):

1. Solve grid from r_b (tau_R = he_ic_rad_taub, default 3000, tau from the mesh top as in
   the pgen) to the outer fine-grid edge. chi = 0 above x1max. The points are uniform in
   the measure ds = dr/(eta r) + dtau/(dtau0 max(1,tau)^0.7), with eta =
   he_ic_rad_dlnr (2e-3) and dtau0 = he_ic_rad_dtau (0.05).
2. VET iteration:
   (a) p-z formal solve. One ray per shell tangent, plus he_ic_rad_ncore (64) core rays
       p < r_b, uniform in p^2. First-order short-characteristic weights (VcolW's),
       log-linear chi. Incoming I = 0 at the top. At r_b the diffusion intensity
       I = E_b + 3 mu F_b/c.
   (b) J, H, K by trapezoid in mu; f = K/J, g_top = H/J at the top.
   (c) The moment equation d(fE)/dr + (3f-1)E/r = -chi F_col/c (trapezoid, log-linear chi)
       is integrated inward from E_top = F_top/(c g_top).
   (d) S = E. Repeat until max|dE/E| < 1e-8 (about 8 passes).
3. Gas (`he_ic_rad_gas`):
   * `column` (default): T keeps the column's value.
   * `re`: aT^4 = E (radiative equilibrium at fixed rho). chi(rho, T) is re-evaluated in
     the loop.
4. Blend to the column: E = w E_sph1d + (1-w) E_col, with w a smoothstep LINEAR in tau
   from 1 at he_ic_rad_taua (30) to 0 at taub. A linear-in-tau blend keeps the flux
   error of the constant deep offset ~1e-4 (a log-tau blend would give ~1 %). With
   gas = re, T = (E/a)^(1/4) over the blended range too, so the deep end is the column's
   T (E_col = a T_col^4 there to 4e-6).
5. Optional HSE (`he_ic_rad_hse = true`, needs gas = re). At fixed T, rho is re-integrated
   outward from r_b on the fine grid: dP_gas/dr = -rho (g - kappa_t F/c), trapezoid,
   secant in ln rho. This is the continuum form of what the WB pair balances, Phi_eff
   with a_ref = kt F/c of the same column. Gamma = kappa_t F/(c g) is 0.2–0.8 here, < 1.
   rho is floored at the column's floor (min rho of the column above r_b, the 1e-16
   plateau) and stays there once it reaches it. The outer loop (E solve <-> HSE)
   repeats until max|d ln rho| < 1e-8.
6. `he_ic_rad_dump = <file>` (rank 0): writes the new column in the 7-column IC format
   (r rho eint F_r E T fmlt) on the union of the file radii and the solve grid. This is
   design (B): the same IC as a file, usable with he_ic_rad = file on any machine.
7. Diagnostics: R_ph old/new, iterations, cost, max |E/E_col - 1|, |T/T_col - 1| and
   |rho/rho_col - 1| over the mesh, and a few tabulated radii.

Not combinable with he_wind_ic or he_ic_rad_thin (fatal).

## 3. Tests

1. 1-D: E_sph1d against the shell means of j.31064135/be at cycles 50 and 200.
2. Reduced AG Car A on CPU (nx2 = nx3 = 32, mb 16x16, 4 ranks, full radial grid, fixed
   dt implicit_cfl = 1251.3, 30 cycles). Arms: default; sph1d; sph1d gas=re;
   sph1d re+hse; sph1d + ic_rad=fs_relax.

## 4. Results

Code: icgen-1011 **0da3369d** (sph1d), **6f8d41f1** (dump covers the fine-grid ends).
Binaries: bin/athena_he_cpu_0da3369d (md5 d8b5da5c), bin/athena_he_a100_0da3369d (md5 2362e021).
Deck: inp/agcar_A_icg.athinput (= icrad's agcar_A_icr + the he_ic_rad keys named, default file).
Reduced runs: cpu_run.sh (1 node, 4 ranks x 18 threads, ~2 s/cycle), analysis ana_ic.py, ana_dyn.py
(output ana_dyn_1.txt), ana_dt.py.

### 4.1 The solve (setup log, AG Car A, 1342 solve points + 64 core rays)

| variant | VET passes | outer | cost (CPU host) | R_ph (tau_R 2/3) |
|---|---|---|---|---|
| gas column | 12 | 1 | 0.42 s | 2.70169e13 (unchanged) |
| gas re | 18 | 1 | 0.75 s | 2.70074e13 |
| gas re + hse | 340 | 29 | 14 s | 2.70236e13 |

Resolution check (dtau 0.025, dlnr 1e-3, 128 core rays): E changes by <= 1e-3 (1.003 R_ph) and
<= 1e-4 elsewhere. g_top = H/J = 0.979 and f_top = K/J = 0.958 at the top.

Change against the file column (gas re + hse; T and rho are unchanged below 0.9 R_ph):

| r/R_ph | 0.95 | 0.98 | 0.99 | 1.0 | 1.003 | 1.006 | 1.01 | 1.02 | 1.05 | 1.2 | 1.5 | 2.0 | 2.9 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E/E_col | 0.999 | 0.963 | 0.859 | 1.032 | 0.608 | 0.465 | 0.420 | 0.397 | 0.383 | 0.398 | 0.476 | 0.601 | 0.752 |
| T/T_col | 1.000 | 0.991 | 0.963 | 1.008 | 0.880 | 0.821 | 0.798 | 0.780 | 0.754 | 0.679 | 0.592 | 0.505 | 0.416 |
| rho/rho_col | 1.000 | 1.018 | 1.004 | 0.713 | 1.71 | 1.96 | 1.56 | 1.24 | 1 (floor) | 1 | 1 | 1 | 1 |

T_col is the isothermal 9051 K above R_ph. T (re) is 7970 K at 1.003 R_ph, 6820 K at 1.05 and
3760 K at 2.9 R_ph. Without hse (gas re), rho is unchanged and E/T differ from the table only in
0.98-1.006 R_ph (E/E_col 0.937 at 0.99, 0.872 at 1.0, 0.509 at 1.003). HSE reaches the column floor
(1e-16) at 1.031 R_ph. Above it, rho stays the floor, which is not hydrostatic in any IC.

### 4.2 1-D check: cycle-0 E against the relaxed default run (j.31064135/be, full size, GPU)

E(cycle 0) / E(relaxed, cycle 50); shell means. The file IC is 2.2-2.6x off.

| r/R_ph | 0.98 | 0.99 | 1.0 | 1.005 | 1.009 | 1.02 | 1.05 | 1.1 | 1.2 | 1.5 | 2.0 | 2.5 | 2.9 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| file (default) | 1.035 | 1.081 | 1.085 | 2.177 | 2.407 | 2.532 | 2.612 | 2.615 | 2.516 | 2.112 | 1.695 | 1.462 | 1.353 |
| sph1d column | 1.003 | 1.026 | 1.097 | 1.045 | 1.009 | 1.002 | 0.999 | 0.997 | 1.000 | 1.007 | 1.014 | 1.017 | 1.020 |
| sph1d re | 1.001 | 1.011 | 0.972 | 1.010 | 1.004 | 1.000 | 0.998 | 0.997 | 1.000 | 1.006 | 1.013 | 1.016 | 1.019 |
| sph1d re+hse | 1.002 | 0.931 | 1.135 | 1.103 | 1.022 | 1.005 | 1.001 | 0.998 | 1.001 | 1.007 | 1.013 | 1.016 | 1.019 |

* Above 1.009 R_ph: within 0.3 % up to 1.2 R_ph, rising to +2 % at 2.9 R_ph. The relaxed code
  itself has E/J_fs = 1.01-1.02 there, and E/J 0.967-1.022 at cycle 30, so this is the code's
  discrete offset, not the solve's.
* 0.99-1.005 R_ph: the relaxed run is NOT a target there. Its photosphere moves (rho -20 % at R_ph
  by cycle 50; E200/E50 = 0.91-1.18). With hse, the IC rho differs from the default's there, so E
  differs as it should.

### 4.3 Reduced AG Car A, CPU (nx2 = nx3 = 32, 30 cycles, dt = 1239.996 s)

All arms rc 0, FATAL 0, NON-CONVERGED 0. Picard mean/max: default 4.23/9; sph1d re 4.43/7;
re+hse 4.53/7.

**L_top / L_top(default, cycle 0)**:

| cycle | default | sph1d column | sph1d re | sph1d re+hse | re+hse + fs_relax | default + fs_relax |
|---|---|---|---|---|---|---|
| 0 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.9999 | 1.0006 |
| 1 | **1.5283** | 1.0066 | 1.0054 | 1.0058 | 0.9998 | 1.0006 |
| 2 | 1.4419 | 1.0062 | 1.0015 | 1.0022 | 0.9998 | 1.0005 |
| 3 | 1.2490 | 1.0054 | 1.0001 | 1.0007 | 0.9998 | 1.0005 |
| 5 | 1.0407 | 1.0022 | 0.9998 | 1.0000 | 0.9998 | 1.0005 |
| 10 | 1.0014 | 0.9997 | 0.9999 | 0.9998 | 0.9998 | 1.0004 |
| 30 | 1.0004 | 1.0000 | 1.0001 | 0.9999 | 0.9999 | 1.0003 |

**E/J at the first build (cycle 0, ic_rad_diag; max over the gd shells of |<E>/<J> - 1|)**:
default 1.05; sph1d column 0.080 (top 0.921); sph1d re 0.023 (top 0.977); re+hse 0.037
(top 0.963). The code's own steady state at cycle 30 has E/J between 0.967 and 1.022.

* gas = column keeps the 9051 K gas in the thin envelope. Its thermal source in the vet_gd formal
  solution (vet_scatter: S = eps B(T_gas) + (1 - eps) J) then gives the 8 % E/J offset. Use re.
* fs_relax on top of re+hse: fs0 starts at max |E/J - 1| 0.037 (default IC: 1.05). The relax
  stops after 7 steps (default IC: 14 steps). L_IC changes by 1.3e-4 (default IC: +6.3e-4), and
  E(cycle 10)/E(cycle 0) = 1.0000 everywhere. On this IC, fs_relax has nothing left to do beyond
  the code's own discrete offset. Its cost is 9.0 s on CPU (default IC: 14.4 s).

**dt convergence** (L_top(dt, c) vs L_top(dt/2, 2c), cycles 1-30, vet_gd_rotate_every 20 at dt/2):
max (rms) relative difference: default 8.1e-2 (2.4e-2); sph1d re 1.3e-3 (2.8e-4);
re+hse 1.1e-3 (2.4e-4). For comparison, DeltaAI fs_relax full size: 2.6e-4 (8.3e-5).

**Hydrostatics (the photosphere, 0.99-1.03 R_ph)**:

| arm | max\|v_r\| c10 / c30 [km/s] | <v_r>(R_ph) c30 | max\|rho/rho0 - 1\| c10 / c30 |
|---|---|---|---|
| default | 0.99 / 2.00 | +0.79 | 0.052 / 0.164 |
| default + fs_relax | 0.96 / 1.97 | +0.79 | 0.054 / 0.166 |
| sph1d column | 0.97 / 1.98 | +0.79 | 0.052 / 0.164 |
| sph1d re | 0.92 / 1.94 | +0.80 | 0.051 / 0.160 |
| **sph1d re+hse** | **0.03 / 0.13** | **-0.02** | **0.001 / 0.019** |
| re+hse + fs_relax | 0.03 / 0.13 | -0.02 | 0.001 / 0.019 |

* The photosphere expansion of the default run (+0.8 km/s, rho -16 % by cycle 30) does not come
  from the radiation IC: fs_relax and sph1d without hse show the same. The file's atmosphere is
  "gravity only", but the code's force includes kappa_t F/c (Gamma = 0.2-0.8 here), so the
  atmosphere is out of balance. sph1d re+hse removes most of it: v 15x and the rho drift 9x lower.
* KE 0.95-1.05 R_ph at cycle 30: default 1.3e37 erg; re+hse 1.6e35 erg. KE < 0.95 R_ph is identical
  in all arms (8.2e41 erg, the convective seed).
* CAVEAT, the floor gas (r > 1.05 R_ph, rho = 1e-16 = the column floor, not hydrostatic). It falls
  in every arm. With gas = re it falls faster (max |v| 17 km/s at cycle 30, against 4.3 km/s for the
  default). kappa_R(3700-6800 K) is 1.3-8x lower than at 9051 K, so Gamma drops from 0.24 to
  0.03-0.18 at once instead of over the ~10 cycles the default gas takes to cool. KE there is 1.1e37
  erg (default 9e35); that is 1e-5 of the total. The relaxed default run's floor gas reaches
  -4 to -7 km/s by cycle 50 anyway.

**Determinism and restarts** (re+hse): two runs are bitwise (53/53 dumps and the hst). A restart from
cycle 10 is bitwise at cycle 30 (dumps; the hst differs only in the cumulative Picard counter of the
last row). The solve re-runs on the restart (14 s) and gives the same column.

**Design (B) round trip**: files/ic_agcar_A_ge_rad.txt (re+hse, md5 ca832bd8) and
files/ic_agcar_A_ge_rad_nohse.txt (re, md5 37167668) were written by he_ic_rad_dump (make_ic_rad.sh).
A run with he_ic_file = the dump and he_ic_rad = file matches inline sph1d: rho 4.5e-4 and E 4.5e-5 at
cycle 0, L_top 4e-6 over 30 cycles. NOTE: the dump covers this deck's radial grid (fine grid incl.
ghost margin), so it is tied to that radial grid.

### 4.4 Full size on GPU (Raven gpudev, 1 node x 4 A100, binary athena_he_a100_0da3369d)

Job 31067719 (60 cycles): g0 default, gh re+hse, ghfr re+hse+fs_relax, gh2 repeat, gr re. Job 31068123: gh, 200 cycles.
Job 31068190: diagnosis, 100 cycles.
- The default (g0) is bitwise the vetacc reference j.31064135/be: all hst columns over 60 cycles.
- gh vs gh2: bitwise (26 dumps + hst).
- L_top/L0 at cycle 1: default 1.5356; gh 1.0103; gr 1.0099; ghfr 0.9999. gh is within 2e-4 by cycle 5.
- E/J at the first build (max over the gd shells): default 1.061; gh 0.018; gr 0.033.
- ghfr: fs_relax took 5 relax steps, 3.3 s, L/L_IC 0.999896. E(c0)/E(c50) for gh is 1.000-1.004 above 1.02 R_ph.
- Cost: solve 14.1 s (re+hse) / 0.75 s (re) on the host. Startup to cycle 0: 40.7 s vs 26.8 s.
  s/cycle median: 0.290 default, 0.322 sph1d.
- Photosphere 0.97-1.028 R_ph, max |v| and max rho drift:

  | cycle | default | sph1d re+hse |
  |---|---|---|
  | 50 | 2.35 km/s, 33 % | 0.14 km/s, 0.9 % |
  | 100 | 2.79 km/s, 80 % | 0.20 km/s, 3.3 % |
  | 200 | 6.42 km/s, 172 % | 4.54 km/s, 63 % (convection/accretion reach it) |

  (The 14 % of the 0.99-1.03 metric at cycle 60 is the floor-gas accretion shock at 1.03 R_ph, present in both.)
- **Picard cost regression with gas = re.** Mean Picard per 20 cycles:

  | arm | 1-20 | 21-40 | 41-60 | 61-80 | 81-100 |
  |---|---|---|---|---|---|
  | default | 4.90 | 3.30 | 3.15 | 3.65 | 3.90 |
  | default + fs_relax | 3.00 | 3.35 | 3.35 | 3.40 | 3.10 |
  | sph1d column | 4.70 | 3.30 | 3.10 | 3.25 | 3.75 |
  | sph1d re | 4.60 | 4.60 | 4.65 | 6.00 | 6.00 |
  | sph1d re+hse | 4.70 | 4.60 | 4.65 | 6.00 | 6.00 |
  | re+hse, blend tau 3-30 only | 5.05 | 4.80 | 4.60 | 5.95 | 6.00 |

  The default reaches 4.0-5.0 only at cycles 100-200. The cause is the floor gas (rho = 1e-16 above 1.03 R_ph). With
  gas = re its T is 3700-6800 K from cycle 0, and its a_ref = kt F/c (from this column) is 1.3-8x lower. It then falls
  3x faster: <v_r>(2 R_ph) is -42 km/s at cycle 100 against -15 km/s for default + fs_relax, which starts from almost the
  same state but keeps the file's a_ref. The deep blend is not the cause (the tau 3-30 arm is the same).
- **Fix: he_ic_rad_gas = re_atm** (commit 3d8861c7). Radiative equilibrium and HSE apply only above the floor; the floor
  gas keeps the column T, hence the file IC's a_ref. Reduced CPU (30 cycles): Picard 3.80/7 (default 4.23/9, re+hse
  4.53/7). <v_r>(2 R_ph) at cycle 30: -4.2 km/s (default -4.1, re+hse -12.0). Photosphere 0.97-1.028 at cycle 30:
  0.09 km/s, 0.4 % (default 2.00 km/s, 16 %). L_top at cycle 1: 1.0058. dt vs dt/2: 1.2e-3 max, 2.5e-4 rms.
  E/J at the first build: 0.037. With fs_relax: 7 steps, L/L_IC 0.99987. Dump: files/ic_agcar_A_ge_rad_atm.txt
  (md5 d5c3f49d).

### 4.4b Full size, 120 cycles, DeltaAI (NOTE-2026-10-11-deltaai-icgen.md, job 3357616, 4 GH200, binary 3d8861c7)

All 6 arms rc 0, NON-CONVERGED 0. Raven's copy 31068342 was cancelled.

| arm | Picard mean (61-80 / 101-120) | s/cycle | L_top c1 | photosphere c120 (max\|v\|, rho drift) | <v_r>(2 R_ph) c120 |
|---|---|---|---|---|---|
| g0 default | 3.82 (3.65 / 4.00) | 0.207 | 1.5356 | 2.92 km/s, 96 % | -17.3 km/s |
| gha re_atm+hse | 3.53 (3.15 / 4.00) | 0.201 | 1.0104 | 0.47 km/s, 9 % | -17.5 |
| ghafr re_atm+hse+fs_relax | **3.37** (3.05 / 3.60) | 0.206 | **0.99997** | 0.46 km/s, 9 % | -17.6 |
| gra re_atm | 3.52 | 0.208 | 1.0100 | 2.84 km/s, 93 % | -17.5 |
| gh re+hse | 5.33 (6.00 / 6.00) | 0.229 | 1.0103 | 0.48 km/s, 9 % | **-51.0** |

- re_atm removes the Picard regression of re: gha is at or below the default in every 20-cycle window.
- gha vs gha2: bitwise (44/44 bin files and both hst files).
- Costs: the sph1d solve takes 5.2 s on the GH200 host; fs_relax on top of it takes 2.2 s (3 relax steps,
  L/L_IC 0.999997). On the default IC fs_relax took 4.3 s and 14 steps.
- The static photosphere comes from hse: gra (no hse) drifts like the default.

### 4.5 Recommendation (production IC)

**ghafr**: `problem/he_ic_rad = sph1d`, `he_ic_rad_gas = re_atm`, `he_ic_rad_hse = true`, plus `rad_m1/ic_rad = fs_relax`
(relax_n 40, relax_tol 1e-4).
- sph1d gives the physical start: E from the transport, gas in radiative equilibrium, and the photosphere in hydrostatic
  balance with the code's own kt F/c. This is what fs_relax alone cannot give: default + fs_relax keeps the
  out-of-balance "gravity only" atmosphere (2.9 km/s, 96 % at cycle 120).
- fs_relax then removes only the code's residual discrete offset (E/J 1.02 in the thin envelope; the +1.0 % step-1 L_top
  of gha). It takes 3 steps and 2.2 s, so L_top is flat from cycle 0 (min 0.99981). It also gives the lowest Picard.
- gha alone (no fs_relax) is an acceptable fallback, with a +1 % L_top transient over ~5 cycles.
- Optionally, the solve can be skipped at every start (5-14 s, restarts included): use the dump
  files/ic_agcar_A_ge_rad_atm.txt with he_ic_rad = file. It is valid for this deck's radial grid only.
- Do not use gas = re (+50 % Picard from the floor-gas infall) or gas = column (8 % E/J at the top).
