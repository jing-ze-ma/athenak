# NOTE h2order2: vet_gd top bath (task B) and the remaining hesdirk2 first-order terms (task A)

Worker L, 2026-10-11/12. Work dir /raven/ptmp/jinma/h2order2_1011.
Branches (fork clone, not pushed): `vgdbath-1011` (worktree wt-vgdbath), `h2order2-1011` (worktree wt-h2order2),
both from prodfinal-1011 df483770. Scratch test branch `l-test-ab` (wt-ltest) = h2order2 + vgdbath (no own commits).

## Task B: vet_gd ignored implicit_ebath_x1max (DONE, commit 9b806790 on vgdbath-1011)

- Fix: rad_m1_vetgd.cpp VetGdSweep, top shell inward rays: `iv = I_bath*ex + wu*S_top + w0*S_0` (was vacuum,
  `wu*S_top + w0*S_0`), exactly the vet_sc downward-ray form; gated by `vet_bc_bath` (same key and default as
  vet_sc: true on a fresh run, a restart file without the key keeps the vacuum top), a MARSHAK outer x1 and
  `implicit_ebath_x1max != 0`. The main field and the twin (shell-mean) field share the code path.
  Key parsed in rad_m1.cpp (GetOrAdd when vet_gd is on), member `vgd_tbath` (rad_m1.hpp).
- Production (ebath_x1max = 0): untouched branch. Reduced CPU AG Car A (480x32x32, 60 cycles, FRESH A0 deck,
  agc/agc.sh, job 31084643): every hydro_w / m1 bin dump data-bitwise vs athena_he_cpu_df483770 (headers differ
  only by the recorded `vet_bc_bath = 1` default line).
- rswt_q (tau 20 top bath, uniform equilibrium; C's decks, matched): max|E/E_b - 1| at t = 2

  | N | before (C, vacuum top) | after |
  |---|---|---|
  | 64 | 7.4e-2 | 6.1e-5 |
  | 128 | ~1e-1 | 1.5e-5 |
  | 256 | 1.7e-1 | 3.8e-6 |

  rho stays at 1 to 7e-5 / 2.5e-5 / 6.6e-6 (was +53 % at N 256). The residual is second order (top cell).
  be (FfM) identical to hesdirk2.
- C's top-bath tests are now usable (binary bin/athena_cpu_b1, runs/B, anaC = C's evaluator on runs/B):
  - rswt X (64/128/256, both dt and dx): dens 1.79, velx 1.40, E 1.31, F1 1.60 (before: -1.3 .. -1.8,
    errors 100x larger).
  - rswtv X: dens 1.31, velx 0.74, E 0.65, F1 0.67 (before: -1.4 .. -2.0).
  - rlatt X (32/64/128): dens 2.46, velz 1.88, E 1.29, F3 1.30.
  - rswt_T (hesdirk2 default start, no vimp_bc): dens 2.49/1.19, v 1.78/0.98, E 0.84/0.25, F1 1.02/0.07,
    errors 1e3-1e4x smaller than before (C: 0.2-0.4 at the 1e-3 level). The remaining T order is task A.
- NOTE: C's rsw / rswv / rlat decks also name implicit_ebath_x1max = 0.03, so those vet_gd runs had the same
  inconsistent vacuum top (smaller effect: tau 2000; rsw_q held the equilibrium to 2e-4).

## Task A: hesdirk2 time order on rsw / rswv (branch h2order2-1011, rebased on prodfinal 1fbe8621)

Commits: 50d87f3b (opt-in `implicit_rad_inertia_solve` + pgen test hook `problem/gas_hold_rho`), a90314fe
(x1 faces only). Binaries bin/athena_cpu_a4 (50d87f3b), a5 (a90314fe), athena_he_cpu_a4/a5.

### A.1 Root cause of the rswv floor: the after-solve split of implicit_rad_inertia (moveq-1011)

- Reproducer without any gas coupling (rswv_cart / mini decks): `problem/gas_hold_v` + `gas_hold_rho` (v = 0.01 sin,
  rho fixed), gas_feedback false, coupling false, kappa_p = kappa_e = 0, closure eddington, implicit_flux central,
  even 1-D transport = implicit_x1: T ladder floor / first order (kappa_s 1000: dL1 5.9 7.3 7.5 6.3 4.3 2.6 1.4e-6,
  orders -0.3 .. 0.86 over levels 1..128; kappa_s 100: 0.62 0.79 0.89; kappa_s 1e4/1e5: errors grow with refinement).
  - It is a function of dt only (N 256 levels l = N 128 levels 2l), ~ v^2 (v 0.001: 100x smaller), identical
    under time2_tableau = trbdf2 (stage order 2), under be it adds to be's own first-order term.
  - v = 0: exactly 2.00. Tolerance (tol 1e-14, lin 1e-15, EW off, lres test), predictor, start mode, closure,
    vet timing, thin switch, enthalpy (upwind/central/plm), source_form ovc, f_source, f0_init, hydro: no effect.
  - NB: H's RVt "tolerance" arm (NOTE-h2order 3.3) never applied the keys (restart decks kept 1e-12).
- A 1-D Python replica of the intended scheme (E row with upwind aE, comoving face F0 with theta, hesdirk2 with
  the code's stage algebra, py/rep.py) is exactly 2.00 for v = 0.01 and kappa 100 / 1000.
- `implicit_rad_inertia = false` restores 2.00 on the pure-transport reproducers and on rsw (all variables
  1.99-2.06); ImplicitRadInertia (after the solve) takes F0 -= theta a (E' - E^ov) on every face AFTER the E rows
  were solved without it: the stored F0 (and its ESDIRK slope) and the E update disagree by theta a dE every solve,
  a first-order splitting error ~ v (a) x v (dE advected), suppressed where theta -> 0 (stiff).
- Fix (opt-in, default bitwise): `implicit_rad_inertia_solve = true` puts dA_c = a_c (E^k_c - E^ov_c) of each Picard
  pass into F0^n through the conservative2 F0^n shift (ImplicitCs2Pass / cs2f*), restored before the write-back;
  the post-solve x1 part is skipped; while it is active the Picard stop takes no contraction estimate (the dA
  lags one pass; with conv_est the kappa 1000 held test still had a 1e-5 floor).
- x1 faces only (a90314fe): the x2/x3 faces keep the after-solve rule (no difference on any test here; G found
  the after-solve rule better on rlat, see A.4).

### A.2 Orders before / after (hesdirk2 production keys, matched, time2_start = slope + implicit_vimp_bc)

Fresh T ladders (levels 1-2, 2-4, 4-8), binaries a4/a5 (both contain the vet_gd top-bath fix; the "D" arms are the
same binary without the key):

| case | D (rad_inertia after solve) | S (implicit_rad_inertia_solve) | rad_inertia = false (diagnostic) |
|---|---|---|---|
| pure transport, held flow (rswv_cart NHvk0, kappa_s 1000) | -0.31 -0.06 0.26 (0.54 0.74 0.86 to l 128) | **2.00 2.03 1.95** | 2.00 2.00 2.00 |
| same, eddington/central/1-D (mini, kappa_s 100) | 0.62 0.79 0.89 | **2.00 2.00 2.01** | 2.00 2.00 2.00 |
| rsw_T (vet_gd, no flow) | E 1.82/1.47, v 1.78/1.29 (H FIN); 1.76/1.40/1.00 at 16 (H Hr0) | **E 1.99/1.97, F1 2.00/1.99, rho 2.03/2.05, v 2.00/1.99** | 2.0 all |
| rswv_T (vet_gd, v_r 0.01) | E 1.53/0.72, v 1.54/0.60, rho 1.20/0.98 (H FIN) | **E 1.93/1.81, F1 1.84/1.82, v 1.89/1.73**; rho 1.31/1.19 | same as S |
| rswv_cart (vet_sc twin) | E 1.63/0.74 (H FIN) | E 1.95/1.90/1.82/1.72, v 1.92/1.84/1.66/1.51 (l 1-32); rho 1.36/1.22/1.22/1.35 | same as S |
| conservative2 (G, prodfinal 1fbe8621) on rsw_T / rswv_T / rswv_cart | rsw 2.0 / rswv ~1.0 (G's old binary, no slope/vimp_bc) | with slope+vimp_bc: rsw 2.0 all; rswv E 1.90/1.77, rho 1.32/1.20 (same as S) | |

Pre-relaxed restart ladders (pre-run to t = 0.5 at dt0/8, then every level from the same state; rst.sbatch):

| case | D | S |
|---|---|---|
| rswv_cartr (vet_sc) | E 0.26 0.14 0.34, rho 0.77 0.25 0.36 (the H "RV floor") | **E 2.06 2.05 1.90, F1 2.01 2.01 1.71, rho 1.99 1.96 1.92, v 2.02 2.00 1.69** (l 1-16) |
| rsw_Tr (vet_gd) | E 1.01 0.35, rho 1.20 0.33 | **E 2.05 2.06, F1 2.02 1.97, rho 2.03 2.04, v 2.02 1.96** |
| rswv_Tr (vet_gd) | E 0.33 0.19, rho 0.72 0.22 | E 1.78 1.79, F1 1.84 1.18, rho 2.09 1.98, v 1.92 1.43 |
| rswv_Tr, vet_gd_rotate_every 0 | | **E 1.86 1.74, F1 1.85 1.72, rho 2.08 1.95, v 1.86 1.73** |

- The fresh-ladder rho ~1.2 in rswv is start-up (IC off the stiff manifold with the slope probe of frac 1): it
  is gone in the restart ladders (rho 1.9-2.1) and improves with time2_start_frac 1e-3 (rho 1.46 1.65 1.79).
- rswv errors at the finest pairs are 2-6e-7 of the pulse; the residual dip of F1/v at 4-8 with rotation 10 is
  the direction-set rotation (step-count dependent), 1.7 with rotation 0.

### A.3 Regression (a4 = a5 on these; S vs D, same binary)

- Reduced CPU AG Car A (480x32x32, 60 cycles, FRESH A0): **key absent: data and file bitwise vs
  athena_he_cpu_1fbe8621** (a4 and a5, every bin). With the key: changes <= 7e-8 (F1 6.7e-7) of max; Picard mean
  3.02 -> 4.18, CPU +11-14 % (the stop without contraction estimate; allowing it from pass 3 brings the held-flow
  floor back: a6 test, reverted).
- rw_t10a T 1-16: D = S to 3 digits (E 2.00 1.99 1.91, F1 1.98 1.84 2.02/2.00, rho/v 2.0); C 32-256 identical.
- rw_p10t100 (tol 1e-14/1e-15): D = S, all 2.0. rw_p10t1: D = S (E 1.79/1.73, F1 1.51/1.04; H's open item is
  NOT rad inertia).
- eq3t (hesdirk2 + slope): D E 1.136750 P_tot 3.002068 (= H); S E 1.136763 P_tot 3.002062.
- rlat_T fresh: D 1.94/2.16; S 1.80-1.82 / 1.66-1.82 (errors 1.4x). Same as G's conservative2 finding.

### A.4 Open

- rlat (lateral wave) prefers the after-solve rule on the fresh ladder (D 2.0 vs S 1.8, G: any in-solve or no
  inertia gives 1.8). Its pre-relaxed restart ladder (rlat_Tr) has a floor in BOTH arms (D 0.1-0.6, S 0.2-0.5,
  errors 2e-5 of the pulse in v): probably the restart does not carry the vet_gd lagged band/rotation state or
  the predictor increment ("cold start" warning) - not resolved.
- rw_p10t1 E/F1 1.0-1.7: not this mechanism (S = D).
- Cost of S: +1.2 Picard passes. A version with dA in the E-row stencil (implicit, no Picard lag) would avoid
  the extra passes; not done.
- Default: S changes AG Car by <= 1e-6 relative; whether to default it (fresh runs) is for the user; rlat
  argues for keeping the x2/x3 faces on the after-solve rule (done) and checking rlat's restart artefact first.

## Files

- ana/: h.py (H's harness, W = this dir), ev.sh, prof2.py, bq.py (bath equilibrium), cord.py (C-mode orders);
  anaC/ (C's evaluator on runs/B); py/ (rep.py Python replica, chk.py stage-algebra checker, rep3/rep4 vs code).
- base/: H's decks + mini.athinput (eddington/central).
- jobs_*.txt / *.sh, logs/, runs/<case>/<arm>_<lev>, agc/ (reduced AG Car, runs/j.<job>).
- bin/: athena_cpu_b1 (vgdbath), a1 (+gas_hold_rho), a3/a4 (50d87f3b), a5 (a90314fe), he_cpu_{b1,a4,a5};
  a2dbg (row dump instrumentation, not committed), ab1/ab2 (test merges).
