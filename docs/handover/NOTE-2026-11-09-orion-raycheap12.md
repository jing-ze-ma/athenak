# NOTE 2026-11-09 ORION -> viper: raycheap12 -- time order of the M1 scheme on clean 1-D benchmarks

Answer to `TASK-2026-11-09-orion-raycheap12.md`. No code changes.

**Interpretation (short)**
- **hesdirk2 IS 2nd order on the clean problems.** EQUIL non-stiff p 1.95 / 1.97 / 1.98 / 1.99 against a python ODE
  reference; diffusion pulse (transport only, scattering, tau_cell 10) p 1.95 / 1.98 / 1.99. be gives 1.00 / 0.99 on both
  (the reference behaves). So the scheme itself has no 1st-order term; the p ~ 1 of the cs test star (raycheap3-11) comes
  from that setup, not from the hesdirk2 stepping of exchange or transport.
- **Stiff EQUIL (k dt0 = 10):** no order reduction visible. Measured at t = 0.2 = 20 relaxation times (the TASK's
  "t = 1 relaxation time" cannot be reached with k dt0 = 10), all schemes damp the transient faster than the
  ODE does at large k dt (p 7 for hesdirk2 / trbdf2, 2.5-2.8 for be), and hesdirk2 / trbdf2 settle at p 2.3-2.6 at the
  finest pairs (k dt 1.25 -> 0.6). trbdf2 equals hesdirk2 to 3 digits at the three finer levels.
- **Marshak (stiff exchange + transport + a cold-medium front):** hesdirk2 p 1.76 / 1.50 / 1.18 from t = 0 and
  2.04 / 1.75 / 1.58 from a common restart at t = 10 (initial layer removed); trbdf2 2.02 / 1.26 / 0.99 and 1.88 / 1.50 / 2.31;
  be 1.02 / 1.01. So Marshak shows a mild, dt-dependent reduction for both 2nd-order tableaus; **trbdf2 does not fix it**
  (not better than hesdirk2). Likely reason (not verified): the discrete front in the cold medium (E ~ 1e-18 ahead,
  T^4) moves ~0.002 per unit time, i.e. one cell (dx = 1/128) per ~4 time units, which our dt range (10 -> 0.625)
  straddles; the hesdirk2 runs also show Newton fallbacks (142-278 per run) and, from t = 0, a few positivity clips.
- **vet_col:** refused on the 1-D mesh (`closure = vet_col: needs <rad_m1>/transport = implicit on a MULTI-D mesh`).
  On a quasi-1-D 2-D mesh (nx2 = 4) Marshak accepts it: hesdirk2 + vet_col p 1.68 / 1.61 / 1.55 (t = 0) and
  1.77 / 1.71 / 1.62 (restart), vs the same mesh with m1 1.52 / 1.49 / 1.10 and 1.69 / 1.77 / 1.31. vet_col does not lower
  the order here, and **vet_source = gas vs relaxed_esc makes no difference to the order** (d within 4-8 %, same p).
  The pulse cannot take vet_col: periodic x1 is refused (`closure = vet_col: a periodic x1 boundary has no top`), and on a
  multi-D mesh the pgen pulse is a 2-D Gaussian in (x, y) (pulseQ runs discarded).
- Refused / skipped: vet_col on 1-D (all); pulse + vet_col (above). Nothing else refused; all 150 + 72 runs rc = 0,
  0 NOT ADMISSIBLE, 0 DT_RETRY.

## Binary
/orion/ptmp/jinma/raycheap11_1109/athena_addb642d0 (head addb642d0, PROBLEM none; the test-only keys are not used).

## Deck changes (from inputs/tests/rad_m1_{equil,thick_pulse,marshak}.athinput; generator gen.py)
- all: `<time>` integrator rk2 (already), nlim -1, tlim (below), **dt_fixed = dt** (constant dt ladder, halving),
  cfl_number 1e6, **dt_retry = false** (otherwise dt_retry lifts the dt_fixed cap: "DT_RETRY recovered ... cap off" and
  the runs took fewer, larger steps; a first submission was discarded for this); `<rad_m1>` transport = implicit,
  subcycle = false, cfl_rad 1e6, implicit_tol 1e-12 (was 1e-8), implicit_maxit 200 (was 30), time_scheme hesdirk2 | be,
  time2_tableau trbdf2 (inserted; absent keys go into deck copies), closure / vet_source for the vet_col sets.
  Outputs: bin m1 + hydro_u at t = 0 and tlim. Checked: every run ends at tlim with cycles = tlim/dt (EQUIL / pulse have one
  extra round-off step of ~1e-15 at L0, L3, L4).
- EQUIL: c = 1, arad = 1, x1 in [0,1] (8 cells, periodic), gas_rho 1, gas_eint 1.5 (T = 1), e_rad 0.1;
  kappa_p = kappa_e = kappa_f = 0.2725 (non-stiff: k = kappa (1 + 4 a T^3/c_v) ~ 1, k dt0 = 0.1; t = 1) or 27.25 (stiff:
  k dt0 = 10; t = 0.2); dt0 = 0.1. Reference: scipy Radau (rtol 1e-13) of dE/dt = c rho kappa (a T^4 - E),
  T = (gamma-1)(E_tot - E)/rho (agrees with the finest level to 7e-6 / 4e-8).
- Pulse: kappa_s 1280 (tau_cell 10), tlim 2.4 (the pgen profile is amp exp(-r^2/w^2), so sigma0^2 = w^2/2 = 0.00125 and
  the variance doubles at t = sigma0^2/(2D) = 2.4 with D = 1/(3*1280)); dt0 0.24.
- Marshak: tlim 100 (front at x ~ 0.36, inside the slab), dt0 10; **implicit_ebath_x1min = 1e-10** (= e_bath: with
  transport implicit the x1min boundary is the implicit Marshak BC, and its bath was 0 in the deck, so no radiation
  entered; a first Marshak set was discarded for this). Quasi-1-D: nx2 = 4 (mesh and meshblock), implicit_allow_multid true.
- Common restart (Marshak): base run per scheme x closure with dt = 10 to t = 10 (+ rst), then every level restarted from
  that file to t = 100 with time/dt_fixed. (A base with a small dt does not work: Mesh::NewTimeStep limits dt growth to
  2x per step, so the coarse levels ramped up for several steps; with the base at the coarsest dt every level's dt is
  exact: 9 / 18 / 36 / 72 / 144 steps.)

## Results (d = |E_k - E_k+1|_1/|E_finest|_1 between successive levels, p = log2 ratio)

### EQUIL non-stiff (k dt0 = 0.1, t = 1 = 1 relaxation time; dt 0.1 ... 0.00625)

| scheme | closure | d (self, 4 pairs) | p (3) | error vs ODE ref (5 levels) / p |
|---|---|---|---|---|
| hesdirk2 | deck | 1.211e-03 3.161e-04 8.084e-05 2.047e-05 | 1.94 1.97 1.98 | 1.636e-03 4.243e-04 1.082e-04 2.734e-05 6.876e-06 |
| hesdirk2+trbdf2 | deck | 1.234e-03 3.222e-04 8.244e-05 2.092e-05 | 1.94 1.97 1.98 | 1.666e-03 4.326e-04 1.104e-04 2.791e-05 6.990e-06 |
| be | deck | 8.092e-03 4.167e-03 2.115e-03 1.066e-03 | 0.96 0.98 0.99 | 1.649e-02 8.410e-03 4.247e-03 2.134e-03 1.070e-03 |

p vs ODE ref: hesdirk2 1.95 1.97 1.98 1.99; hesdirk2+trbdf2 1.95 1.97 1.98 2.00; be 0.97 0.99 0.99 1.00

### EQUIL stiff (k dt0 = 10, t = 0.2 = 20 relaxation times; dt 0.1 ... 0.00625)

| scheme | closure | d (self, 4 pairs) | p (3) | error vs ODE ref (5 levels) |
|---|---|---|---|---|
| hesdirk2 | deck | 2.163e-02 1.449e-04 9.588e-07 1.475e-07 | 7.22 7.24 2.70 | 2.178e-02 1.461e-04 1.145e-06 1.862e-07 3.872e-08 |
| hesdirk2+trbdf2 | deck | 2.063e-02 1.327e-04 9.588e-07 1.475e-07 | 7.28 7.11 2.70 | 2.077e-02 1.338e-04 1.145e-06 1.862e-07 3.872e-08 |
| be | deck | 1.104e-02 1.960e-03 2.756e-04 4.116e-05 | 2.49 2.83 2.74 | 1.333e-02 2.289e-03 3.285e-04 5.299e-05 1.184e-05 |

p vs ODE ref: hesdirk2 7.22 7.00 2.62 2.27; hesdirk2+trbdf2 7.28 6.87 2.62 2.27; be 2.54 2.80 2.63 2.16

### Diffusion pulse (kappa_s 1280, tau_cell 10, t = 2.4 = variance doubled; dt 0.24 ... 0.015; 1-D periodic, m1)

| scheme | closure | d (self, 4 pairs) | p (3) | L1 error vs exact diffusion Gaussian (5 levels) |
|---|---|---|---|---|
| hesdirk2 | deck | 5.573e-04 1.438e-04 3.654e-05 9.216e-06 | 1.95 1.98 1.99 | 2.179e-03 1.622e-03 1.478e-03 1.442e-03 1.433e-03 |
| hesdirk2+trbdf2 | deck | 5.573e-04 1.438e-04 3.654e-05 9.216e-06 | 1.95 1.98 1.99 | 2.179e-03 1.622e-03 1.478e-03 1.442e-03 1.433e-03 |
| be | deck | 4.063e-03 2.090e-03 1.060e-03 5.339e-04 | 0.96 0.98 0.99 | 9.713e-03 5.650e-03 3.560e-03 2.500e-03 1.966e-03 |

### Marshak, from t = 0 (t = 100; dt 10 ... 0.625), 1-D, m1

| scheme | closure | d (self, 4 pairs) | p (3) |
|---|---|---|---|
| hesdirk2 | deck | 3.681e-04 1.086e-04 3.848e-05 1.697e-05 | 1.76 1.50 1.18 |
| hesdirk2+trbdf2 | deck | 2.070e-04 5.115e-05 2.129e-05 1.075e-05 | 2.02 1.26 0.99 |
| be | deck | 1.130e-02 5.550e-03 2.741e-03 1.362e-03 | 1.03 1.02 1.01 |

### Marshak, from t = 0, quasi-1-D 2-D mesh (nx2 = 4): m1 vs vet_col (vcR = relaxed_esc default, vcG = vet_source gas)

| scheme | closure | d (self, 4 pairs) | p (3) |
|---|---|---|---|
| hesdirk2 | deck | 4.286e-04 1.498e-04 5.342e-05 2.500e-05 | 1.52 1.49 1.10 |
| hesdirk2 | vcR | 6.925e-04 2.159e-04 7.095e-05 2.428e-05 | 1.68 1.61 1.55 |
| hesdirk2 | vcG | 7.185e-04 2.270e-04 7.460e-05 2.554e-05 | 1.66 1.61 1.55 |
| hesdirk2+trbdf2 | deck | 2.445e-04 7.368e-05 4.138e-05 2.840e-05 | 1.73 0.83 0.54 |
| hesdirk2+trbdf2 | vcR | 2.918e-04 1.004e-04 3.828e-05 1.532e-05 | 1.54 1.39 1.32 |
| hesdirk2+trbdf2 | vcG | 3.169e-04 1.132e-04 4.146e-05 1.658e-05 | 1.49 1.45 1.32 |
| be | deck | 1.135e-02 5.573e-03 2.762e-03 1.367e-03 | 1.03 1.01 1.01 |
| be | vcR | 1.142e-02 5.581e-03 2.754e-03 1.363e-03 | 1.03 1.02 1.01 |
| be | vcG | 1.127e-02 5.501e-03 2.713e-03 1.342e-03 | 1.03 1.02 1.02 |

### Marshak, COMMON RESTART at t = 10 (all levels from one rst; then t = 100), 1-D m1

| scheme | closure | d (self, 4 pairs) | p (3) |
|---|---|---|---|
| hesdirk2 | deck | 3.504e-04 8.545e-05 2.536e-05 8.495e-06 | 2.04 1.75 1.58 |
| hesdirk2+trbdf2 | deck | 2.547e-04 6.919e-05 2.441e-05 4.931e-06 | 1.88 1.50 2.31 |
| be | deck | 1.123e-02 5.531e-03 2.730e-03 1.354e-03 | 1.02 1.02 1.01 |

### Marshak, COMMON RESTART at t = 10, quasi-1-D: m1 vs vet_col

| scheme | closure | d (self, 4 pairs) | p (3) |
|---|---|---|---|
| hesdirk2 | deck | 3.909e-04 1.213e-04 3.552e-05 1.437e-05 | 1.69 1.77 1.31 |
| hesdirk2 | vcR | 6.421e-04 1.880e-04 5.744e-05 1.867e-05 | 1.77 1.71 1.62 |
| hesdirk2 | vcG | 6.677e-04 1.986e-04 6.085e-05 1.976e-05 | 1.75 1.71 1.62 |
| hesdirk2+trbdf2 | deck | 2.876e-04 8.926e-05 2.746e-05 1.730e-05 | 1.69 1.70 0.67 |
| hesdirk2+trbdf2 | vcR | 2.966e-04 8.250e-05 2.621e-05 9.982e-06 | 1.85 1.65 1.39 |
| hesdirk2+trbdf2 | vcG | 3.147e-04 8.983e-05 2.857e-05 1.108e-05 | 1.81 1.65 1.37 |
| be | deck | 1.127e-02 5.546e-03 2.747e-03 1.356e-03 | 1.02 1.01 1.02 |
| be | vcR | 1.136e-02 5.562e-03 2.746e-03 1.360e-03 | 1.03 1.02 1.01 |
| be | vcG | 1.120e-02 5.482e-03 2.706e-03 1.339e-03 | 1.03 1.02 1.02 |

## Jobs and files
Jobs 228087-228102 (1-2 exclusive nodes, 1-rank runs, all concurrent per job; whole study ~40 min). /orion/ptmp/jinma/raycheap12_1109: gen.py, decks/, decksR/, run*.sbatch, restart.sbatch, ana.py, ana.out, anaR.out, runs*/*/run.log. bins and rst deleted.
