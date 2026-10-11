# NOTE: rad_m1 implicit moving-medium fixes (moveq-1011)

Branch `moveq-1011`. Code commit 7a9af8bc. Final binary `athena_cpu_mq3`, built from that commit.
Work dir: `/raven/ptmp/jinma/moveq_1011` (runs/, ana/, ladder/, fig/).
Baseline binary: `athena_cpu_base`, built from a61c1e75.

## 1. Analytic check (IDORT 4.3.1 case 3, moving thermal equilibrium)

Setup: rho = 1, gamma = 5/3, P = rho T, a_r = 1, c = 10, v0 = 3, E0 = T0 = 1, uniform and periodic.
The radiation thermalises in the comoving frame. Total lab momentum rho v + F/c^2 = 3 and total
energy are conserved. Source: `ana/eqana.py`.

| model                              | v       | T       | E (lab) | F/c^2   |
|------------------------------------|---------|---------|---------|---------|
| relativistic                       | 2.95677 | 1.00017 | 1.12850 | 0.04323 |
| code O(beta^2) mixed frame         | 2.95537 | 1.00018 | 1.13263 | 0.04463 |
| O(v/c), E = E0                     | 2.95686 | 1.02279 | 1.09432 | 0.04314 |
| gas not decelerated (old behaviour)| 3.00000 | 0.97685 | 1.03473 | -       |

The old code reached E = 1.0325 and v = 3.0011, so the gas was not decelerated, and the total
momentum P_tot went to 3.0425. Radiation momentum was created out of nothing.

## 2. Root causes

- (a) **Cold start of the face flux.** F0 = 0 at the start, so the cell flux is F = F0 + aE = aE.
  The radiation therefore starts with comoving-advected momentum that no gas momentum paid for.
- (b) **Missing -dA/dt term.** The F0 face equation lacks the -dA/dt term, where A = aE. As E
  changes, the advective part of F changes with no momentum exchange.
- (c) **Momentum term not given to the gas.** The chat*dt*v*g0 term is removed from F but never
  added to the gas momentum.
- (d) **gas_feedback = false.** The work v.dm is skipped, so a diffusing pulse in a moving medium
  is advected at (4/3) v by the enthalpy flux alone, not at v.

## 3. Fix (7a9af8bc, `src/rad_m1/rad_m1_implicit.cpp`, `rad_m1_impl_wb_kernel.hpp`, `rad_m1.hpp`)

New keys:

- `rad_m1/implicit_f0_init` = auto | zero | lab. `auto` starts F0 from the lab flux consistent
  with F, which fixes (a).
- `rad_m1/implicit_rad_inertia` (default true). Adds the -dA/dt radiation-inertia term, including
  the a.W/c^2 work momentum, and gives the removed chat*dt*v*g0 momentum to the gas. This fixes
  (b) and (c).
- `rad_m1/implicit_nofb_work` (default true). With gas_feedback = false, it deposits the work
  v.dm at the fixed v, which fixes (d).
  - `ImplicitNofbWorkRow` puts this work into the E row, Picard-lagged, instead of adding it
    after the solve.
  - A test build with the row switched off (mq4) gave **bitwise-identical** ddn outputs, and
    mq2 = mq3 bitwise on eq3t. The row is therefore inert in every test run here.
  - It is kept as committed, but it is unvalidated and could be removed.

## 4. Results: before (base, a61c1e75) and after (mq3)

### 4.1 Moving equilibrium, eq3t, t = 1

Target values (code frame): E 1.1326, v 2.9554, P_tot 3.

| run                    | before E / v / P_tot       | after E / v / P_tot        |
|------------------------|----------------------------|----------------------------|
| cen be                 | 1.0325 / 3.0011 / 3.0425   | 1.1333 / 2.9551 / 2.99984  |
| cen hesdirk2           | 1.0325 / 3.0011 / 3.0425   | 1.1333 / 2.9551 / 2.99984  |
| tJe be                 | 1.0325 / 3.0011 / 3.0425   | 1.1333 / 2.9551 / 2.99984  |
| tJe h2+vimp+vstage     | 1.0325 / 3.0011 / 3.0425   | 1.1333 / 2.9551 / 2.99985  |
| hr be                  | 1.0286 / 3.0030 / 3.0445   | 1.1261 / 2.9585 / 3.0033   |
| cen ovc                | 1.0696 / 3.0023 / 3.0453   | 1.1771 / 2.9534 / 2.99993  |
| cen nofb (v fixed)     | 1.1367 / 3 / 3.0456        | 1.1365 / 3 / 3.0455        |

- The intermediate build mq1 lacked the a.W/c^2 work momentum. Its cen be run gave
  P_tot 3.0046 and E 1.1221.
- `hr` still leaks 1e-3 of momentum.
- In the nofb run the gas has infinite inertia, so P_tot is not conserved by construction.

Figure: `/raven/ptmp/jinma/moveq_1011/fig/moveq_equil.png` (made by `ana/fig_eq.py`).

### 4.2 Dynamic diffusion, ddn (gas_feedback = false) and ddf (feedback), t = 16

The columns are the peak relative to the analytic value, the centroid minus v t, and the
variance-slope ratio D_eff/D.

| run          | before peak / xc-vt / D_eff/D | after peak / xc-vt / D_eff/D |
|--------------|-------------------------------|------------------------------|
| ddn cen be   | 0.885 / **+5.333** / 1.290    | 0.891 / -0.0000 / 1.272      |
| ddn cen h2   | 0.998 / **+5.333** / 1.001    | 0.984 / -0.0000 / 1.033      |
| ddf cen be   | 0.875 / +0.0001 / 1.233       | unchanged                    |
| ddf cen h2v  | 0.970 / -0.0001 / 0.984       | unchanged                    |
| ddf hvy be   | 0.784 / -0.0001 / 1.635       | unchanged                    |

- ddn: the pulse used to drift at (4/3) v (+5.33 at t = 16). It now sits at v t.
- ddn be with cfl 0.1 gives a peak of 0.970 at t = 4, against 0.902 at the default cfl. The low
  be peak is therefore BE time error.
- ddf: unchanged at 4 digits. The low be peak (13 %) is BE error plus the real gas response.
  hesdirk2+vimp+vstage gives 0.970.

### 4.3 Static and Marshak data arrays, base vs mq3

Compared with `ana/bincmp.py` (.bin data arrays) and on .tab data lines.

| case                                     | result                                           |
|------------------------------------------|--------------------------------------------------|
| nltek1/4/8 (cen, tJe), sdifk (cen, tJe)  | bitwise                                          |
| marsh c1 (cen, hr, tJe)                  | final dump differs at about 1e-10 relative (E, F, eint); v ~ 1e-6 there |
| eq1 (cen, hr, tJe), v = 0                | E and eint bitwise; only noise-level v ~ 1e-17 and F ~ 1e-12 differ at roundoff |
| eq2 (cen, hr, tJe)                       | gas develops v ~ 1e-3; F differs at 1e-4 relative, E and eint at 1e-7 relative |

All of these differences come from the new inertia and work terms acting on the small velocities
present in these runs.

### 4.4 Order ladders, rw_t10a and rw_t1000a, 120 rows

Arms: {cen, tJe} with be, and {cen, tJe}Pvs with hesdirk2. Studies X, C and T.
Tables: `ladder/lad_run_old.txt` and `ladder/lad_run_new.txt`.

- 15 of 120 rows differ, all in the third significant digit.
- Every order is preserved. hesdirk2 stays at about 2 in X, C and T.

For example, rw_t10a cenPvs hesdirk2 C e: L1 1.30e-03 becomes 1.17e-03 at N = 128, and the
orders go from 2.15/1.74/2.03 to 2.16/1.88/1.87.

### 4.5 Linear radiation-modified acoustic waves (c = 10, Ma+2025 4.3.2)

Tables:

- before: `/raven/ptmp/jinma/littests_1011/ana/rw_table.txt`
- after: `ana/rw_table_mq3.txt`, made with `ana/rw/run_radwave.py` and `ana/rwtab.py`

Results:

- At P = 0.01, every row is identical except the 6th digit of w_re.
- At P = 10, L1 changes by at most 1 % and the orders are unchanged.
- For P = 10, tau = 100, cen, the damping rate goes from -0.43669 to -0.43704. The reference is
  -0.43666.
- The `hr` failures at N >= 128 (P = 10, tau = 0.01 and tau = 1) were present before and are
  unchanged.

## 5. Open items

- `hr` keeps a 1e-3 momentum leak in eq3t.
- `ImplicitNofbWorkRow` is bitwise inert in all tests; validate it or remove it.
- The BE D_eff/D of about 1.27 in ddn/ddf be is not produced by the work split: the row on and
  off give the same result. It is most likely BE time error of the M1 flux relaxation.
  hesdirk2 gives 1.03.
