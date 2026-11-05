# Athena++ BSG inner-wall collapse: 1-D implicit DO toy (11-04)

## Summary
- **The toy reproduces the Athena++ collapse.** The toy is the converged Athena++ DO discretisation coupled to vl2/PLM/HLLC hydro, with his IC, BCs and opacities. With his BC (arm A0), the wall is drained:
  - inner-face L is -100..-125 L, into the bath;
  - dE_tot is -3.9 / -9.8 / -13.5e45 erg at 1 / 2 / 3 d. Athena++ "fix" (converged sweep=2, his BC tc4, ledger "boundary") gives -3.2 / -11.5 / -18.9e45.
  - The wall cells form the hot/cold pattern: cell 0 hot (1.2e6 K); cell 1 cold and dense (T 1.06e6 -> 9.4e5 -> 8.2e5 -> 6.3e5 K, rho x2).

  With a closed wall (mirror radiation ghost + reflecting hydro, arm A2, like Athena++ Bd/Be):
  - no energy is lost (dE ~ 0);
  - the cell-scale wall instability appears: cells 0-1 hot and rarefied, cell 2 cold and dense (T 4e5 K, rho x3 at 2 d);
  - the run dies at ~2.5 d. Athena++ Bd1 converged dies at 2.7 d with the same pattern.
- **Mechanism: the Doppler (O(v/c)) terms of the radiation solve take the WRONG VELOCITY.**
  - Athena++ hands the implicit radiation solve the stage velocity v* = u(IM1)/rho of the hydro-updated state (`GetTgasVel`: "using the quantities from the partially updated u"). That is v^n + dt * a_hydro, BEFORE the radiation force of the same step.
  - In the thick base the gas is held up partly by radiation. Each step the hydro update accelerates the wall gas inward by dt * a_rad, about 2.3 km/s per 212 s step. The radiation momentum source then gives it back: v* = -2.95 km/s, but v^{n+1} = -0.63 km/s at 1 d (`diag_wall.py`).
  - The solve sets the lab flux to F = F_diff + (4/3) v* E/c. That flux is 100% Doppler: cell 0 carries -270 L, and F_diff is only a few L.
  - The spurious part is dF = -(4/3) E dt a_rad, inward. Its ratio to the true flux is Lambda = (4/3) kappa E_rad dt / c. That is about 66 at 1.2e6 K (kappa 0.45 cm^2/g, dt 212 s).
  - The near-central DO face flux (taucell tau-factor HLL, f -> 1/tau) puts half of the cell-centred flux on each face:
    - at the bath face, half of it leaves the star;
    - at a closed face it piles energy into the wall cells and drives the cell-scale hot/cold instability.
  - Proof: the drain is proportional to dt (below).
- **Recommended mimic:** none. It is a first-order operator-splitting error (it scales with dt), not physics. AthenaK does not have it (see the AthenaK section), so it is correct that AthenaK does not reproduce the collapse.
  - If an Athena++-like run is wanted anyway, the compact recipe is an opt-in key that feeds the v-dependent radiation terms the post-hydro, pre-radiative-kick velocity v*.
  - An effective wall-sink term is the cheaper alternative: L_sink ~ -(2/3) 4 pi r_w^2 E_w |v*_w|, with |v*_w - v_w^n| ~ dt kappa F/c.
  - Neither a DO boundary layer nor DO angles are needed.

## The toy (code in this directory)
- `toy.py`: radiation and gas energy only, with the gas frozen.
- `toyhd.py`: the radiation-hydro version. It follows Athena++ (053bc73; its radiation directory is unchanged since 05-14, so it equals e47fabb7):
  - 1-D spherical, nx1 384, uniform 20-80 Rsun, his units;
  - angle set: 2*nmu Gauss-Legendre (nmu 4: 8 angles, w/2);
  - per-angle HLL face flux with the tau_scheme-1 taufactor f_r/f_l, tau = taucell (dx sigma_L + dx sigma_R);
  - the implicit part uses the velocity (mu c - adv), plus the explicit adv_vel I flux (donor in stage 1, PLM in stage 2, GetTaufactorAdv);
  - the comoving source exactly as AbsorptionScattering (tran_coef, wmu_cm, Planck/Rosseland split);
  - lab-frame billing of dE and dF to the gas (GetHydroSourceTerms, AddSourceTerms);
  - vl2 structure: stage 1 is hydro donor + radiation solve over dt/2 from I^n; stage 2 is hydro PLM + radiation solve over dt from I^n;
  - point mass GM, cfl 0.3 (dt ~ 200-212 s);
  - his IC: profile linear in r, the +-2.5 % isobaric rho seed, I+- = (E/sum w+ +- F/sum mu w+)/2;
  - his ghost bath, including ComToLab with the ghost velocity;
  - his HLLC zero-contact hydro ghost (iterated on the PLM face state) and his outer BCs;
  - Ma2026 Rosseland/Planck tables, kappa_es = 0.34.
- **The solve:** full Newton on all (I, T) unknowns. That is the converged fixed point, the same as the viper sweep=2 "fix" runs, not the unconverged default Jacobi. Not implemented: angle_flag 1, the pre-solve RHS clip, the Jacobi iteration.
- **Running it:** all switches are command-line options of `toyhd.py`. `batch.sh lists/<file>` runs a list in a held allocation (1 srun step, nice 19).
- **Analysis files:**
  - `summ.py`: the tables;
  - `plot.py`: the figures;
  - `diag_wall.py`: the per-term wall decomposition from `state/*.npz`.
- **Run data:** logs and json in `runs/`.
- **Cost:** about 1 min CPU per star-day.

## Results (3 d unless noted; dE in 1e45 erg; wall L in units of the profile L, > 0 into the star)

**Base and the Doppler arms.** All arms except A2 and B1 use his BC.

| arm | change | dE 1 d | dE 3 d | T1 at 3 d | wall L | verdict |
|---|---|---|---|---|---|---|
| A0 | Athena++-like base | -3.9 | -13.5 | 6.3e5 | -125 | collapse (as Athena++) |
| A2 | closed wall (mirror + reflect) | -0.06 | dies 2.5 d | (T2 4e5 at 2 d) | 0 | wall instability (as Bd1) |
| A8 | Doppler (tran_coef) off | +3.2 | +4.2 | 1.12e6 | +5..+15 | stable (10 d: B2) |
| B1 | closed + Doppler off | -0.07 | -0.08 | 1.07e6 | 0 | stable |
| B4 | Doppler off in cells 0-2 only | +3.2 | +3.8 | 1.15e6 | +0..+18 | no drain |
| B6 | Doppler on in cells 0-19 only | -3.6 | -14.4 | 6.1e5 | -126 | drains = the wall layer |

**Velocity handed to the radiation solve.**

| arm | change | dE 1 d | dE 3 d | T1 at 3 d | wall L | verdict |
|---|---|---|---|---|---|---|
| V1 / V5 | the solve uses v^n instead of v* | -0.3 | +0.7 (10 d: +2.6) | 1.04e6 (10 d 1.07e6) | +8..+15 | **cured**, 10 d stable |
| V3 | closed + v^n | -0.06 | -0.09 | 1.07e6 | 0 | **cured** |
| V2 / V4 | v* + explicit radiative-kick prediction | NaN at once | | | | explicit kick unstable (chi c dt huge) |

**Null tests: these do not matter.**

| arm | change | dE 1 d | dE 3 d | T1 at 3 d | wall L | verdict |
|---|---|---|---|---|---|---|
| A1 | no rho noise | -3.6 | -12.3 | 6.4e5 | -131 | same |
| A3 | no ghost ComToLab boost | -6.0 | -14.8 | 9.7e5 | -105 | same |
| A5 | no radiation advection flux | -3.8 | -16.1 | 4.3e5 | -158 | same / worse |
| A6 | taucell 1 | -3.9 | -13.5 | 6.4e5 | -124 | identical |
| A9 | Rosseland in the exchange | -3.9 | -13.5 | | | identical |
| A10 | nmu 8 | -3.9 | | | | identical |
| A11 | bath + reflecting hydro | -2.4 | -9.8 | 9.1e5 | -73 | same |
| B7 | nx1 768 | -2.7 | | | -88 | same |
| A7 | pure upwind DO | | | | +5e4 | upwind DO is unusable in the thick base |

**dt scan** (his BC, dE at 1 d; `fig_dt.png`). The drain is proportional to dt:

| dt [s] | 424 | 212 | 106 | 53 | 26.5 |
|---|---|---|---|---|---|
| dE at 1 d | -8.2 | -4.1 | -1.9 | -0.96 | -0.57 |

The closed-wall arms with dt/2 and dt/8 (C4, C5) stay stable to 1 d; A2 is already unstable by then.

**Frozen gas** (`toy.py`, radiation and gas energy only, rho and v frozen):
- the wall is stable for 1 d;
- the bath feeds +1 L;
- this holds with and without noise and with the mirror wall.

=> Hydro motion is required, through the Doppler terms. The DO transport, taucell, angle count and gas-radiation exchange alone are inert.

**Wall decomposition** (A0 at 1 d, one step, `diag_wall.py`):

| cell | stage velocity v* | v^{n+1} | cell flux | (4/3) v* E | transport (diffusion) part | momentum the radiation gives back per step |
|---|---|---|---|---|---|---|
| 0 | -2.95 km/s | -0.63 km/s | -266 L | -270 L | +3 L | +2.3 km/s |
| 1 | -2.55 km/s | | -85 L | -88 L | | |
| 2 | -1.52 km/s | | -81 L | -79 L | | |

## Coordinator's hypothesis: radiation enthalpy through a zero-mass-flux face
**(a) Measured face-0 flux** (A0, 1 d): Lw = -125 L. It splits into:
- **explicit advective part: -78 L.** v_face = mean(ghost +0.65, stage cell -2.95) = -1.15 km/s, so (4/3) E v_face A = -104 L.
- **implicit DO part: -47 L.** It contains half of the cell-0 Doppler flux, (2/3) v* E = -134 L, plus the tau-factor dissipation and the ghost flux.

So yes, radiation enthalpy leaves through a face that carries no mass. But the velocity in it is v*, not the physical cell velocity.

**(b) Velocity set to zero only at the wall face:**

| arm | change | dE 1 d | dE 3 d | wall L |
|---|---|---|---|---|
| A4 | explicit advection v_face = 0 | -3.3 | -14.6 | -126 |
| W1 | v_face = 0 for advection and the ghost boost | -5.9 | -15.2 | -100 |
| W2 | W1 + v^n in the solve | -4.4 | -8.5 | -30..-60 |

=> **Zeroing the face velocity does NOT cure the drain.** The leak moves into the DO part: the near-central face flux averages the cell-centred Doppler flux (4/3) v E onto the wall face whatever the face velocity is. W2 is worse than V1: there the face advection with the ghost velocity had been partly cancelling the residual (2/3) v^n E leak.

The cure is the velocity value. With v^n (V1, V5) the cell Doppler flux is small, and the cell-0 velocity relaxes to -0.1..-0.4 km/s. With v* it is (4/3) E dt a_rad, roughly 65x the diffusive flux.

The BC inconsistency is real, but it is the outlet, not the source. The closed wall (A2) has no outlet and still blows up, from the same v*.

**(c) AthenaK, same face** (`wt_lhx1104`, `rad_m1_implicit.cpp` ~14630-14700):
- The fin4 he_bc_linhao interior-face formula does NOT zero the wall advective flux. It uses vf = face average of the ghost (HLLC zero-contact) and cell velocities in `ch*dt*vf*g0f` and the enthalpy stencil.
- The linhao-prod-1104 halfrange option adds (2/3) v_g E_g + (2/3) v_is E_is, also not zeroed.
- AthenaK does not drain (K1: +2.5 L) because of the VELOCITY its v-terms see:
  - the lagged stage-start v^n (time2 "lag"), the same as toy V1;
  - or, under implicit_vimp, the solve's own iterate including the radiative kick (time2_vstage), i.e. the self-consistent end-of-stage velocity.
- It never uses the post-hydro, pre-kick v*. Caveat: implicit_vimp is silently off on 1-D meshes, with non-hesdirk2 schemes and with non-bicgstab solvers; those fall back to lag = v^n. That is still correct here.

## Mimic recommendation
1. Do not mimic. The collapse is an O(dt) splitting artifact of Athena++'s v* in the Doppler terms. It vanishes as dt -> 0 and with v^n or a self-consistent velocity. Linhao's wall collapse and envelope deflation driven by the bath drain are numerical.
2. If an Athena++-reproduction arm is wanted:
   - an opt-in rad_m1 key that evaluates the v-dependent terms (comoving correction, F = F0 + (4/3) v E/c, the wall face) at v* = v^n + dt a_hydro, i.e. the gas velocity after the hydro stage without the radiation kick;
   - no DO layer;
   - expected: a wall drain of about -(2/3) 4 pi r^2 E_w dt kappa F_w / c per unit L (about 100 L at dt 212 s) and the hot/cold wall pattern on a closed wall.
3. Do not chase taucell, angle count, bath shape, the ghost boost, noise or the clip: all are inert in the toy.

## Caveats
- The toy has angle_flag 0 only, and no Jacobi or clip (it matches the converged "fix" runs, not default Jacobi).
- r(rho=1e-10) also dips to 46 Rsun at 2 d and rises to 62-66 Rsun at 3 d in the cured V1 arm. That envelope oscillation is not the wall drain. It is absent only when Doppler is off. Possibly af0 + atmosphere, not studied.
- The toy's A0 crashes (NaN) at about 2.7-3 d in some realisations. Athena++ af0 arms die at about 5 d.
- Figures:
  - `fig_main.png`: base, closed wall, Doppler off, v^n;
  - `fig_nulls.png`: the null tests;
  - `fig_wallv.png`: the wall-velocity tests and the dt extremes;
  - `fig_dt.png`: the dt scan.
