---
name: red-giant-deep-onset-is-rcb-pileup
description: "The 10x-MLT deep motion in the red giant is a STARTUP TRANSIENT: L switched on into a non-convecting envelope piles heat into the first cells above the RCB, sup grows linearly from t=0 (3.5e-4 per 1e6 s), convection ignites at t~7e5 from 200x the MLT-equilibrium driving and overshoots. Not the wall, not numerics, not waves."
metadata:
  type: project
---

Measured 2026-09-09 on `grey_prod` (61 dumps, t = 0..3e6), scripts `onset.py` and the
pile-up check in the scratchpad.  This RETRACTS the "checkerboard numerical instability"
reading of earlier the same night and the "wall cannot shed the flux" reading before it.

## The sequence, in the data

| t [s] | sup(i=5) | sup(i=6) | v_rms(i=10) | v/v_mlt surface (i=210-270) |
|---|---|---|---|---|
| 0 | +1.9e-6 | +9.9e-7 | 7.5e3 (seed) | 0.06-0.23 (seed) |
| 3e5 | -5.5e-5 | +7.3e-5 | 8.7e3 | 0.03-0.09 (seed decaying) |
| 6e5 | -1.6e-4 | +1.7e-4 | 4.0e3 | 0.15-0.47 growing |
| 9e5 | -2.9e-4 | +3.4e-4 | 3.0e4 growing | 0.28-0.82 |
| 1.5e6 | -1.1e-3 | +7.9e-4 | 4.3e4 | 0.54-0.78 SATURATED |
| 3.0e6 | | | 2.9e5 = 10 v_mlt | 0.5-0.7 |

- A superadiabatic/subadiabatic DIPOLE forms at the first cells above the RCB
  (i = 5/6, r = 0.044-0.047 R) and grows LINEARLY from t = 0, before any velocity grows.
  That is the injected luminosity piling up where radiation stops carrying it and
  nothing else has started to.
- Deep velocity grows only from t ~ 7e5, once sup(6) ~ 2e-4; e-fold ~4e5 s, which is the
  convective growth rate for sup ~ 5e-4 (sigma^2 ~ g sup/H_p).  Starting from ~200x the
  MLT-equilibrium superadiabaticity it overshoots to v ~ sqrt(200) ~ 10-20 v_mlt.  Then
  it overmixes (sup(9-11) = -3e-3 at 3e6), spreads outward as a front, and the KE peaks
  at ~3e45 near t ~ 8e6 in `prod_topre` and declines -- a relaxation transient.
- The SURFACE convects first and correctly: 0.5-0.8 v_mlt by t ~ 1e6, saturated after,
  from 0.97 to 1.06 R.  It carries only 0.33 L because the photosphere has 0.08
  horizontal cells per H_p ([[red-giant-flux-deficit-is-spinup]] has the fluxes).
- <v_r>/v_rms deep stays < 0.01: no radial pulsation.  Wall cells stay within 0.2 % of
  the initial column in p and T (v_rms(i=0) = 600 cm/s at 3e6): the wall is quiet.
- The seed is one smooth mode per MeshBlock: v1 ~ sin(3*2pi*x), v2 ~ sin(2*2pi*xi)
  (period 4 cells in an 8-cell block), v3 ~ sin(2pi*eta).  The period-4 anti-correlation
  at lag 2 in j that looked like a checkerboard IS THE SEED, and the growing mode rides
  on it until t ~ 2e6 and then cascades to large scales (corr(j,j+1) = +0.8 by 3e6).
- Total energy in the hst is conserved to 8e-6 over 3e6 s.  The `face budget` print's
  "inner gain = -7000 L" is NOT a real leak (that would be -4e-4); do not chase it.

## What to do about it

Do not switch L on impulsively into an envelope with no convective flux.  The remedy
under test is `problem/mlt_alpha = 3.0` (= mlt_alpha_ic), the pgen's own subgrid MLT
flux: at t = 0 the column is then in flux equilibrium F_rad + F_mlt = L everywhere, no
pile-up, and the resolved flow grows from the real (small) superadiabaticity.  It also
carries the surface flux the grid cannot resolve.  Runs R0_base/R1_mlt/R2_llf/R3_nowb
(8x8 per panel, 1.1 R domain, tlim 3e6) in /orion/ptmp/jinma/Athenak/red_giant/.
Related: [[rt-transparent-cell-cancellation]] (the separate ambient-medium failure).
