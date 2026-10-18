# NOTE 2026-10-18 raven: accretor 3-D theta-wall pile-up - cause found and fixed (tests only)

Follow-up to NOTE-2026-10-17-raven-accretor-3d-physics.md. USER approved (10-17 night). Branch **accretor-3d-1017 @ 89df6f64**
(pushed, NOT merged; merges cleanly with 62119e93). Design doc section **16.5**. Dir /raven/ptmp/jinma/acc3d_1017 (NOTE.md).
Setup of all arms: static, stream off, donor_irr off, 420 x 58 x 512, +-20.4 deg, phys_3d; comparisons at equal t.

## Cause: the reflecting RADIATION wall
- Not resolution: theta resolution x2 gives the same growth (x6.0 vs x5.8 at t 0.0075) in the same ~2-deg angular layer.
- A wall effect growing with the tilt: a +-30 deg wedge gives x7.9 at its wall, nothing at 23 deg.
- At the wall rows the dense gas is heated (T/T0 1.027 vs 1.0095 inside, t 0.01) and rises (v_r 3.0 vs 0.9 km/s). The tilted
  equipotentials' equilibrium (von Zeipel) flux has a poleward component; F_theta = 0 at a reflecting face makes it converge there.
- Not the hydro wall type (reflect / hse / outflow +- sponge all the same), not Coriolis/centrifugal, not the WB rows.

## Radiation theta-face options (rad_m1, both opt-in, read only when named; 2-D gates BITWISE A100 A/C and MI300A A)
| <rad_m1>/implicit_bc_x2min/max | wall-row growth (row 0, thin class, t 0.0075) | Picard | luminosity / budget |
|---|---|---|---|
| reflect (default) | x5.8 (hse), x5.7 (outflow) | 2.0 | ok |
| outflow (new: zero-gradient outgoing F, lagged) | x2.9 | **15.9** | Lmeas -13 % by t 0.0085, unbooked leak -> rejected |
| **flux (new: imposed equilibrium lateral flux)** | **x0.85 (hse), x0.79 (outflow); flat in latitude to t 0.04** | 2.0 | leak 0.9 % of L, booked in Ebnd; Lmeas unchanged |
`flux`: the pgen fills F_theta = F_col(psi) dPhi/(r dtheta)/|grad Phi| (true 3-D Roche Phi), only inside the photosphere (psi > 0;
without that cut it imposed 1.4 L sideways and broke at the walls).
s3dHfx continued to t 0.0396 on 4 A100: 0 NC, 0 retries, dt 1.39e-5 -> 1.23e-5 (thin atmosphere gas). At t 0.035 the all-phi thin
class is FLAT in latitude (rows 0/1/2/4/10/eq 15.2/15.0/15.3/15.7/16.0/15.1; reflect 43.8/32.1/22.5/19.4/17.1/15.4).

## Other measurements
- Poleward drift in the interior (~-1.75 km/s, rho0 1e-5..1e-3): the same with every wall; it is the thin atmosphere's evolution
  from the IC (the band's lateral flows), not a wall effect.
- Mass through the walls (deep envelope, net inflow): hse |dM/M| 6e-6 by t 0.012 (2e-5 by 0.04), outflow 5e-5 by 0.014.
- Energy: residual dEtot - (Ebnd + Erel - Emsk) stays constant with `flux` (the leak is booked).

## Recommendation
**theta_bc = hse + <rad_m1>/implicit_bc_x2min = implicit_bc_x2max = flux** (deck acc3d_1017/run/s3d_hse3fx.athinput) for
static/pre-impact 3-D runs; switch the hydro walls to outflow (books Mout) once the splash reaches the walls (the stream stays
inside +-12 deg).

## Open (USER)
1. After t ~0.02 the thin 3-D atmosphere (rho0 < 1e-3, ~1e-7 of the star) densifies ~x15 by t 0.035 while the 2-D band's drains
   (x0.19 at t 0.03), in every 3-D arm. The envelope agrees with the band. Not a wall effect; not yet understood.
2. `flux` uses F_col(psi) without a von Zeipel g factor; the bottom flux stays uniform (env_lin_band_factor 0.963286).
