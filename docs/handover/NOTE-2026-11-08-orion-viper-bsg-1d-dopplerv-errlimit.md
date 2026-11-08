# NOTE 2026-11-08 Orion -> viper: 1-D BSG, Doppler velocity x error limit

**From:** Orion session (USER asked to send). **For:** viper; FYI Caltech.
**Action:** none requested. This is a result report.

## Result
With a converged implicit solve and the fixed Doppler velocity, the 1-D BSG
envelope does not contract. The shrink to ~43 Rsun seen at error_limit 1e-6 is
numerical. Most of it comes from the unconverged solve leaking energy, and a smaller part from
the over-kicked Doppler velocity v*.

Figure: `bsg-1d-orion-1108/rph1d_six.png`. Left: r_photo(tau=1) vs time. Right: E_tot(t)/E_tot(0)-1,
with an inset for 1e-8 / 1e-10 in %.

| error_limit | over-kicked v (stock) r_photo @200 d | fixed v r_photo @200 d | over-kicked dE_tot | fixed dE_tot | mean iter/solve |
|---|---|---|---|---|---|
| 1e-6  | 43.4 | 43.8 | -13.3% | -10.9% | ~1 (0-30 d), ~5 (30-200 d) |
| 1e-8  | 48.6 | 50.6 | +0.28% | +0.11% | 7.2 |
| 1e-10 | 49.5 | **51.7** (flat after ~150 d) | +0.22% | **+0.01%** | 14.5 |

All runs start at r_photo = 52.3 Rsun.
- 1e-6: both versions shrink 52 -> ~42 Rsun between ~50 and ~130 d, then sit near 43. The fixed velocity does not help here.
- Over-kicked v at converged solve: ~-2.8 Rsun / 200 d remains at dt ~ 4 s.
- 1e-8 is not fully converged, because going to 1e-10 still cuts the drift: over-kicked -3.7 -> -2.8 Rsun, fixed -1.7 -> -0.65 Rsun.
- nlimit was the default 100. Only 7 of ~8e6 solves hit the cap in each 1e-10 run, and none at 1e-6 or 1e-8.

## Terms
- **Over-kicked v (stock, vel_mode 0):** the Doppler terms of the implicit solve use v*. This is the
  velocity after the stage's hydro update (gravity and pressure kick already applied), taken before the radiation
  kick that the solve itself applies.
- **Fixed v (vel_mode 1):** a one-line `w <-> w1` swap before `GetTgasVel` in
  `nr_radiation/implicit/rad_iteration.cpp`. Stage 1 uses v^n, and stage 2 uses v^{n+1/2} after stage 1's radiation kick.

## Why an unconverged solve loses energy
In `IMRadiation::Iteration` (Jacobi), each cell's face flux uses its own new I and its
neighbours' previous-sweep I. So the two sides of a face disagree until the iteration converges. When the
loop stops (`tot_res < error_limit`, where tot_res = sum|I^{k+1}-I^k| / sum|I|), the last iterate is kept.
There is no final conservative flux pass. The gas-radiation coupling is conservative cell by cell, so the
leak is in transport. A per-solve error of order error_limit over ~4e6 steps is O(1) at 1e-6,
consistent with the -11 to -13% we see.
Caveat: E_tot is a volume integral. We have not split the drift into a boundary-luminosity mismatch versus an
internal leak with a boundary-flux ledger. The scaling with error_limit points to the solver.

## Setup
Linhao's BSG deck (col4x4) made truly 1-D:
- nx1 = 256 (his x1rat = -1 faces), nx2 = nx3 = 1;
- x2 in exactly [pi/3, 2pi/3] (1-D requires x2 symmetric about pi/2);
- npsi = 0 (1-D forbids npsi > 0), nzeta 5;
- cfl_rad 0.025, 0 to 200 d;
- Orion, 32 ranks, p.shared.

The deck for the fixed-v 1e-10 run is `bsg-1d-orion-1108/athinput.lhx_vf_e10`. The other runs differ only in
error_limit and vel_mode. Analysis scripts are `rph1d.py` (r_photo from Sigma_s+Sigma_a, tau integrated from the
outer face) and `plot1d.py`. Data on Orion: `/orion/ptmp/jinma/athpp_repro3d/run1d/lhx{,_vf}{,_e8,_e10}`.

## Implication for the paper's 3-D run (inference, not tested in 3-D)
The paper's 3-D r_photo is ~43 Rsun at 40 d, down from 1-D MESA ~51. That run uses error_limit 1e-6 with the
over-kicked v. The same two numerical effects could account for part or all of that contraction. Our 1-D
cannot say whether 3-D convection would contract the envelope by itself. The 1-D 1e-6 runs reach 43 Rsun
at ~110-120 d, about 3x slower than the paper's 3-D.

Related 3-D status on Orion: vfs23 (fixed v, cfl_rad 0.025, 1e-7, 23 nodes) has run past the ~1.1 d stall
point with no iteration caps so far.
