---
name: cs-angular-momentum
description: Why spherical polar conserves Lz to machine precision (Mignone flux-recombination, an IGNORABLE coordinate) and why the same trick CANNOT transplant to the cubed sphere; the drift is SECULAR at 6.1e-05 per rotation at nx2=64; the one route that would work is evolving rho*ell_z as a SCALAR
metadata:
  type: project
---

## Why spherical polar gets it exactly, in THIS code

`Coordinates::SrcTermsSphericalPolarHydro` (coordinates.cpp): the phi-momentum source is

```
src3 = -y_ov_rC*(F_theta A_theta + F_theta' A_theta')/V
       - factor*(F_r A_r + F_r' A_r')/V
```

with `y_ov_rC = (sin_r - sin_l)/(sin_r + sin_l)` and `factor = (r_r - r_l)/(r_r + r_l)`.
**It has NO PRESSURE TERM AT ALL**, and is built purely as a recombination of the SAME
face fluxes and areas the divergence already uses -- the Mignone (2014) form. That makes
the phi-momentum update algebraically identical to a conservative update of rho*R*v_phi,
so the ring sum telescopes to round-off.

**The enabling property is that phi is an IGNORABLE COORDINATE OF THE GRID**: the grid maps
onto itself under a phi-shift, the lever arm R = r sin(theta) is constant along phi, and
dp/dphi therefore lives entirely in the flux divergence. It buys exactly ONE axis -- even
in spherical polar Lx and Ly are not machine precision.

## Why it cannot transplant, and the one-line diagnostic

The gnomonic grid has NO continuous symmetry, only the discrete octahedral group; no
coordinate line has a constant lever arm. The signature is visible in the code:
`SrcTermsGnomonicEquiangleImpl` has

```
src3 = y_ov_rC*(ptot + rho*v2^2*sine2 - bb_2)   <-- a PRESSURE term
```

Both tangential momentum equations carry one, because neither direction is ignorable.
**That pressure source is exactly the term that cannot telescope.** No choice of
geometric source term fixes this while momentum stays the evolved variable.

## The route that WOULD work: change the variable, not the source

`ell = (r x v).zhat = x*vy - y*vx` is a genuine SCALAR with a conservation law
`d_t(rho ell) + div(rho ell v + p R phihat) = 0` valid on ANY grid. Evolve `rho*ell` in
place of one tangential momentum component, forming each face's flux with THAT FACE'S
lever arm; both neighbours then use identical (R_face, F_face, A_face) and the sum
telescopes exactly. The boundary torque vanishes identically because the cubed sphere's x1
faces are exact surfaces of constant r, so `phihat . nhat = 0` there.

Honest costs: it PRIVILEGES ONE AXIS (breaking the six-panel symmetry -- the +/-z panels
become special); the ell <-> v inversion degenerates as R -> 0, i.e. at the two poles of
the chosen axis, which on a cubed sphere sit at the CENTRES of two panels; reconstruction
then happens on ell rather than v_phi, which changes the answer near the axis. The seam
exchange actually gets SIMPLER for that component, since ell is a scalar and the current
tangential momentum is a covariant pair in the local basis.

## MEASURED: the drift is SECULAR, and here is the budget

Rigid rotation, reflecting radial BCs, nx2 = 16, relative Lz drift vs `time/nlim=0`:

```
  t= 1   -1.2684e-04
  t= 4   -5.0215e-04    x3.959 for 4x the time
  t=16   -2.0201e-03    x4.023
  t=64   -8.1147e-03    x4.017
```

**Exactly linear in t** -- a constant-rate sink, always the same sign, not a bounded
oscillation. The RATE converges at 3rd order in h (see [[cs-hydro-validation]]): 1.268e-04
per unit t at nx2=16, 1.938e-06 at nx2=64, ratio 65 for 4x resolution. With omega = 0.2
(period 31.4):

```
  nx2=16   3.98e-03 of Lz lost per rotation   -- 1% in 2.5 rotations
  nx2=64   6.09e-05 per rotation              -- 1% in 164 rotations
```

So at production resolution it is a **percent-level systematic sink over a few hundred
rotations**. Whether that matters is a physics question -- compare it against the real
torques in the problem -- but it is NOT negligible by inspection, and it is one-signed.

## RE-MEASURED 2026-09-06 (evening) on the current binaries -- rigid rotation AND the dhj runs

Rigid rotation about z (inputs/tests/cubed_sphere_rigidrot.athinput, omega 0.2, so t=1 is
0.0318 rot; ix1_bc = user, the TEST's own wall, so the fefe6a17 `reflect` mirror does not
apply), Lz at nlim=0 vs t=1, identical on the 03:00 (4990eb41) and 09:51 (build_cs) binaries:
    nx2   dLz/Lz per rot   dM/M per rot
     8      -2.96e-2         -8.3e-4
    16      -4.60e-3         -2.5e-4      (old memory: 3.98e-3)
    32      -6.86e-4         -6.6e-5
Ratios 6.4x and 6.7x per doubling (~2.7th order). Part of the Lz loss is the mass leaving
through the user radial wall (its own 2nd-order leak); the genuine source-term drift at
nx2=32 (the PRODUCTION angular resolution) is ~6e-4 per rotation of the INERTIAL-frame Lz.

dhj production runs ($S/lz.py: Lz = sum rho u R dV in the ROTATING frame, gnomonic cell
areas, StretchRPoly radii; volume sum 1.000148 x exact): the rotating-frame Lz is only
~1e-2 of the frame's Omega*I, and it evolves IDENTICALLY on cs and sp -- in units of the
frame L: rot 2 1.207e-4 vs 1.120e-4, rot 10 8.31e-3 vs 8.66e-3, rot 20 1.106e-2 vs
1.168e-2, rot 100 -4.291e-2 vs -4.295e-2 (cs vs sp). So the cs-specific drift in the dhj
runs is below ~1e-4 of the frame L (below ~1 % of the rotating-frame Lz) over 100
rotations: the gnomonic source error acts only on the rotating-frame flow, and the
rigid part is carried by the exact Cartesian Coriolis/centrifugal source. The big Lz
swings (+1e-2 to -4e-2 of the frame L) and the mass loss (-0.8 % by rot 100, same on both
grids, non-monotonic afterwards on cs) are the SETUP (open top, floors), not the grid.
cs MHD bcfix: +3.7e-3 of the frame L over rot 8-20, in line with the hydro runs.
