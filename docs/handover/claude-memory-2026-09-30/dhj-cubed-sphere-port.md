---
name: dhj-cubed-sphere-port
description: WORKING -- the deep_hot_jupiter_rt pgen runs the FULL PRODUCTION PHYSICS (ck RT, tabulated EOS, point-mass gravity, eos resistivity) on the cubed sphere and tracks spherical polar. Spherical polar is BITWISE unchanged throughout. FOUR real bugs found and fixed on the way, three of them pre-existing
metadata:
  type: project
---

## What works

`src/pgen/deep_hot_jupiter_rt.cpp` now supports `mesh/use_cubed_sphere`:

- **`CSCellAngles(panel, x2v, x3v, theta, lam, phi)`** -- a new device-callable helper.  On
  the cubed sphere x2/x3 are PANEL coordinates on [-1,1], not angles, so there is no
  expression in them alone: go through `cubed_sphere::PanelToCart`, then theta = acos(z),
  lam = pi/2 - theta, **phi = atan2(-cy,-cx)**.  That last one matters -- the
  spherical-polar branch uses `phi = x3v - M_PI`, so the SUBSTELLAR point is at coordinate
  longitude pi (along -x), and atan2(-cy,-cx) reproduces both the origin and the (-pi,pi]
  interval with no wrapping.
- **15 angle sites** patched (9 in UserProblem spelled `lam = -x2v+...`, 6 in SourceFunc
  and the three RT routines spelled `lam = -theta+...` -- two different spellings, so a
  single pattern MISSES half of them).  The RT ones are where `mu0 = sin(theta)*cos(phi)`,
  the stellar beam.
- **12 `if (use_spherical_polar) x1v -= ap;`** sites now include the cubed sphere: x1 is
  the RADIUS on both grids, so radius->height must convert the same way.
- **`grav_point_mass` is allowed on the cubed sphere** (it only needs a radial x1).  The
  **stellar tide is still refused** -- it writes components in the spherical (r,theta,phi)
  basis and the gnomonic tangent pair is neither that nor orthonormal.

**SPHERICAL POLAR IS BITWISE UNCHANGED**, verified against a stashed build; the cubed-sphere
ideal-MHD resistive test is bitwise unchanged too; CPU and HIP both compile.

## TWO REAL BUGS FOUND ON THE WAY, both silent

1. **The general/tabulated EOS cache** -- see [[cs-general-eos-stale-cache]].  6.3% at the
   midlines, 15.4% at the panel corners.
2. **The cubed sphere never FILLED `x1v`, `xx1f`, `x2v`, `xx2f`, `x3v`, `xx3f`.**  They are
   ALLOCATED (the constructor guard covers both curvilinear grids) but `CoordGnomonicEquiangle`
   wrote none of them, so every consumer read ZEROS.  The pgen's hydrostatic outer BC forms
   `(x1f[i+1]-x1v[i])/(x1f[i+1]-x1f[i])` = 0/0 and turned the WHOLE DOMAIN to NaN in ONE
   cycle -- while the run reported a normal exit.  Now filled in the gnomonic setup: x1 with
   the same Mignone volume-centroid definition spherical polar uses, x2/x3 with the raw
   panel coordinates (documented as NOT angles).

## IT WORKS -- but the comparison is HYDRO+RT ONLY, see the B = 0 section below

REPRODUCED BITWISE at HEAD 6afde979 on 2026-09-01 (both dumps, both grids, against the
recorded run).  The dump-1 column is **20 cycles, not 150** -- the earlier "150 cycles"
label was wrong; the 150-cycle run exists on the CUBED SPHERE ONLY (t=1.6294e3, eint
1.4759 .. 1.0993e9, nan=0) and has no spherical-polar twin.

```
                      dump 0 (t=0)                dump 1 (20 cycles, t=2.17e2)
  spherical polar   1.729681 .. 1.099120e9      1.677038 .. 1.098951e9   nan=0
  cubed sphere      1.731221 .. 1.098744e9      1.679470 .. 1.098866e9   nan=0
```

## THE MHD PATH: was SILENTLY B = 0, now PORTED (d25d8792)

`pgen_b0` was a single `if (use_spherical_polar) { ... }` with NO else, so on the cubed
sphere it wrote nothing: `b0`/`bcc0` kept their zero-initialised values and a run with
`bbot = 3` was PURE HYDRO.  Nothing warned -- `<mhd>` was constructed and
`ohmic_resistivity = eos` ran, on a zero field.  **The history file is the gate that
catches it**: 1/2/3-ME all 0.000e+00 at t=0 on cs, against 4.745e30 / 2.372e30 on
spherical polar from the same input.  So the eint agreement above validated the
hydro + ck RT + tabulated EOS + point-mass-gravity path, and NOT MHD.

**d25d8792 ports it.**  Same dipole A_phi = 0.5*bbot*r0*sin(theta)/(r/r0)^2, built as the
DISCRETE CURL of A on the very edges mhd_ct.cpp uses (circulation over `dxedge`, divided
by `area`), so div B is round-off not truncation -- CT preserves whatever divergence it is
handed.  Follows cs_test iprob = 12.  Two cs-specific points:
  - **no angle is ever formed**: sin(theta)*phihat = (-y,x,0)/r, so
    A = 0.5*bbot*r0/(rf/r0)^2 * (-qy,qx,0) with q straight from `PanelToCart`, regular on
    the axis.  A.rhat = 0, so only the two TANGENTIAL edge components are needed -- but
    unlike spherical polar, where phihat IS x3, BOTH panel tangents see the field.
  - **bcc must use the plain average** (the non-spherical-polar branch of
    general_mhd.cpp), not the volume-centroid weights, or the magnetic energy added at
    setup and the one the first ConsToPrim subtracts disagree.
  - `xx1f` on cs is the STRETCHED radius, so use it rather than `LeftEdgeX`.

Measured (production physics, nx1 = 234 + poly stretch, 16x16 per panel):

```
  grid              div B         max |B|    total ME at t=0    eint at t=0
  cubed sphere    1.71857e-15     2.9904 G      7.155e30          1.729698
  spherical polar 3.93641e-16     2.8205 G      7.117e30          1.729681
```

Total ME agrees to 0.5%; the per-component split differs because the cs components are on
the PANEL TANGENT basis, not (r,theta,phi).  eint now agrees to FIVE digits where it
agreed to four with the field missing.  150 cycles clean with the field (nan=0, ME decays
0.06%).  Spherical polar BITWISE unchanged; cs bitwise at np=1 and np=2.

**A permanent setup gate now prints BOTH div B and max |B|** -- two numbers, because div B
alone is satisfied by a ZERO field, which is exactly what hid this ([[cs-mhd-blast]]).

## THE LAST BUG: a SHADOWED DECLARATION, pre-existing, that only the cubed sphere reached

`SourceFunc`'s coordinate block read

```cpp
Real x1v, x2v, x3v;
if (use_spherical_polar) { x1v = x1v_(m,i); ... }
else { Real x1v = CellCenterX(...); ... }   // <-- REDECLARES, shadows, discards
```

The else branch declared NEW locals that shadowed the outer ones and were thrown away at
the closing brace, so on **any** non-spherical-polar grid x1v/x2v/x3v were read
UNINITIALISED a few lines later.  Nothing had ever exercised that path; the cubed sphere
is the first, and the whole domain went NaN on the first source-term application.  Eight
other sites in the same file ASSIGN correctly -- only this one redeclares.

**HOW IT WAS FOUND, and this is the method that matters:** bisect the physics against a
control known to be clean.  `user_srcs = false` was the only clean variant (rt_ck off,
ohmic off and bbot = 0 all still failed), which pointed straight at SourceFunc.

## STRETCH PORTED, and the A/B THAT MADE IT FINDABLE

`use_grid_stretch_r` and `use_grid_stretch_r_poly` now work on the cubed sphere:
CoordGnomonicEquiangle applies the same 1-D map CoordSphericalPolar does (it is a function
of the radial coordinate ALONE and knows nothing about the angular grid, which is why it
transplants unchanged).  `use_grid_stretch_theta` stays refused -- theta is not a
coordinate of this grid.  The cs regression tests are BITWISE unchanged (they use no
stretch).

**With a VALID grid on both sides -- nx1 = 234 + the polynomial stretch, production
physics (ck RT, tabulated EOS, point-mass gravity, eos resistivity), angular resolution
reduced only:**

```
                       dump 0 (initial)            dump 1
  spherical polar   nan=0  eint 1.730 .. 1.099e9   nan=0  eint 1.677 .. 1.099e9
  cubed sphere      nan=0  eint 1.731 .. 1.099e9   ALL NaN
```

**The initial conditions agree to four digits (1.731 vs 1.730), so the port's angles,
gravity and hydrostatic setup are CORRECT.**  min eint ~1.73 is the NORMAL production
value, not an underflow -- my earlier "the upper atmosphere underflows" reading was wrong,
and it came from comparing against a different (ideal-gas) setup rather than a control.

So there IS a genuine cubed-sphere bug in the EVOLUTION, now isolated with a matched,
valid control.  Two further constraints on it:
  - the ideal-gas + two-stream-RT cs configuration ran 300 cycles CLEANLY, so it is
    specific to the tabulated EOS and/or the correlated-k path;
  - it is not the ck table lookups themselves any more (those are NaN-safe as of 99d1bbb9,
    which turned the segfault into a clean NaN).

**NEXT: bisect the physics on the STRETCHED grid** -- user_srcs off, rt_ck off (grey
instead), ohmic off, bbot=0 -- each against the spherical-polar twin, which is now known
to be clean.

## SUPERSEDED: the earlier "reduced grid" runs, all invalid on BOTH grids

Tried to run the production physics on the cubed sphere -- ck RT, tabulated EOS,
point-mass gravity, eos resistivity (bench/ck_mhd_b3's input) -- reduced to nx1 = 96
uniform and 16 per panel so it fits a login-node CPU.  It crashed.  **Then I ran the SAME
reduced setup on SPHERICAL POLAR as a control, and it fails too:**

```
  config (nx1=96 UNIFORM, no radial stretch)      spherical polar        cubed sphere
  full extent x1max 2.0556e10, point mass      dump0 ok, dump1 NaN    same IC, SEGFAULT
  short extent x1max 1.2e10,   point mass      dump0 ok, dump1 NaN    same IC, SEGFAULT
  full extent, CONSTANT g (point mass off)     dump0 ALREADY NaN      dump0 ALREADY NaN
```

**Every reduced grid I built is invalid on BOTH grids, so none of them tests the cubed
sphere at all.**  What breaks them is dropping `use_grid_stretch_r_poly`: the production
input's polynomial radial stretch is not cosmetic, it is what resolves the upper
atmosphere.  Without it the hydrostatic column underflows (min eint 1.66 against ~5e5 on a
short extent) and the tabulated EOS, which is **NaN below 71 K** ([[eos-table-dump]]),
returns NaN.  Constant g is WORSE, not better: over r/ap = 2.18 it understates H by ~2.7x
([[grav-point-mass-flag]]) and the initial condition is NaN before a step is taken.

**SO THE NEXT PIECE OF WORK IS PORTING `use_grid_stretch_r_poly` TO THE CUBED SPHERE**,
which mesh.cpp currently makes a hard FATAL.  The radial map is 1-D and independent of the
angular grid, so it is contained: apply it in CoordGnomonicEquiangle exactly where
CoordSphericalPolar applies it, and the pgen's ApplyRStretch calls then line up.  Until
that exists there is no valid production-like cubed-sphere grid to test on.

## A SEPARATE, REAL cs-SPECIFIC ROBUSTNESS BUG

When the state does go bad, the two grids fail DIFFERENTLY: spherical polar produces NaN
and keeps running (the Debug build completes 20 cycles with no bounds error), while the
cubed sphere indexes the correlated-k Planck table with **INT_MIN** and segfaults --
`Kokkos::View ERROR: out of bounds access label=("ck_pf") with indices [-2147483648,0]`.
A NaN temperature cast to an int gives INT_MIN.  Worth a clamp on that index regardless of
this port, since it turns a recoverable NaN into a crash.  Found with a **Debug build
(Kokkos bounds checking)** -- the same tool that found [[dhj-ideal-input-oob]].

## EARLIER, on the ideal-gas smoke config: NaN between cycle 300 and 400

Clean at nlim = 2, 20, 50, 100, 200, 250, 300; NaN by 400.  Smoke config: the
`deep_hot_jupiter_rt_ideal_xe` input on 6 panels, 16^3 per panel, one block per panel,
grid stretch removed (a hard FATAL on cs), polar BCs -> `panel`.

**READ THE DUMPS, NOT THE EXIT STATUS.**  Every one of these runs "completed normally";
the failure is entirely silent, exactly as [[cs-shocks-through-seams]] warns.  I reported
"it runs" once before checking the data and was wrong.

Not yet investigated: whether this is a genuine instability of this configuration (a very
coarse grid, and an input tuned for a stretched spherical-polar mesh), the cube-vertex
corner interacting with the RT column solve, or another unported assumption.  Bisect the
physics the same way it was bisected before (reflect vs user BC, srcs off, bbot=0) but at
nlim ~ 350 where it first bites.
