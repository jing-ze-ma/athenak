---
name: cs-nonorthogonal-audit
description: The "treats non-orthogonal components as orthonormal" bug class on the cubed sphere -- both live instances (history.cpp KE/ME, newdt fast speed) FIXED in e415b91a; derived_variables clean, solar_convection latent
metadata:
  type: project
---

**2026-09-03.** After [[cs-mhd-c2p-floor-corrupts-ue]], swept the code for the same class:
**a quantity formed by summing squares of, or dotting, components that are NOT mutually
orthogonal on the gnomonic tangent basis.** The tangent basis has e_xi.e_eta = cos_cell
(up to -0.5 at a cube vertex) and the face normals have nhat_xi.nhat_eta = -cos_cell, so
any bare sum of squares of a tangential pair is wrong at O(cos_cell).

**Both live instances are now FIXED, committed as e415b91a** (see the status line at the
end). The diagnosis below is kept because it is what the fix is gated on.

## 1. history.cpp -- REAL, measured, diagnostic only  [FIXED e415b91a]

`LoadHydroHistoryData` / the MHD twin know about the metric for the VOLUME
(`curvi_ ? vol_(...)`, fixed in an earlier session) but NOT for the energies:

    KE:  vol*0.5*SQR(u0(IM1))/u0(IDN)  ...       <- 0.5 m_i^2/d, but m_i is COVARIANT,
                                                    so KE is 0.5 m_i v^i, with a cross term
    ME:  vol*0.25*(SQR(bx2f(j+1)) + SQR(bx2f(j))) <- face-NORMAL projections, non-orthogonal

**MEASURED on `cs_prod_hyd` (dump 00063, t = 3.84e7, max |cos_cell| = 0.4755): the history
formula overstates the true kinetic energy by 0.94%.** (Cell-summed without the volume
weight, so indicative, not exact.) The correct form is
`0.5*rho*(v1^2 + v2^2 + v3^2 + 2c v2 v3)`.

Affects REPORTED energies only, not the evolution -- but it means **cs history KE/ME are
not trustworthy to better than ~1%**, and the MHD ME is worse since it uses the
non-orthogonal face triple directly. Note the sp/Cartesian paths are unaffected (c = 0).

## 2. newdt -- the fast speed is understated by 1/sin_cell  [FIXED e415b91a]

This REFINES, and partly corrects, the "dx2 is missing sin_cell" note in
[[cs-mhd-dhj-blowup]]. `dx2_ = r*dth_xi` is the arc length along e_xi, and the true face
separation is `dx2*sin_cell`; but `w0(IVY)` is contravariant on the unit basis, and the
speed normal to that face is `sin_cell * v^xi`. **Both scale by sin_cell, so in the
ADVECTIVE ratio dx2/|v2| they cancel exactly -- the kinematic branch is right as written.**

What does NOT cancel is the wave speed: the correct constraint is
`dt2 = dx2 / (|v2| + cf/sin_cell)`, so `cf` is understated by 1/sin_cell -- up to **13.7%
at a cube vertex**, zero on the panel axes. So the timestep is too large there by at most
that, which the CFL scan already showed is NOT what kills the dhj run (halving the CFL
moved the failure 0.09%). Still a real defect; fix it on the wave-speed term, not by
rescaling dx.

## 3. Clean

`src/outputs/derived_variables.cpp` -- no bare component-square sums.
The resistive path has its own gnomonic form (`resistivity_gnomonic.cpp`, uses cos_cell).
The flux rotations, the EMF rotation and the geometric source terms were all verified
algebraically correct in [[cs-mhd-low-beta-divergent]].

## 4. Latent, only if ever run on cs

`src/pgen/solar_convection.cpp` builds the magnetic energy as
`0.5*(SQR(bcc IBX)+SQR(bcc IBY)+SQR(bcc IBZ))` in several places. That is correct ONLY
because cs stores bcc in an orthonormal frame (`GnomonicEquiangleRaiseVelMHD`); it is the
pattern to watch, not a live bug. Any NEW pgen adding magnetic energy on cs must use the
orthonormal bcc, never the raw face averages.

## 5. The fix, e415b91a

`history.cpp` KE now `0.5*m_i*v^i` = `0.5*m2*(m2 - c*m3)/(d*(1-c^2))` per angular slot;
ME now `(a_i + 0.5*c*p1*p2)/(1-c^2)`, the cross term split evenly between the two angular
slots so the PAIR sums correctly however apportioned. The radial slot is orthogonal to
both angles and is untouched. `mhd_newdt.cpp` / `hydro_newdt.cpp` use
`|v| + cf/sin_cell` on the two angular directions.

All branches gated on `use_cubed_sphere`, and `sin_cell`/`cos_cell` are dereferenced only
inside the taken branch. **VERIFIED bit-identical on a Cartesian MHD blast** (32^3,
tlim = 0.02) against the stashed pre-change code -- an A/B, not an assertion -- and the
cs_test MHD gates still pass (per-step cube-vertex order +2.00 at beta = 0.01,
dp/p / b0c^2 = 2.0076e-04).

**Consequence for old data: every cs history file written before e415b91a has KE/ME
overstated at the ~1% level.** Do not compare a pre-fix and a post-fix cs history file
digit-for-digit and read the difference as physics.
