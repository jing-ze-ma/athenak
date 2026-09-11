---
name: sp-geometric-source-residual
description: MEASURED 2026-09-05 -- the spherical-polar geometric source leaves a spurious force on a uniform field that is 14x larger in the POLAR ROWS than in the interior (1 part in 630 of |T| at the production resolution, the cube-vertex level); first order in my normalisation, exactly linear in B^2. sp_test pgen + inputs/tests/spherical_polar_uniform_field.athinput
metadata:
  type: project
---

**Instrument:** src/pgen/sp_test.cpp (build_sp_test, `cmake -D PROBLEM=sp_test`), uniform gas
at rest + uniform Cartesian field b0 zhat set as FACE AVERAGES (discrete div B exactly 0),
nlim=1, bin dumps every step, F = rho dv/dt from two dumps (scratchpad sp_src/measure.py).
Polar MHD needs >= 2 MeshBlocks in phi ("divisible by two in x3" check).

    |F| / ((b0^2/2)/r), one RK2 step, raw HLLD
    grid                 interior mean   interior max   polar rows (mean = max)
    sp 16 x 32           1.14e-2         2.47e-2        1.26e-1
    sp 32 x 64           4.65e-3         1.24e-2        6.41e-2
    sp 64 x 128 (prod)   2.33e-3         6.12e-3        3.23e-2
    cs 16/panel wb off   7.41e-3         5.68e-2 (vertex)  edge rows 1.50e-2
    cs 32/panel wb off   2.79e-3         3.71e-2           5.42e-3
    cs 64/panel wb off   1.60e-3         1.81e-2           2.36e-3
    cs 64/panel wb ON    6.19e-4         1.80e-2           1.85e-3
    b0 = 0.1 / 1 / 10 on sp 32x64: identical to 4 digits -> pure geometry, linear in B^2.

In this normalisation a 2nd-order cancellation error eps (relative to |T|) appears as
eps*N, so eps = F*dtheta: sp interior ~1e-4, sp POLAR ROWS 1.6e-3 (1 in 630) at 64x128;
cs interior 1.4e-4, cs cube vertex 2e-3 (1 in 500, = the old memory number). So sp's
polar rows carry the SAME class and SIZE of defect as the cube vertex, uniform along the
whole row. Low-beta criterion: the spurious force exceeds the physical pressure-gradient
scale when beta < ~eps; the dhj top reaches beta 1e-4, i.e. 16x past it in the polar rows.

**Cause (by construction of SrcTermsSphericalPolarMHD):** the off-diagonal stresses already
use the face FLUXES, and the diagonal coefficients are exact area differences, but the
anisotropic diagonal stresses (rho v_th^2, rho v_ph^2, B_th^2, B_ph^2 in m_ii_h / m_pp) are
CELL-CENTRE values -- and in the polar row the cell-centre trig differs most from the face
averages (the pole face has zero area). The face-sum (well-balanced) form of
[[cs-wb-source-cached]] would remove it; on sp the face basis is just (rhat, thhat, phhat)
at the face angles. NOT built yet. Whether it changes any sp answer is unmeasured; the sp
polar blow-up was the EMF ([[sp-polar-field-blowup]]), not this.

## BUILT AND MEASURED (2026-09-06 early): the face-sum source does NOT touch the polar rows

`<mhd>/sp_wellbalanced_src` (commit after 037b1554; SrcTermsCurvilinearWB, BuildWBGeometry sp
branch), DEFAULT OFF. cs bitwise unchanged. sp: interior 2.33e-3 -> 1.56e-3 at 64x128 (the
r-component gets WORSE, theta better), polar rows 3.23e-2 -> 3.21e-2. Then the diagnosis:
the polar residual is independent of reconstruction (dc/plm/ppm4 identical to 4 digits), of
the energy definition (sp_test problem/econsistent), and of the source form. It is the pole
itself: B_theta is linear through zero across the polar cell, so the outer theta-face's
normal field is 3/2 the cell's centroid-interpolated bcc (bcc on sp IS interpolated to x2v,
ideal_mhd.cpp:93 -- my "arithmetic mean" hypothesis was WRONG, the code already does this),
the geometric terms scale as cot(theta) ~ 1/dtheta there, and their O(dtheta^2) relative
cancellation error leaves 0.65 dtheta in units of (B^2/2)/r. Second order in the local
terms; the sp analogue of the cube vertex's 1-in-500. No source-term rewrite removes it.
**RETRACTED:** my claim that this residual is "the sp equivalent of what the WB source removes
on cs". It is not. [[measure-impact-before-claiming]]

## 2026-09-06: CARTESIAN polar-row momentum update BUILT (34882bd7, `sp_cart_polar_momentum`, default OFF)
Face fluxes rotated to Cartesian with each face's basis + exact vector-area correction
T~.(Avec - A n); source masked in the polar rows; pressure closure exact (null test clean).
Result on sp_test (32x64, polar-row means, units (B^2/2)/r):
    axial field (bdir=2):  r 4.0e-3 -> 5.1e-3, theta 6.40e-2 -> 6.37e-2  (NO change)
    x-field    (bdir=0):   r 4.18e-1 -> 4.7e-3 (90x!), theta 3.45e-2 -> 1.76e-2, phi 0.11 same
So the radial spurious force of a TRANSVERSE field at the pole (the bigger one, 0.42) is
removed; the axial-field theta residual (0.65 dtheta) survives EVERYTHING tried: solver
(hlld/hlle/llf), reconstruction (dc/plm/ppm4 identical to 4 digits -- which contradicts the
j-flux reconstructing bl/br from bcc, UNEXPLAINED), energy definition, source form,
Cartesian update. It is fixed by the face-normal fields alone. The cs vertex oracle
([[cs-vertex-oracle-halo-innocent]]) says the cs vertex residual is likewise not the halo.
Next discriminator if resumed: instrument uflx at the polar row's outer theta-face and
compare with the exact T.n (sp_test knows the exact field).
