---
name: sp-polar-row-reconstruction
description: THE polar-row force residual on sp is EXPLAINED (2026-09-06) -- the polar cell is reconstructed FLAT because the mirror ghost equals the cell for every even component, so the limiter clips it as an extremum; fixed 7.7x by an extremum-preserving quadratic (6410be99, mesh/polar_quadratic_recon, default OFF, shock-safe with the 3-triple curvature bound). Also -- sp ALWAYS uses GridPiecewiseLinearX2, <mhd>/reconstruct is IGNORED on sp
metadata:
  type: project
---

**Instrument:** sp_test `problem/flux_probe=1` (+ `problem/user_srcs=true`): prints the code's
momentum flux at the polar cell's faces next to the exact T.n.  Verdict on the axial field:
80 % of the polar-row residual is the OUTER theta-face flux; HLLD has no momentum
dissipation for v=0, so the flux is the L/R mean of p_tot -- the R state (from cell j+1)
is 0.99534 ~ exact 0.99519, the L state is 0.997592 = the polar cell's MEAN: the polar cell
is FLAT. Cause: the mirror ghost across the pole holds the same value as the cell for every
component EVEN across the pole (B_r of an axial field, B_phi of a transverse one, rho, p);
minmod(0, x) = 0; every TVD limiter is first order at an extremum and the pole IS an
extremum of every even profile. This is why dc/plm/ppm4/wenoz were identical -- AND because
on sp the flux kernels always call GridPiecewiseLinearX2 (non-uniform PLM) regardless of
`reconstruct` (mhd_fluxes.cpp:393, hydro_fluxes.cpp:329). **`<mhd>/reconstruct` does
nothing on spherical polar.** Worth knowing for every sp run ever made.

**Fix (6410be99):** in the polar rows only, when the cell is an extremum, the quadratic
through (mirror ghost, cell, neighbour) at their centroids, curvature bounded by the left
AND right triples (all same sign, C = 1.25, Colella-Sekora), face values clamped to the
three cells' range; monotone -> old PLM. `<mesh>/polar_quadratic_recon`, DEFAULT OFF.

    sp_test axial, |F| in (B^2/2)/r:   polar row      neighbour row    worst row
    nx2=32   old -> new                6.29e-2->8.2e-3  1.06e-2->2.84e-2  2.2x better
    nx2=64                             3.16e-2->4.1e-3  5.41e-3->1.45e-2  2.2x better
    (the neighbour worsens: the shared face's old error had cancelled part of its own)
    x-field phi residual 0.11: UNCHANGED (B_phi is even but uniform in theta -- different origin, open)
    hydro blast on the axis (inputs/tests/spherical_polar_pole_blast.athinput): same cycle
    count, mass/E to 6 digits, no NaN, diffs confined to the polar rows. The UNLIMITED
    quadratic went NaN; the centre+right-only bound gave a 3.7x lower min density -- the
    three-triple bound is REQUIRED.

Together with [[sp-geometric-source-residual]] (Cartesian update: transverse-field radial
residual 90x) the sp pole is now: axial theta 7.7x, transverse radial 90x, transverse phi
open. NOT yet tried on the dhj run; both flags default off. Next: a dhj sp arm with
polar_quadratic_recon (+ sp_cart_polar_momentum) vs sp_mhd_diss.

**dhj ARM LAUNCHED (2026-09-06 01:20): bench/sp_mhd_pole2** = sp_mhd_diss's input + BOTH flags
(polar_quadratic_recon, sp_cart_polar_momentum), from scratch, HIP binary of 6410be99
(md5 0e591d78...), jobs 11428062-63 chained 24 h apu1. Compare against sp_mhd_diss at matched
rot: 1-ME (axis energy), horizontal KE, dt, and the polar-row B_r/B_phi profile from restarts.
