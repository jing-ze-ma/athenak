---
name: cs-mhd-instability-characterized
description: The cs+MHD blow-up is a low-beta instability that is NOT gnomonic-specific -- spherical polar does it too, at the same beta, with cs worse by ~10x. Refinement DELAYS onset on both grids. Two retractions: the vertex order deficit (shrinking bins) and "the rate rises with resolution" (window-dependent fits)
metadata:
  type: project
---

**2026-09-03. Commits a6401406 (cs diagnostics) and df2edc8e (the sp twin, `sp_lowbeta`).**
Supersedes the "skew cube-vertex operator" diagnosis in [[cs-mhd-minimal-reproducer]] and
[[cs-mhd-dhj-blowup]].

## THE HEADLINE: it is not gnomonic-specific

`src/pgen/sp_lowbeta.cpp` + `inputs/tests/spherical_polar_mhd_lowbeta.athinput` is the
SPHERICAL-POLAR twin of cs_test iprob = 8: same shell r = 1..1.5 at nx1 = 16, same uniform
medium at rest (d0 = 1, p0 = 0.1), same tilted uniform Cartesian field from A = (1/2)BxR,
same gamma/plm/hlld/rk2/cfl, starting dt equal to three digits. theta is a WEDGE with exact
analytic ghosts (the polar axis has no cs analogue); phi is the full 2 pi, periodic.

**Spherical polar is unstable too, at the same beta.** So the feedback is NOT in the
gnomonic non-orthogonal machinery; it is in the structure sp and cs SHARE -- a curvilinear
geometric source built from the cell-centred `bcc` against fluxes built from reconstructed
face states, which cannot cancel for a force-free field the way the uniform-pressure state
does. The cubed sphere is worse, but by a factor, not in kind.

Fit-free gate (KE per unit volume, so the cs full shell and the sp wedge compare; the
exact solution is at rest, so KE IS the error):

    beta = 0.02      time to reach KE/V = 1e-6        KE/V at t = 0.8
      n=16    cs 0.059   sp 0.084                cs 5.2e-5   sp 5.7e-6
      n=32    cs 0.329   sp 0.849                cs 5.1e-6   sp 9.2e-7
      n=64    cs never(*) sp 2.769               cs 4.4e-7   sp 1.1e-7
    (*) the cs n=64 record only runs to t = 1.0.
    Reaching KE/V = 1e-4: cs 1.50 (n16), 2.20 (n32);  sp 2.00 (n16), NEVER at n >= 32.

cs crosses any level 1.4-2.6x sooner and gets an order of magnitude further. The RADIAL
MONOPOLE control (`iprob = 2` on sp, `iprob = 1` + `b0r` on cs) is clean on BOTH grids --
KE rings at ~4e-6 and decays -- so on both grids the instability needs TANGENTIAL field.

## RETRACTION 1: "the growth rate rises with resolution"

**Not supported.** Every rate in this thread (mine, and the 13.8-vs-17.3 pair in
[[cs-mhd-minimal-reproducer]]) is a least-squares slope of ln(KE) over a window that
DIFFERED FROM RUN TO RUN, because each run was fitted over the last half of whatever record
it had before dying. The growth is not a single clean exponential, so that slope moves by
2-3x with the window: cs n=64 at beta 0.02 reads +3.64, +1.02 or +1.53 depending on it.

Fit-free, refinement DELAYS onset at every point of the cs ladder:

    beta     n=16 t(1e-6)   n=32 t(1e-6)     KE/V at t=3 (n=16 -> n=32)
    0.200    never          never            5.1e-7 -> 4.2e-8
    0.102    0.616          never            2.0e-6 -> 1.9e-7
    0.050    0.245          1.284            1.1e-5 -> 2.0e-6
    0.020    0.059          0.329            (both dying)
    0.006    0.013          0.101            (both dying)

That is what a fixed instability fed by a SMALLER TRUNCATION SEED looks like. It also kills
the companion claim that the threshold beta rises with resolution. What survives: growth
from rest sets in below beta ~ 0.1 and is violent by beta ~ 0.02, on both grids.

**The lesson, and it is general:** never read a growth rate off an exponential fit whose
window is set by when the run happened to die. Use a fit-free crossing time. Cf.
[[validate-the-instrument]], [[measure-impact-before-claiming]].

## RETRACTION 2: the cube-vertex order deficit

The region bins in `CSTestConvErrors` were a fixed TWO CELLS wide, so they shrink
physically under refinement -- the trap already recorded in [[cs-seam-order-limiter]].
`problem/conv_nband` is now an input. Per-step L1(v_t), b0c = 1:

    region        fixed 2 CELLS      fixed PHYSICAL width      nband = 1
    interior       1.65 / 1.81         1.86 / 1.96            1.73 / 1.91
    seam           1.43 / 1.57         1.81 / 1.92            1.46 / 1.31
    CUBE VERTEX    1.33 / 1.17         1.79 / 1.91            0.94 / 1.15

One picture fits all of it: an error density with an integrable 1/distance profile at the
vertex. **The operator is second order in L1 everywhere.** A perfect halo
(`exact_panel_ghosts=1`) reproduces every one of these to 3 digits.

## The RHS split, the measurement that was queued

`<mhd>/cs_diag_no_coordsrc` and `cs_diag_no_divf` drop one half each of the momentum RHS.
(dt is RADIAL-bound in this test, so it is the same dt at every nx.) Each half is O(1) and
resolution-independent, and their SUM is second order:

    region      -div F     source     their sum (nx=64)   cancellation
    interior    4.18e-3    4.19e-3    2.31e-6             1 : 1800
    seam        4.81e-3    4.86e-3    6.19e-6             1 :  780
    vertex      7.52e-3    7.67e-3    1.49e-5             1 :  500

`cs_diag_no_magsrc` is NOT a usable probe: the inconsistency it introduces is O(1) and the
state goes ballistic (KE ~ t^2 at 1e-2 within a few steps), drowning the effect.

## Other eliminations, all by direct A/B

* **The halo is not it, for stability either.** `exact_panel_ghosts=1` blows up
  indistinguishably from the baseline (stratified reproducer, b0c = 3.162, nx = 16).
* **Reconstruction is not it.** ppm4 changes the seed, not the behaviour.
* **It is NOT vertex-localised.** Per-region L1(v_t) at t = 0.05/0.20/0.40 grows by the
  same factor in interior, seam and vertex; the vertex amplitude is only ~1.5x the
  interior's. The mode fills the panel interior.
* **Gravity and stratification are not needed.** The unstratified iprob = 8 does it all.

## THE FORCE BUDGET -- why "the real field is not force-free" does not rescue it

Objection: the tests use a force-free field, real fields carry current.  Answer, measured.

The scheme's SPURIOUS acceleration on a uniform field (b0c = 3.162, B^2 = 10, rho = 1,
per-step at nlim = 1, a = L1(v_t)/dt):

    nx      interior      seam       CUBE VERTEX
    16      2.38e-2     4.63e-2      7.74e-2
    32      6.59e-3     1.28e-2      2.19e-2
    64      2.00e-3     3.62e-3      6.05e-3
    order    1.85/1.72                          <- second order, as it should be

Now the PHYSICAL force a static equilibrium can carry: |J x B| = |grad p| <~ p/L =
beta*B^2/(2L), which is 10*beta here (L = 0.5).  Setting the two equal gives the beta below
which the scheme's error EXCEEDS the entire physical force budget:

    nx=16   beta = 2.4e-3 (interior)   7.7e-3 (vertex)
    nx=32          6.6e-4              2.2e-3
    nx=64          2.0e-4              6.1e-4

The stratified reproducer reaches beta = 1.3e-4 at the top of the domain.  At nx = 64 that
is BELOW the interior crossover: no equilibrium there is representable at all.

**Direct confirmation.**  cs_test iprob = 11 is a genuinely current-carrying equilibrium
(uniform `bunif` plus azimuthal `b0c`*(-y,x,0), whose curl is 2*b0c*zhat, with the pressure
balancing a real J x B).  At bunif = 3.162, p0 = 0.2, beta = 0.02..0.06, the current
b0c = 0.01 gives |J x B| = 0.063 -- SMALLER than the nx = 16 vertex spurious force above.
Against its own b0c = 0 control:

    b0c=0 (force-free)   n=16 t(KE/V=1e-6) 0.059, dies t=2.06;  n=32 0.360
    b0c=0.01 (current)   n=16               0.057, dies t=2.06;  n=32 0.328

**Adding a real current changes nothing.**  That is the point: at low beta an equilibrium
field is force-free to within O(beta) BY NECESSITY -- the Lorentz force it exerts can only
be as large as the pressure gradient -- so the physics is a small residual of two large
cancelling terms, and our error in that cancellation is currently the same size as the
residual.  Note the parameter trap: on iprob = 11 `b0c` is the AZIMUTHAL amplitude and
`bunif` the uniform one, and `bazi` is not read at all.  Getting that wrong makes p0
hugely negative and the run is garbage -- check the t = 0 magnetic energy against
0.5*B^2*V = 49.74 before believing any iprob = 11 result.

## THE WELL-BALANCED SOURCE: BUILT, AND IT DOES NOT FIX IT

**c3992145, `<mhd>/cs_wellbalanced_src`, default OFF.** Builds the geometric source as the
face sum it is meant to cancel, S_i = (1/V) sum_f A_f (T_cell . nhat_f).e_i(x_f), with
CARTESIAN vectors -- the cell's v and B assembled on the cell's own triad, every face's
normal and basis from the panel map at THAT face's (xi,eta). For a constant Cartesian
stress the divergence theorem makes this identically equal to the flux divergence, so the
two cancel to round-off whenever the face states equal the cell state, on any grid, with no
metric identity required of the areas. It is a strict generalisation: for the ISOTROPIC
part it reduces algebraically to the old coefficients, so the uniform-pressure balance is
carried over, not re-derived. Verified against spherical polar, where the same sum
reproduces the textbook Christoffel source term by term.

    null (b0c=0)            1.1e-17, unchanged -- the isotropic balance survives
    hydro rigid rotation    unchanged to 4 digits -- the rho*v*v part is right
    convergence             ~1.9, so it is consistent
    per-step L1(v_t), PLM   2.07e-05 -> 2.33e-05 (interior).  NO GAIN.
    per-step L1(v_t), PPM4  9.93e-07 -> 6.12e-07 at nx=64.  1.6x.
    STABILITY, t(KE/V=1e-6) 0.245 -> 0.243 (n16, beta .05);  0.059 -> 0.048 (n16, .02)
                            1.284 -> 1.081 (n32, beta .05).  NO CHANGE.

**So the hypothesis that the source/flux quadrature mismatch drives the feedback is
REFUTED.** With PLM the residual is dominated by the reconstruction and by the CT face
average of B.nhat, not by the source's quadrature -- which is exactly why the gain only
appears once PPM4 has shrunk the reconstruction error. It also costs 1.68x wall clock
(five PanelTangents per cell per stage; removable by caching the triads per (m,k,j), since
they depend on the two ANGLES alone).

## WHERE TO GO NEXT

Neither the truncation MAGNITUDE nor the source FORM changes the growth: ppm4 shrinks the
residual and does not slow the instability, and the well-balanced source does not either.
A growth rate is set by the operator's response to a perturbation, not by the size of the
seed, so what is left is the part of the operator that acts ON perturbations: the Riemann
solver's dissipation on the TANGENTIAL field components (recall the radial monopole is
stable on both grids, so the tangential components are the ones that matter) and the CT /
corner-EMF coupling. The next measurement is a dissipation A/B at fixed everything else:
hlld vs hlle vs llf on the beta ladder. If a more dissipative solver stabilises it, the
scheme is under-diffusive for the tangential magnetic modes on a curvilinear grid; if none
of them does, the feedback is in CT.

Gate any fix on the fit-free table above, at matched beta, resolution and time -- not on a
growth-rate fit.

## How to run either grid

    # spherical polar (build_sp, -D PROBLEM=sp_lowbeta)
    build_sp/src/athena -i inputs/tests/spherical_polar_mhd_lowbeta.athinput -d . \
      mesh/nx2=$n mesh/nx3=$((3*n)) meshblock/nx2=$n meshblock/nx3=$((3*n)) \
      problem/b0c=$b time/tlim=3.0
    # cubed sphere (build_cs, -D PROBLEM=cs_test)
    build_cs/src/athena -i inputs/tests/cubed_sphere_mhd_strat.athinput -d . \
      problem/iprob=8 problem/user_srcs=false mesh/nx1=16 meshblock/nx1=16 \
      mesh/nx2=$n mesh/nx3=$n meshblock/nx2=$n meshblock/nx3=$n \
      problem/b0c=$b time/tlim=3.0 output1/dt=0.01
    # cs volume 9.948, sp wedge 4.974;  beta = 2*p0/b0c^2 = 0.2/b0c^2
