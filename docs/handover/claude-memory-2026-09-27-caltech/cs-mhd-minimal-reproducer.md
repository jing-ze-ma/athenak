---
name: cs-mhd-minimal-reproducer
description: THE MINIMAL REPRODUCER of the cs+MHD blow-up -- cs_test iprob=13, a hydrostatic atmosphere AT REST carrying a force-free field, dies in 2.5 min on one CPU core. Growth rate rises with resolution. Halo, stratification and shear all EXONERATED
metadata:
  type: project
---

**SUPERSEDED IN PART, 2026-09-03 by [[cs-mhd-instability-characterized]]:** three things
below are RETRACTED -- the vertex order deficit (a shrinking-bin artifact: the bins were 2
cells wide), the mode being vertex-localised, and **"the growth rate RISES with
resolution"**, whose 13.8-vs-17.3 pair is a least-squares slope over run-dependent windows;
fit-free, refinement DELAYS onset. The defect is an INSTABILITY, not a consistency loss,
and it is not gnomonic-specific -- spherical polar does it too. The reproducer, the
tangential-velocity instrument and the four eliminations all stand.

**2026-09-03. Committed as f3d35a96.** `cs_test` iprob = 13 +
`inputs/tests/cubed_sphere_mhd_strat.athinput`. This SUPERSEDES the iprob=9 reproducer in
[[cs-mhd-low-beta-divergent]], whose defect was fixed by a46b75e0 while the production run
kept dying.

## The problem

Constant inward gravity (a user source term -- `<srcterms>` is refused on cs, but a purely
RADIAL source needs no gnomonic form because rhat is orthogonal to BOTH angular basis
vectors), an isothermal profile `rho = d0 exp(-(r-r0)/H)`, and the SAME uniform Cartesian
field iprob = 8 uses. The field is force-free at any beta, so **the exact solution is the
atmosphere sitting still forever and v is exactly zero in it.** p falls by exp(-5) over the
domain while |B| does not, so ONE RUN sweeps beta from 2.0 at the base to 1.3e-2 at the top.

## IT BLOWS UP -- from truncation noise alone

nx1=16, b0c=3.162 (beta 0.02 base, 1.3e-4 top), NOTHING else: no RT, no resistivity, no
rotation, no shear, no general EOS, no grid stretch, ideal gas, uniform grid, at rest.

    nx2=nx3=16   KE 1.2e-05 -> 1.8e-03 over t = 0.02..0.20;  dt 2.6e-4 -> 1.0e-4;
                 DEAD (dt = 0) at t = 0.2013;   amplitude growth rate 13.8
    nx2=nx3=32   KE 9.7e-07 -> 6.1e-05 over t = 0.02..0.14;  dt falling;  rate 17.3

**The GROWTH RATE RISES WITH RESOLUTION** -- a grid-scale numerical INSTABILITY, matching
the production run's "refining makes it fail sooner". The finer run starts from a smaller
seed, so which arm dies first depends on seed AND rate; do not read death times alone as a
resolution trend.

**~2.5 minutes on one CPU core.** That is the whole point.

## The instrument: the TANGENTIAL velocity norm

`CSTestConvErrors` now prints `TANGENTIAL ONLY  L1(v_t)` per region. Essential: the radial
hydrostatic imbalance is the largest error and is IDENTICAL in all three regions, so it
buries the vertex signal -- interior/seam/vertex agreed to 3 digits until the norms were
split. Clean NULL: at b0c = 0, L1(v_t) is 4e-18..2e-17 everywhere, so the entire signal is
magnetic. Cf. [[validate-the-instrument]].

Per-step (nlim=1) L1(v_t), b0c = 1:

    nx      interior      seam        CUBE VERTEX
    16     2.5144e-05   4.9540e-05   8.4388e-05
    32     8.0831e-06   1.8381e-05   3.3684e-05
    64     2.3098e-06   6.1945e-06   1.4933e-05
    order    1.65/1.81    1.43/1.57    1.33/1.17   <- vertex ~1.2 and FALLING

## Four things ELIMINATED by direct A/B, not by inference

1. **The HALO is not the limiter.** `problem/exact_panel_ghosts=1` now has an MHD form that
   replaces the panel ghosts -- state AND both angular face fields -- with analytic values
   at the ghost's own chart-continued (xi,eta), i.e. a PERFECT halo. Orders go 1.30/1.19,
   unchanged. (It DID take effect: seam and vertex values move while the interior stays
   identical to 5 digits -- exactly the right pattern, and the check that it is not a no-op.)
2. **Stratification is irrelevant.** iprob = 8, no gravity, same beta: vertex 1.34/1.18,
   interior 1.65/1.82. The same orders. So this is a plain vertex defect in the ideal
   operator, NOT something stratification introduces.
3. **Shear/rotation/RT/EOS/resistivity/stretch** -- none are present at all, and it still dies.
4. **Reconstruction carries PART of it.** ppm4 (nghost=4) lifts the vertex to 1.28/1.64 and
   the seam to 1.60/1.91, and halves the absolute vertex error. Part, not all -- the vertex
   still trails the interior. A valid substitution because ppm is MORE accurate than the
   effect being chased; cf. [[cs-mhd-low-beta-divergent]]'s retracted donor-cell family.

The error is LINEAR in |B| over b0c = 0.1..3.162, so the ACCELERATION error is QUADRATIC --
as a Lorentz-force cancellation error must be, with dt ~ 1/B supplying the other power.

## How to run it

    build_cs/src/athena -i inputs/tests/cubed_sphere_mhd_strat.athinput -d . \
      mesh/nx2=$n mesh/nx3=$n meshblock/nx2=$n meshblock/nx3=$n time/nlim=1 problem/b0c=1.0
    # ... and read the "TANGENTIAL ONLY" rows.  For the blow-up:
    #   mesh/nx1=16 meshblock/nx1=16 problem/b0c=3.162 time/tlim=2.0 output1/dt=0.02

## WHERE TO GO NEXT

The defect is in the INTERIOR OPERATOR at skew vertex cells (halo exonerated, reconstruction
only partial). The suspect that fits every fact: on a curvilinear FV grid the momentum flux
is rotated back to covariant slots using EACH FACE's OWN trig, while the geometric source
term uses the CELL-CENTRE trig -- the two must cancel, and the basis rotation per cell is
largest exactly at a cube vertex. Both were verified ALGEBRAICALLY correct
([[cs-mhd-low-beta-divergent]]); what has NOT been checked is whether their DISCRETE
cancellation is second order there. The measurement to make: evaluate the momentum RHS on
the exact initial state and split it into flux divergence vs geometric source, per region,
and check the order of their SUM.

Note hydro is clean in this test only because the hydro state is a function of r alone, so
every tangential term vanishes by symmetry -- not because the operator is proven good.
