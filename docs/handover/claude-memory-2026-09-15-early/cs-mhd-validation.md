---
name: cs-mhd-validation
description: IDEAL and RESISTIVE MHD on the cubed sphere by region. RE-MEASURED on HEAD 999458a3: ideal MHD conserves to machine precision and is 2nd order EVERYWHERE (the old 1.84/1.76 was a fixed-radial-resolution floor); the resistive energy "drift" is PHYSICAL Ohmic heating supplied through the radial boundary, matching eta*J^2*V*t to 0.1-0.3%, NOT a bug
metadata:
  type: project
---

Measured 2026-08-31, `build_cs`, rk2 + plm + hlld. Companion to [[cs-hydro-validation]].

## THE TRAP THAT COST THE FIRST ROUND: a uniform field THREADS a reflecting wall

My first MHD conservation test put reflecting radial walls on a state with a UNIFORM tilted
field. Mass grew 6%. That is **not a cubed-sphere defect and not a code bug**: a field with
B.rhat != 0 at the wall is not a closed system in MHD, and the ill-posed wall Riemann
problem leaks mass at O(B^2). Verified: the leak scales as B^2 (2.29e-04, 2.18e-06,
-6.2e-15 for b0c = 0.1, 0.01, 0) with **ZERO floor events**, and the SAME behaviour
reproduces on a PLAIN CARTESIAN grid with AthenaK's own `inputs/mhd/blast_mhd.athinput`
(mass +7.5e-04 weak field, +1.6e-02 strong, exactly 0.0 at b_amb = 0).

**The fix: `iprob=11` with `problem/bunif=0` is purely AZIMUTHAL, B.rhat = 0 EXACTLY, so
reflecting walls pass no flux and a global budget is meaningful.** `CSTestResistCheck` now
runs the conservation/seam/level gates BEFORE the resistivity guard, so dropping
`<mhd>/ohmic_resistivity` turns it into the IDEAL MHD closed test.

## CRITERION 1 -- CONSERVATION (closed test, tlim = 0.2)

```
  ideal MHD          mass -1.65e-14   energy  4.79e-14     <- MACHINE PRECISION
  resistive eta=0.5  mass  1.77e-14   energy  8.84e-04
  resistive eta=0.05 mass -1.39e-14   energy  1.16e-04
  resistive eta=5e-3 mass  8.85e-15   energy  1.22e-05
```

Mass is always round-off. The resistive ENERGY drift is proportional to eta, linear in t,
and **does NOT converge with resolution** (1.164e-04 at nx2=16 vs 1.193e-04 at nx2=32,
order -0.04). Two readings, NOT yet separated: a real Poynting flux out of the reflecting
wall (E_t = eta*J_t != 0 there, which would be physical given the BC and would not
converge), or a non-conservative resistive energy update. **Discriminator not yet run:**
dilute it -- grow the radial extent and see whether the relative loss falls like
area/volume (boundary) or stays (bulk). Ideal MHD with the same walls is exact, so
whatever it is, it is the eta*J channel.

**Still missing: a CLOSED MHD SHOCK test.** The blast uses a uniform field, which threads
the walls. It needs an azimuthal-field blast (A_z = -0.5*b*(x^2+y^2) added to iprob=12's
vector potential) -- a small pgen addition, not yet made.

## CRITERION 2 -- CONVERGENCE, BY REGION (this is the question that matters)

`CSTestConvErrors` and the evolved-field gate now split every norm into **panel INTERIOR /
panel SEAM / CUBE VERTEX**, classified by ANGLE (within 2 cells of |xi| = pi/4 and/or
|eta| = pi/4) so it is independent of how many MeshBlocks span a panel. A domain norm
cannot answer "second order everywhere": the seam is codim-1 and the vertex codim-2, so a
first-order region a few cells wide holds a vanishing share of L1.

```
                        L1 order (16->32, 32->64)      Linf order
  HYDRO      interior      2.48, 2.31                   2.02, 2.11
             seam          1.99, 2.18                   1.79, 1.90
             CUBE VERTEX   1.58, 2.08                   1.85, 2.29     <- 2nd order EVERYWHERE

  IDEAL MHD  interior      2.11 (dt->0 extrapolated)    1.32
   (dt->0)   seam          1.84                         1.15
             CUBE VERTEX   1.76                         1.22

  RESISTIVE  interior      1.89, 2.03                   0.82, 0.90
   (eta=1e-3) seam         1.29, 0.99                   0.48, 0.32
             CUBE VERTEX   0.68, 0.64                   0.52, 0.54     <- NOT 2nd order
```

**IDEAL MHD MUST BE EXTRAPOLATED TO dt -> 0 FIRST.** At a fixed CFL the raw order reads
~1.1-1.35 in ALL THREE regions equally -- a BULK dt-proportional error, the known radial-BC
ghost-clock phase lag ([[cubed-sphere-mhd-convergence]]), not a seam defect. Halving the
CFL at nx2=32 changes L1(B) by 32%. Use `E0 = 2E(dt/2) - E(dt)`.

**RESISTIVITY IS THE REAL GAP:** second order in the panel interior, ~first order on the
seam, and ~0.65 at the cube vertex, with Linf barely converging anywhere. The cube-vertex
fill ([[cs-wire-fill-wip]]) lowers the vertex error 9-17% but **does not change the
ORDER** (0.68/0.64 on, 0.73/0.72 off).

## METHOD NOTES

- **Pick eta so dt is NOT diffusion-limited** (the user's point): at eta = 0.5 the
  diffusive dt ~ h^2 makes nx2=64 a 3-hour run; at eta = 1e-3 the hyperbolic CFL binds and
  the same series takes minutes. The gate normalises by |eta*2*b0c|, so the test stays
  meaningful.
- Refine EVERY direction together, or the fixed radial resolution is a floor
  ([[cs-hydro-validation]]).
- `pkill -f` from inside the Bash tool kills the calling shell too (exit 144). Write the
  script first, launch with `setsid nohup ... < /dev/null &`.

## RE-MEASURED on HEAD 999458a3, 2026-09-02.  Bench `/viper/u2/jinma/ATHENAK/bench/cs_reval`

Everything above reproduces **digit for digit** except the ideal-MHD region orders, which
were WRONG for a protocol reason.  Hydro L1(v) 3.3954e-04 / 6.5117e-05 / 1.2279e-05 and
the resistive evolved field 7.6807e-07 / 1.9097e-07 / 4.7433e-08 (global order 2.01/2.01)
came back identical to the recorded values, as did every d0..d5 bin.

**CORRECTION -- IDEAL MHD IS 2nd ORDER EVERYWHERE, INCLUDING THE SEAM AND THE VERTEX.**
The recorded 2.11 / 1.84 / 1.76 was measured with **nx1 FIXED at 8** while the tangential
axes refined.  Refining all three jointly (nx1 = nx2/2) and extrapolating dt -> 0:

  L1(B), dt->0        16->32   32->64
    panel INTERIOR      2.49     2.36
    panel SEAM          2.12     2.18
    CUBE VERTEX         2.06     1.99

**This is NOT the co-location fix.**  The A/B settles it: the pre-fix binary
(`bench/base_wt`, 5bceca9d) run through the SAME joint protocol gives 2.49/2.36,
2.12/2.19, 2.06/2.01 -- indistinguishable.  The old seam/vertex numbers were the fixed
RADIAL resolution acting as an error floor, i.e. exactly the trap already written up in
[[cs-hydro-validation]], which had bitten this table itself.  I nearly reported the
improvement as a win for the seam fix; the A/B is what stopped it.
Cf. [[measure-impact-before-claiming]].

The resistive series does NOT have that confound: refining radially too (nx1 8->16 with
nx 32->64) moves it 2-7% (global 2.01 -> 2.04, seam 1.60 -> 1.64, vertex 1.37 -> 1.46),
confirming the radial floor already refuted in [[cs-seam-order-limiter]].

## THERE IS NO RESISTIVE ENERGY BUG.  THE DRIFT IS PHYSICAL OHMIC HEATING

**RESOLVED, and my own "bulk non-conservation" reading above is RETRACTED.**  The state
is exactly static -- curl B is the UNIFORM vector 2*b0c*zhat, so curl(eta J) = 0 and B
never moves -- yet Ohmic heating eta*J^2 is deposited in every cell.  That energy has to
come from somewhere, and it comes in as a Poynting flux through the RADIAL boundary.  The
expected gain over the run is therefore exactly `eta*J^2*V*t`, and the code delivers it:

    eta      BC        dE/E measured   eta*J^2*V*t/E    ratio
    0.5      user        2.6975e-03     2.6985e-03      1.000
    0.05     user        2.6910e-04     2.6985e-04      0.997
    0.005    user        2.6405e-05     2.6985e-05      0.979
    0.5      reflect     8.8360e-04     2.6985e-03      0.327
    0.05     reflect     1.1643e-04     2.6985e-04      0.431
    0.005    reflect     1.2218e-05     2.6985e-05      0.453

and it holds under refinement and a 9x change of volume (nx2=64: 0.998; x1max=4, nx1=24,
V = 263.9: 0.999).  **A `reflect` wall is a perfectly CONDUCTING wall, E_t = 0**, which is
inconsistent with E_t = eta*J_t != 0 here, so it under-delivers the influx by 0.33-0.45 --
that deficit is the entire "non-conservation" this file used to report.  The exact-state
`user` BC gets it right to 0.1-0.3%.

The update is strict flux form (`AddResistiveFluxes` runs INSIDE `MHD::Fluxes`, i.e.
BEFORE SendFlux, so the seam reconciliation sees it), and the gates confirm the interior
telescopes: SEAM FLUX MISMATCH dE/dt = 1.0e-18, LEVEL-BOUNDARY clean.  With those two
exact, the only faces that can source energy are the physical radial walls -- which is
what the eta scan then measured.

**`CSTestResistCheck` now PRINTS `eta*J^2*V*t`** next to the conserved sums, with the
explanation in the source, so the gate can never be read as a defect again.

**HOW I GOT IT WRONG, and the lesson.**  I first ran a DILUTION test (grow the radial
extent, see whether the relative drift falls like area/volume) and read "it did not fall,
therefore bulk".  **That test is VACUOUS for this setup and the reasoning was wrong:**
with p0 = 1 and b0c = 0.1 the energy is THERMAL-dominated, so total energy goes as r^3,
and the wall term, the seam term and the bulk term ALL go as r^3 too -- every hypothesis
predicts a flat ratio, so a flat ratio discriminates nothing.  Check that the candidate
hypotheses actually predict DIFFERENT things before running the discriminator.
Cf. [[localise-by-dilution]], [[confounded-tests-rejected]], [[measure-impact-before-claiming]].
The test that DID work was comparing the drift against an ANALYTIC prediction over three
decades of eta -- an absolute number, not a scaling.

## PROTOCOL: THE CLOSED TEST NEEDS `mesh/ix1_bc=reflect` ON THE COMMAND LINE

Every shipped cs input except the blast has `ix1_bc = user` -- the exact-state DIRICHLET
BC, which is **not closed**.  Run the ideal-MHD conservation gate as shipped and it reads
mass -2.2e-10, energy -3.5e-7; add `mesh/ix1_bc=reflect mesh/ox1_bc=reflect` and the same
run gives +2.2e-15 / +2.6e-14.  Five to seven decades, purely from the BC.  I read the
first as a regression for several minutes.  Also: an arm that changes the DOMAIN needs its
OWN `time/nlim=0` reference, or the drift is nonsense (I got mass +8.0).
