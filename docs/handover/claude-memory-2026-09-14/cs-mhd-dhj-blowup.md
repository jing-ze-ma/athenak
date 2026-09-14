---
name: cs-mhd-dhj-blowup
description: The cs hot-Jupiter MHD blow-up. SUPERSEDED as the working front by cs-mhd-minimal-reproducer (f3d35a96), which reproduces it in 2.5 min on one core; this file is the elimination trail
metadata:
  type: project
---

**SUPERSEDED IN PART, 2026-09-03 by [[cs-mhd-instability-characterized]]:** "at a cube
vertex" is RETRACTED -- the growing mode fills the panel interior. It is a low-beta
instability that needs a TANGENTIAL field, and it is NOT cubed-sphere specific: spherical
polar does the same thing at the same beta, only ~10x weaker.

**2026-09-03: SUPERSEDED as the place to work.** The failure now reproduces in a 2.5-minute
single-core test with none of this problem's physics -- see [[cs-mhd-minimal-reproducer]]
(cs_test iprob = 13, f3d35a96), and go there first. It is a grid-scale INSTABILITY whose
growth rate rises with resolution, in the interior operator at skew cube-vertex cells.
Everything below stays valid as the ELIMINATION TRAIL for the production run; the dhj run
itself is still unfixed (a46b75e0 moved it only 0.2023 -> 0.2129 rot).

**OPEN, 2026-09-02.** Found by checking `cs_prod_mhd` (job 11311755, cancelled). The
cubed-sphere deep-hot-Jupiter run blows up at **~0.2 rotations with MHD**, on HEAD
`a4d30d43` -- i.e. WITH the RaiseVel floor fix (5bceca9d) and the seam co-location fix
(999458a3), both of which postdate the earlier cs failures.

## THE THREE-ARM DISCRIMINATOR: it is NOT the resistivity

Bench `/viper/u2/jinma/ATHENAK/bench/cs_etatest`, 0.3-rotation arms on apudev, hst every
0.01 rotation:

    arm                       outcome
    max_eta = 1e13            DIED at 0.2000 rot
    max_eta = 1e12            DIED at 0.1927 rot
    IDEAL MHD (no <mhd>/ohmic_resistivity at all)   DIED at 0.2023 rot

**All three die within 4 % of the same time, including ideal MHD.**  So the resistivity is
innocent and the eta cap is innocent: I had raised `max_eta` 1e12 -> 1e13 to match
`sp_dhj_ctl`, and blamed that change; **that hypothesis is RETRACTED.**  The defect is
cubed sphere + MHD on this problem, with or without eta.

**The control is clean:** `cs_prod_hyd` (same grid, same binary, same everything but the
`<mhd>` block) is healthy past **11.2 rotations**, 0 NaN, mass drift -0.71 %, dt 17.5.
And spherical polar + the SAME MHD physics at max_eta = 1e13 (`sp_dhj_ctl`) ran 95
rotations clean.  So: grid x physics -- only (cs, MHD) fails.

## Symptoms

* A FLOOR STORM, not an EOS failure: `eos_fail` stays 0 and `c2p_it` 0 throughout, while
  `eos_efloor` runs away by four orders of magnitude and eventually INTEGER-OVERFLOWS the
  counter (-1135976192).  Energy is being driven non-positive and floored, repeatedly.
* dt is squeezed BEFORE the collapse: 28.98 -> 9.15 by cycle 2100, while the hydro arm
  holds 22-25 at the same times.  Then dt freezes (10.4466) and the job grinds on with a
  NaN state -- the [[cs-dhj-long-run]] signature, and it will burn its whole allocation.
* e12 and e13 have IDENTICAL mass to 6 digits (3.43859e26) at the last common output, so
  the eta cap barely perturbs the solution before the blow-up.  It is a sudden event.

## What to do next

1. Dense output through 0.15-0.21 rot on one arm (bin dumps every ~0.005 rot) and find
   WHERE it starts -- panel interior, seam, or cube vertex.  The whole failure is ~2000
   cycles, a couple of minutes of wall time, so this is cheap.
2. The ideal-MHD arm is the one to debug: it removes the entire resistive path and still
   fails, so the suspect list is the ideal cs MHD operators (mhd_fluxes gnomonic
   rotation, corner EMF, the fc seam halo) interacting with a STRATIFIED atmosphere --
   note every cs MHD validation so far ([[cs-mhd-validation]], [[cubed-sphere-mhd-seam]],
   [[cubed-sphere-seam-emf]]) used an unstratified test problem.
3. Check whether the field is even dynamically important here: `ck_mhd_b3` recorded that
   the field "doesn't bind" dt at bbot 3 G.  If B is passive, an O(1) error in it can
   still wreck the energy through the ME term in the conserved total.

**Do NOT re-blame max_eta, and do not re-run the eta scan** -- three arms already settled
it.  Cf. [[measure-impact-before-claiming]].

## A trap this run produced

At an intermediate check the arms read "e12 CLEAN to 0.190 rot, ideal CLEAN to 0.150" and
I nearly reported that the higher cap was the culprit.  They were simply **not yet at the
failure time**.  A partial arm is not a result; wait for every arm to reach the same
simulated time before comparing.  Cf. [[validate-the-instrument]].


## 2026-09-02, later: LOCALISED, and four mechanisms eliminated

**Where.** Dense dumps every 250 s through the onset (bench `cs_onset`, restart at 0.180
rot). Max Alfven speed split into corner / seam / interior bands (2 cells wide):

    t=5.53e4   corner 4.88e6   seam 5.51e6   interior 5.62e6   <- corner HEALTHIEST
    t=5.85e4   corner 7.75e6   seam 6.81e6   interior 7.01e6   <- corner takes the lead
    t=5.90e4   corner 1.39e7   seam 6.56e6   interior 7.31e6   <- 2x away and running

Ignition is at a **cube VERTEX, radial index i ~ 80**, at t = 5.85e4; the whole domain is
NaN by 6.008e4, 0.008 rotations later. The seam band tracks the interior throughout, so
this is NOT a general seam defect.

**RETRACTED as symptoms.** The "floor storm" is normal -- the healthy cs hydro run and the
95-rotation sp control log the same 1e8-scale `eos_dfloor` counts. So is the 33x field
amplification: it is identical in all bands including the deep interior, and sp's
magnetic energy grows comparably.

## ELIMINATED

1. **Resistivity** -- ideal MHD dies too (the earlier three-arm test).
2. **The vertex-fill ALGORITHM.** Same restart, only `mesh/cs_vertex_fill`:
   real-data widened buffers die at 6.008e4, quadratic extrapolation at 5.880e4. The real
   fill is better, in the expected direction, by **2%**. Not the cure.
3. **The TIMESTEP.** A CFL scan from the same restart: cfl 0.300 -> first NaN 6.00760e4;
   **cfl 0.150 -> 6.01314e4, a shift of 0.09%**. (cfl 0.075 ran out of wall time at
   5.98e4 without reaching the failure.) The dt collapse begins at the SAME simulated time
   ~5.85e4 in every arm and the arms' dt ratio is exactly the CFL ratio, so the same
   physical speed limits all of them. **The runaway is dt-independent: spatial, not a
   stability margin.**
4. **CT / magnetic monopoles.** With `mhd_divb` fixed to the finite-volume form (see
   below), divB through the whole ignition window is round-off and the vertex is
   INDISTINGUISHABLE from the interior:

       t=5.50e4  vertex 5.0e-21  interior 4.4e-21
       t=5.85e4  vertex 6.7e-21  interior 4.7e-21     <- ignition
       t=5.95e4  vertex 6.0e-21  interior 5.9e-21

   Flat while the corner Alfven speed doubles. **Constrained transport is exact at the
   cube vertex.**

   **The control makes this emphatic.** The same diagnostic on the HEALTHY sp run
   (`sp_dhj_ctl`, 95 rotations clean) gives max |divB| = **6.16e-14**, and it sits exactly
   ON THE POLAR AXIS (theta index 0), that grid's own geometric singularity, where the
   cell volume goes to zero. So the cs cube vertex is **seven orders of magnitude cleaner
   in divB** than the spherical-polar pole -- and the polar run survives. divB is not the
   discriminant, in either direction. So divB is the wrong gate here too, but for a different reason than at
   the seam ([[cubed-sphere-seam-emf]]): the EMF is single-valued, it is just not proven
   accurate. A wrong-but-consistent EMF still gives divB = 0.
5. The static uniform-field corner gate (`cs_test` iprob=8) is clean: dp/p = 2.0e-4
   (the documented O(dx^2)), v = 0 and EMF = 0 exactly. Too easy a test to discriminate --
   a uniform field is the case the halo transform handles best.

## TWO REAL DEFECTS FOUND ALONG THE WAY, NEITHER THE CAUSE

* **The cs timestep uses the wrong cell width.** `dx2_ = r_c*dth_xi`, `dx3_ = r_c*dth_eta`
  are arc lengths ALONG the two panel tangents, but those are non-orthogonal
  (`e1.e2 = -sin(xi) sin(eta)` = `cos_cell`). The face separation is `dx*sin_cell`, and
  `sin_cell` is 1.000 at a panel centre but **0.880 in the corner cell**, so dt is
  overestimated by up to 13.7% -- effective CFL 0.341 against the requested 0.300 --
  with the error peaking EXACTLY at the cube vertex. Seductive, and the CFL scan rules it
  out. `dx1_/dx2_/dx3_` are read by nothing but `hydro_newdt.cpp` and `mhd_newdt.cpp`, so
  the fix is contained. **NOT APPLIED** -- it changes every cs answer slightly.
* **`mhd_divb` was meaningless on every curvilinear grid** -- it differenced face fields
  over `size.dx1`, the block's uniform CARTESIAN extent. FIXED to
  `sum(area*B_face)/volume` for cubed-sphere and spherical-polar, which is what CT
  conserves. The Cartesian form is nonzero everywhere by geometry and would have hidden
  exactly the monopole the diagnostic exists to find. This fix is what made test 4
  possible.

## 6. The FIELD SEAM TRANSFORM is exact, into the corner

`cubed_sphere::TransformFieldToDstNormals` unit-tested standalone (strip the athena.hpp
include from the header and compile it alone -- it needs only `Real` and
`KOKKOS_INLINE_FUNCTION`). Project a known Cartesian field on the SOURCE panel's face
normals, push it through, compare against the direct projection on the DESTINATION
normals at the same physical point. Swept from the seam midpoint into the cube VERTEX,
where the inter-chart shear is O(1), for a uniform, a dipole and a sheared field:
**exact to 1e-16 everywhere, corner included.** The math in that header is right.

## 7. The MeshBlock DECOMPOSITION is irrelevant -- the block halo is EXONERATED

Same bbot = 3, ideal MHD from t = 0, only `meshblock/nx2,nx3` changed:

    1 MeshBlock per panel   (6 blocks)   DIED at t=6.169770e4 = 0.202288 rot
    2x2 MeshBlocks per panel (24 blocks) DIED at t=6.170160e4 = 0.202300 rot

**0.006% apart**, and not bit-identical (different t), so these are genuinely different
runs converging on the same failure. Removing EVERY intra-panel block boundary changes
nothing. Verify the override took: the log must say `Root grid = 6 x 1 x (1 x 1)` and
`Total number of MeshBlocks = 6`. **So the defect is in what BOTH decompositions share --
the panel seam and cube vertex geometry itself -- not in the block-to-block halo.**

## 8. bbot THRESHOLD SCAN: between 1 and 3 G

cs, ideal MHD, from t = 0, tlim 0.3508 rot (bench `cs_bbot`):

    bbot = 3.0 G   DIED at 0.2023 rot   (reproduces the eta-test arm to 4 digits)
    bbot = 1.0 G   clean to 0.3150 rot  (stopped by wall clock, not tlim)
    bbot = 0.3 G   SURVIVED the full 0.3508 rot
    bbot = 0.1 G   SURVIVED the full 0.3508 rot

**Do NOT read this as "marginality, no bug".** The user pushed back on exactly that and was
right: spherical polar has the SAME beta ~ 1e-3 atmosphere, the same field, the same input
file, and runs 95 rotations -- a 475x difference in survival. A physical margin does not
care which chart it is written on. A cs bug injecting an error proportional to |B| at the
vertex produces this identical curve, so the threshold measures SEVERITY, not cause.

## 9. The HALO ASSEMBLY is clean, and the GENERAL EOS is exonerated

**Along-seam resample map**: `atan(tan(a)*tan|n|)` unit-tested against the charts
themselves -- take the ghost's physical point, read its along-seam angle in the SOURCE
chart, compare. Seam midpoint to CORNER, both ghost layers: **exact to 2e-15 cells**.
The angle->index conversion by a constant `dang` is exact too, not a linearisation,
because the grid is equiangular.

**Staggering** (`angles`, bvals_fc.cpp): `vv==1` (b.x2f) takes the xi FACE and eta CENTRE,
`vv==2` the reverse, `vv==0` both centres -- right per component; and the stencil bound
`bhi = ke_ + (vv==2 ? 1 : 0) - 2` extends by one exactly where the along-seam axis is
face-sampled. Self-consistent.

**Vertex EMF**: `AveragePanelCornerEMF` already does a genuine THREE-way average
`(own + bj + bk)/3` over the three panels meeting at a vertex, with the wrong-diagonal
generic pairing deliberately skipped. Correct.

**The EOS is NOT it.** The gate that certifies the cs magnetic-energy inversion
(`cs_test` iprob=8) runs with `eos = ideal`, while the failing run uses `eos = general` +
table -- so general EOS x cs x MHD was UNGATED, and that cache once carried a cs-specific
bug. Tested by re-running the same grid, same bbot=3, same ideal MHD with an IDEAL GAS
(bench `cs_idealgas`; verify with "ideal EOS, via <units>/mu" in the log and ZERO
"General EOS: tabulated" lines):

    general EOS (table)   DIED at 0.2023 rot
    IDEAL GAS             DIED at 0.1604 rot   <- dies EARLIER

Dying anyway is the decisive branch of that test, so **the tabulated EOS, the cached
p/Gamma_1 and the general-EOS magnetic-energy inversion are all exonerated.**

## 10. STRIP-DOWN batch 1: radial stretch eliminated, grid asymmetry survives

All arms ideal gas + ideal MHD + bbot = 3 (bench `cs_strip`, `cs_idealgas`):

    cs, stretched radial grid    DIED at 0.1604 rot
    cs, UNIFORM radial grid      DIED at 0.1506 rot   <- stretch ELIMINATED
    spherical polar, ideal gas   SURVIVED the full 0.3508 rot

The ignition sits at i ~ 80 where the polynomial stretch is active, so the stretch was a
natural suspect; removing it entirely moves the failure 6%, and dying anyway is the
decisive branch. **And the cs-vs-sp asymmetry is not an EOS artifact**: sp with an ideal
gas is clean through the whole run while cs dies at 0.15.

Strip axes deliberately NOT used first, because they are CONFOUNDED -- removing rotation
or lowering bbot merely stops the field winding up, so surviving proves nothing. Prefer
axes that change the DISCRETISATION at fixed physics: grid stretch, angular resolution,
reconstruction order.

## 11. STRIP-DOWN batch 2: the failure gets WORSE with resolution

Base = the most stripped failing setup (cs, ideal gas, ideal MHD, uniform radial grid,
bbot = 3, 1 block/panel), which dies at 0.1506 rot.

    nx2=nx3=32 (base)   DIED at 0.1506 rot
    nx2=nx3=16 (half)   DIED at 0.1901 rot   <- COARSER survives LONGER

**Refining makes it fail sooner.** That is the opposite of an under-resolved feature
converging away, and it rules out "not enough cells" as the cause. It is what an
INCONSISTENT operator, or a grid-scale mode whose growth rate rises with 1/dx, looks like.
Only two points, so treat the 26% as a direction, not a scaling law.

**The donor-cell arm is CONFOUNDED -- do not lean on it.** `mhd/reconstruct=dc` was clean
to 0.2588 rot, but (a) it stopped on the WALL CLOCK, not tlim, so it is not a confirmed
survival, and (b) it is not the same problem: at t = 9.15e3 its radial KE is 3.03e35
against the base's 1.23e32, a factor 2500, because first order cannot hold the
hydrostatic balance. Surviving proves little when the solution itself has changed that
much. The override DID take effect -- that KE ratio is the proof -- so the arm is valid,
just uninformative. The decisive branch would have been dying.

## WHAT IS LEFT

It is a **cs bug in the core IDEAL-MHD operators**, not in the halo, not in the EOS, not
in CT, not in the timestep. Everything communicated or transformed has been verified
exact; what has NOT been verified on a STRATIFIED, low-beta state is the interior
operator chain at a skew cell: `mhd_fluxes` gnomonic rotation, the cs corner EMF
(e63b571a), and `GnomonicEquiangleRaiseVelMHD`'s orthonormal bcc (which feeds
reconstruction and the fast speed, and carries 1/sin_cell = 1.14 at the vertex).

**NEXT: build a MINIMAL REPRODUCER.** The ideal-gas arm is the vehicle -- it fails FASTER
(0.16 rot) and is far cheaper than the tabulated EOS. Keep stripping the problem (RT off,
simpler atmosphere, uniform gravity, no rotation) until it stops failing; the last thing
removed is the cause. A failure this robust to everything tried so far should survive a
long way down, and a small reproducer is what makes the operator bug findable.
Cf. [[validate-the-instrument]].
