---
name: cubed-sphere-smr
description: f14e075d + 608d43e1 -- STATIC refinement works on the cubed sphere for HYDRO with RADIAL level boundaries: conservative to round-off, and verified on CPU, MPI (1/2/3/6 ranks) and GPU (1 and 2 MI300A). Records why no coarse-mesh geometry was needed, two cross-panel bugs, and exactly what is still refused
metadata:
  type: project
---

**WORKS as of 2026-08-29, f14e075d + 608d43e1.** Hydro static refinement, level
boundaries RADIAL. `inputs/tests/cubed_sphere_smr.athinput` is the documented example.
Verified on **CPU serial, MPI at 1/2/3/6 ranks, and GPU on 1 and 2 MI300A**.

## Measured (rigid rotation, reflecting radial BCs = exactly zero physical flux, t=0.25)

| | with the fix | unweighted operators |
|---|---|---|
| relative mass drift | **-9.6e-14** | 9.1e-09 |
| relative energy drift | **1.0e-13** | 1.5e-08 |

Round-off, and 5 orders better than the plain averages -- that A/B is what says the
weighting is doing the work. The UNREFINED case is bit-identical to before.
Lz is NOT exact and cannot be (geometric source terms; not a telescoping invariant).

## The key insight: NO coarse-mesh geometry is needed

`Coordinates` computes volume/area only on the fine mesh, and I first scoped adding
coarse versions. Not necessary:

* **Flux correction**: send the AREA-WEIGHTED sum `sum(A_f F_f)`; the coarse side divides
  by **its own** face area. The coarse update then multiplies by that same area, so the
  coarse cell removes exactly what the fine cells removed. **This works even though the
  coarse area does NOT equal the sum of the fine areas** -- the gnomonic solid angle is
  not additive under tangential refinement.
* **Restriction**: volume-weight and normalise by the **sum of the eight fine volumes**,
  not by a coarse volume.

Both use only the fine-mesh geometry each block already holds.

## The bug that blocked it: cross-panel neighbour search stopped at the root

`FindNeighborCrossPanel` called `panel_trees[nb]->FindNeighbor(lloc_other, 0,0,0)`. Called
on that tree's ROOT, FindNeighbor stops at the root node, so on a refined mesh it returned
a COARSE ANCESTOR -- and a seam between two blocks that were BOTH refined looked like a
level boundary. Fixed by descending with `FindMeshBlock(lloc_other)`.
`TransformToPanel` was already level-aware and was NOT the problem.

## THE MPI HANG, and the second cross-panel bug (608d43e1)

SMR ran on one rank and HUNG on two. The hang is a blocking `MPI_Wait` in
**`ClearRecv`**, which runs BEFORE `RecvU` in `InitBoundaryValuesAndPrimitives` -- so
instrumenting `RecvAndUnpackCC` shows nothing. **Instrument ClearRecv** for this class.

Cause: cross-panel EDGE and CORNER slots are **not reciprocal**, in general. AthenaK
already skips one family (the cube vertex, `IsCubeVertexCorner`); the generic diagonal
was left alone because with ONE MeshBlock per panel almost every cross-panel diagonal IS
a vertex. Refining a panel makes many ordinary ones appear -- same level AND mixed level
-- and they hang. New `IsSkippedPanelDiagonal` skips every cross-panel edge/corner slot,
send and receive together, **gated on `multilevel`** so the validated unrefined path is
bit-identical.

I first tried skipping only MIXED-LEVEL diagonals; that was not enough (the same-level
ones hang too). Do not repeat that half-step.

## The seam check is on FACE buffers ONLY -- and that is deliberate

An edge/corner neighbour can be across a seam AND across a level boundary even when every
level boundary is radial: the block at the outer edge of a refined shell meets, diagonally,
the coarse block of the next panel. Refusing those rules out radial refinement altogether.
Justified because a DIMENSIONALLY SPLIT sweep reconstructs along one axis at a time and
never reads those diagonal blocks -- the same ground the halo already stands on
([[cubed-sphere-seam-interp]]).

**It is the SPLITTING that makes this safe, not the order of reconstruction** -- PPM and
WENO still read one axis, just a longer line, so "a higher-order scheme would need them"
is WRONG and was corrected in 22dcaf94. What genuinely reads diagonally is the **corner
EMF of the CT update** (O(1) garbage there once cost 1.3% of the magnetic energy on cycle
1 -- see FillPanelCornersCC) and the **cross-derivatives of viscosity and conduction**.
All three are startup FATALs on the cubed sphere today, so nothing that survives the
guards looks at these cells -- but **whoever lifts one of those guards must deal with
this**, and the failure mode is a run that completes looking plausible.

## MHD REFINEMENT: attempted, THREE PIECES IN, STILL FAILS (fbc1f6f5, 68ee6663, 373d589c)

All three curvilinear operators are written and committed as groundwork. They are NOT
enough -- the guard stays up.

| piece | commit | what |
|---|---|---|
| RestrictFC | fbc1f6f5 | area-weighted, normalised by sum of fine areas |
| EMF restriction | 68ee6663 | LENGTH-weighted (CT consumes dxedge*E); keyed on `!same_level` so the seam path is untouched |
| FC prolongation | 373d589c | both halves on the face FLUX A*B; Toth & Roe rewritten with flux accessors |
| FillCoarseInBndryFC | ecb178e8 | a SECOND restriction path, also area-weighted |

**Measured with the guard temporarily lifted** (iprob=8, a static uniform field that
must not evolve, inner radial shell refined):

    nlim=0   ME = 4.888, 4.411, 5.365   <- setup CORRECT, matches unrefined
    nlim=1   ME = 2.8e3, 1.4e3, 2.2e3   <- ONE cycle destroys it
    nlim=2   NaN

The three pieces together improved cycle 1 by a factor of a few and changed nothing
qualitative. **Cartesian MHD+SMR control passes** (`inputs/tests/cpaw3d_smr.athinput`,
120 blocks) with all of this in, so the machinery is sound and the failure is
cubed-sphere-specific.

**WHERE IT FAILS (localised).** The per-block energy check on iprob=8 alternates:

    1.75  5.5e5  2.37  1.3e4  1.80  1.7e3  2.37  4.0e2   (max|dp/p|, panel 0)

The catastrophic blocks are exactly those with **lx1 = 1** -- the fine blocks whose +x1
neighbour is the COARSE block. **The failure is at the RADIAL coarse/fine interface, not
at a seam.**

**HYPOTHESES RULED OUT -- do not repeat these:**

* **"The coarse buffer does not get the seam transform."** The GAP IS REAL: the `do_cs`
  branch of `PackAndSendCC` has an EMPTY `if (neighbor at coarser level)` case where the
  normal branch packs from `ca`. But it is **NOT REACHED** -- a host audit counts ZERO
  cross-panel buffers at coarser or finer level surviving `IsCubeVertexCorner` and
  `IsSkippedPanelDiagonal`, because seam faces at mixed level are refused at startup and
  the diagonals are skipped. **Latent bug, worth fixing if seam-crossing refinement is
  ever attempted, but not the cause here.**
* Areas zero in ghost zones (would make the `/A` blow up): no, `CoordGnomonicEquiangle`
  covers the full ncells range including ghosts.
* All FIVE weighting fixes: RestrictFC, the EMF restriction, both halves of FC
  prolongation, and `FillCoarseInBndryFC` (a SECOND, independent restriction path in
  `prolongation.cpp` -- easy to miss). The last moved cycle 1 from 2.79161e3 to
  2.79154e3, i.e. an O(h) change where an O(1) one is needed.
* `IsSkippedPanelDiagonal` applied to the FC halo.

**Cartesian MHD+SMR control PASSES** with all of this in, so the machinery is sound and
what remains is cubed-sphere-specific at a radial level boundary. Next suspect worth
examining: b.x1f/x2f/x3f are face-normal projections on a NON-ORTHOGONAL basis whose
orientation varies with (xi,eta), so interpolating those components between coarse and
fine locations may need basis handling even WITHIN a panel.

## Still refused, each with its own startup message

* **MHD + refinement**: `RestrictFC` and the Toth & Roe divergence-preserving FC
  prolongation are derived for a uniform Cartesian mesh; on a gnomonic grid the constraint
  is `sum(area*B)`, not `sum(B)`. This is real numerics, not plumbing -- the next piece of
  work if MHD refinement is wanted.
* **Adaptive**: refinement flags are walked through a single MeshBlockTree; a cubed sphere
  has one tree per panel. The STATIC path already loops over `panel_trees`, so the tree
  side is further along than the old guard text suggested.
* **A level boundary lying on a seam**: prolongation/restriction have no seam path.

## Measured under MPI and on GPU

| | |
|---|---|
| CPU MPI 1/2/3/6 ranks | all run; global mass/energy agree to ~1e-13 |
| conservation, 2 and 6 ranks | rel mass drift 2.6e-14 and 1.5e-15 |
| GPU 1 and 2 MI300A | run; 0 of 8 repeats abort; 1-vs-2-rank sums agree to 3.5e-15 |
| GPU vs CPU | 1.4e-13 |
| unrefined cubed sphere | unchanged, serial bit-identical and on 1/2/4 ranks |

`CheckCubedSphereRefinement` also needed `sync_host()`, not `sync<HostMemSpace>()`,
which does not compile for GPU -- the same trap as [[gpu-this-capture-device-lambda]].

## Gotcha when writing a refined_region for a cubed sphere

The region must lie inside the root mesh, so x2/x3 must be exactly [-1,1]; the index
loops then cover the full tangential extent correctly. `<refined_region>` refines the same
logical location on EVERY panel, which is what makes a radial shell come out as a sphere.

## THE GEOMETRY WAS NOT ADDITIVE ACROSS LEVELS -- FIXED, 151da303

Found 2026-08-30 by asking whether the coarse/fine interpolation needed the
non-orthogonal basis. It did, but not as a basis transform: the metric factor is
evaluated at the MIDPOINT.

`area.x1f` and `volume` both carried the angular factor `dth_xi*dth_eta*sin_cell`, an
O(h^2) quadrature of the solid angle that is **not additive**: four fine faces overshoot
their parent by **3.1e-4 relative at nx2=32**, converging as h^2 (1.15e-3 / 3.11e-4 /
8.06e-5 at 16/32/64). Replaced by the EXACT solid angle of the spherical quadrilateral,
`w(xr,yr)-w(xl,yr)-w(xr,yl)+w(xl,yl)` with `w(X,Y)=atan(XY/sqrt(1+X^2+Y^2))`, X=tan(xi).
A corner difference, so it telescopes to round-off (2e-14) and reproduces 4pi/6 exactly.
**Spherical-polar never had this defect** -- `|cosl-cosr|*(phir-phil)` is already additive.

**Why a uniform grid cannot see it.** The radial momentum balance cancels for ANY area
values, because the source term IS `z_ov_rE = (area1r - area.x1f)/volume`. At a level
boundary the flux correction hands the coarse cell `sum(A_fine*F)` while the source still
uses `A_coarse`; the residue `p*(sum(A_fine)-A_coarse)/V` is a spurious radial force. With
the exact solid angle Omega cancels out of `z_ov_rE` altogether.

**MEASURED, hydro iprob=1 static uniform, ONE cycle, radial level boundary:**

| | before | after |
|---|---|---|
| spurious 1-momentum | -1.82e-05 | **-4.0e-18** |
| kinetic energy | 8.61e-11 | **5.3e-32** |

Unrefined control was and stays exactly 0. No uniform-grid regression: MHD L1(B) at
nx2=32 5.8666e-5 -> 5.8921e-5 (0.43%), 16->32 ratio unchanged, seam conservation
5e-15/1e-14, SMR conservation -3.5e-14/-3.7e-14, halo still 2nd order, and seam face
areas now agree ACROSS PANELS to 3.5e-13.

**THE LESSON: conservation is not balance.** Hydro SMR telescoped to round-off from the
first day and that is exactly why this hid -- the bookkeeping was perfect while the
equilibrium was not. Test a STATIC state for zero acceleration, not just for closed
budgets. `iprob=1` + the `.hst` 1-momentum column is a two-minute check needing no code.

**MHD is NOT fixed by this** (nlim=1 still 1e38) -- a separate, larger, FC-side defect at
the radial interface. Guard stays up. But these are now ruled out for it:
* ghost VALUES: seam ghosts 5e-5, radial ghosts 1e-2, both small
* MONOPOLES: ghost div B 2.5e-5 == the initial condition's own 2.1e-5 floor, so
  `ProlongFCInternalCurvi` (Toth & Roe on Phi = A*B) is doing its job
* CC prolongation: iprob=8's `u0` is a global CONSTANT, so it is exact by construction
* the failure is at `i = ie` on EVERY block, spread over all (j,k) -- the radial
  interface, NOT the seams, and it hits blocks whose radial ghosts are exactly clean,
  which points at the interface flux/EMF correction rather than at the halo

New diagnostics in `cs_test.cpp`: radial-ghost error, per-cell div B vs exact, and the
(i,j,k) of the worst energy cell -- the last is what localised this in one run.

## THE MHD DEFECT IS THE SKIPPED DIAGONALS -- FOUND, baf23918 (not yet fixed)

`IsSkippedPanelDiagonal` (`bvals.hpp`) drops EVERY cross-panel edge/corner slot when
`multilevel`, because they are non-reciprocal and hang `ClearRecv`. Its own comment
already warned that **the corner EMF of the CT update is the one thing that reads
DIAGONALLY**. That warning is now measured.

iprob=8, |B| = 1e-2, radial level boundary. The exact EMF is IDENTICALLY ZERO (v = 0
everywhere, ghosts included: u0's momentum is a global zero and prolongation of zero is
zero). Measured:

| where | EMF |
|---|---|
| interior of a block | 1e-10 |
| interior of the interface FACE | 1e-14 |
| level boundary MEETS a tangential block boundary | **1.3e-2** = O(\|B\|) |

Forcing the predicate to `false` (serial, no hang) drops it to **1.2e-8 -- a factor 1e6**,
and the induction is then healthy (interface only 4.1x interior, smooth in b0c). **The fix
is to make those diagonals RECIPROCAL, not to skip them.** Guard stays up.

**RULED OUT, do not re-try:** the FC flux correction (disabling it makes the EMF WORSE,
2.5e-2, and breaks same-level edges too -- it is what makes shared edges single-valued);
prolongated ghost VALUES (1e-4 = 1% of B, exactly linear in b0c); the interface face
`b.x1f` itself (exact, 0.0, both sides); ghost MONOPOLES (div B 2.5e-5 vs the initial
condition's own 2.1e-5 floor); CC prolongation (u0 is a global CONSTANT here).

**A SECOND defect remains**, on the MOMENTUM side not the induction: with the diagonals
restored, dp/p is still ~100x the unrefined baseline (which is itself 2.2e-4 -- a uniform
Cartesian field is NOT exactly preserved on this grid) and grows faster than B^2. Likely
the magnetic analogue of the [[cubed-sphere-smr]] pressure-balance fix 151da303.

## Two measurement traps that cost real time this session

1. **Scanning an EMF component outside its own index range.** `efld.x1e` runs ALONG x1,
   so it is a CELL index in i and a CORNER index in j,k, cyclically for the others.
   Scanning all three to +1 in every index reads memory the update never wrote and prints
   ~1e12 -- indistinguishable from a catastrophic bug. I nearly reported it as one.
2. **`fmax(NaN, x)` returns `x`.** A run that had already gone NaN reported dp/p =
   5.08e-05 and looked HEALTHY, better than the working case. Gates now count non-finite
   cells and say so. Check `grep -c nan` on the log before believing any max.

## DEFECT #2 IS A CONSERVATION FAILURE, not an accuracy one -- 6b21255c

After restoring the diagonals the induction is healthy but the run still dies. Measured
with REFLECTING radial BCs (exact physical flux = 0), iprob=8, b0c=1e-2, nx2=32:

| | mass | energy |
|---|---|---|
| hydro SMR | -3.5e-14 | -3.7e-14 |
| MHD, 2 cycles | +8.0e-08 | +1.3e-07 |
| MHD, 5 cycles | +1.7e-07 | +2.9e-07 |

**Seven orders of magnitude, through the SAME cell-centred flux correction**, secular, and
INDEPENDENT of defect #1 (+1.7e-7 skip on vs +2.1e-7 skip off). No floors fire. It is an
INSTABILITY: nx2=16 survives to t=0.02 (dp/p 1.1e-4), nx2=32 goes non-finite in <10 cycles
-- worse with resolution. **Mass is the sharp clue**: the mass flux is rho*v.n with no
magnetic term, and hydro telescopes to round-off with real velocities through that same
correction. Next probe: instrument level-boundary flux telescoping the way
`CSTestSeamFluxCheck` already does for seam faces.

**RULED OUT with numbers:**
* **FC prolongation.** Ghost error is 2nd order in **L1** (6.62e-4 / 1.66e-4 / 4.17e-5 at
  16/32/64, ratios 3.98, 3.99). Linf is 1st order (2.01, 1.98) but that is a **limiter
  clipping at an extremum** -- ordinary. **Never read the Linf number alone**; I nearly
  reported "the prolongation is first order" as the defect.
* Floors (event counters empty). The b0c=0 control is **WEAK** -- v stays identically zero
  so the mass flux is zero and the correction is never exercised.
* The unrefined baseline as an explanation: MEASURED, scales exactly as b0c^2 (2.456e-10,
  2.456e-8, 2.465e-6, 2.443e-4), so refined really is ~70x it, rising to ~1900x at b0c=0.1.

## TELESCOPING INSTRUMENTED -- 07a680d6. **THE COARSE/FINE FACE IS EXACT.**

New **`CSTestLevelFluxCheck`** in cs_test.cpp: pairs every block-boundary face BY GEOMETRY
in ALL THREE directions, buckets by (panel, direction, position at the COARSEST block's
width) so a coarse face and its 4 children land together, sums signed area*flux. Zero =
telescopes. Independent of the buffer machinery it audits.

**It OVERTURNED the hypothesis.** iprob=8, b0c=1e-2, reflecting walls:

| | same-level faces | coarse/fine faces |
|---|---|---|
| unrefined hydro | **0.0 exactly** | -- |
| unrefined MHD | **0.0 exactly** | -- |
| refined hydro | 2.7e-19 (round-off) | 1.1e-22 (round-off) |
| refined MHD | **1.0e-09 = 14% of the flux** | **1.7e-24 EXACT** |

**The coarse/fine interface is NOT the problem -- it telescopes to round-off.** What
breaks is the ORDINARY SAME-LEVEL block boundary, and ONLY once refinement is on. The user
framed it right: mass flux is rho*v.n with no B in it, so something hydro does is not being
done in MHD. Same-level faces telescope for a DIFFERENT reason than corrected ones -- the
ghosts are exact COPIES, so both sides reconstruct the same state and the Riemann solves
agree bit for bit. **So something about a multilevel run stops the MHD ghosts being exact
copies at same-level boundaries.** That is the next thing to chase: diff what
`multilevel` changes in the FC/CC same-level ghost exchange (bvals_fc / bvals_cc), since
that is the only switch that differs between the two working rows and the broken one.

NOT the skipped diagonals: disabling them makes it WORSE (1.8e-5 against a 1.1e-5 scale =
170%, vs 14%). **Defect #2 is independent of defect #1.**

**Two gate traps, both caught by the HYDRO NULL CONTROL** (always validate the instrument
on a case known to conserve):
1. **The DIRECTION must be in the bucket key.** Without it x1 and x2 buckets collide and
   pool unrelated faces -- unrefined hydro then reported a 100% mismatch.
2. **A bucket needs BOTH signs** to be a shared face. One-sided buckets are PHYSICAL
   domain boundaries where nothing should cancel; counting them reports the wall's own
   flux as a telescoping failure (it showed a fake -4.6e-5 "same-level" source).


## DEFECT #2 IS CLOSED -- 23201e2c. TWO multilevel-only boundary defects, both fixed.

**The control that cracked it: a mesh with `refinement = static` but NO level boundary
anywhere.** 48 blocks all at one level, `multilevel` true. MHD already fails there, hydro
does not. So the defect never needed a coarse/fine interface -- only the `multilevel`
switch -- which reduced the search to what that flag changes in the SAME-LEVEL exchange.

| config (iprob=8, b0c=1e-2, scale 2.8e-9) | x1 radial | x2/x3 tangential |
|---|---|---|
| multilevel FALSE, identical 48-block mesh | 0.0 exact | 0.0 exact |
| hydro, multilevel TRUE | 0.0 exact | 0.0 exact |
| MHD, multilevel TRUE (as committed) | 7.3e-10 | 4.6e-10 |
| MHD, diagonals restored only | **4.5e-16** | 5.5e-9 WORSE |
| MHD, expansion dropped only | 7.3e-10 | **0.0** |
| MHD, BOTH | **0.0** | **0.0** |

Perfectly separable: one defect breaks the RADIAL faces, the other the TANGENTIAL.

**Defect 2a -- the FC buffer expansion.** `buffs_fc.cpp` widens the SAME-LEVEL edge and
corner buffers by one shared face when `multilevel` ("always include the overlapping
faces"). That layer is REDUNDANT even in Cartesian: the same (i,j,k) is already written
by the FACE buffer of the block that owns it, and which source wins is only the order of
the unpacks. Cartesian values agree; across a seam they do not (the face buffer applies
the basis transform and the along-seam resample, the diagonal is a plain copy from
another chart) and the diagonal CLOBBERS the correct value. Skipped for use_cubed_sphere.

**Defect 2b -- the blanket diagonal skip (= defect #1).** Auditing the non-reciprocal
slots ONE AT A TIME shows they are only ever (i) CUBE VERTICES or (ii) MIXED-LEVEL
cross-panel diagonals. **Every same-level cross-panel diagonal pairs correctly.** So skip
only the mixed-level ones. The earlier "skip only mixed level was not enough" is
**EXPLAINED, not contradicted**: `IsCubeVertexCorner` covered the x2x3 EDGES [40..47] but
NOT the 3-D CORNERS [48..55], a family that cannot exist until the mesh splits the RADIAL
direction (with one block spanning x1 a corner slot has no neighbour at all). Extending
it to 48..55 is what makes the half-step work.

The gate for MPI safety is not "is the table reciprocal" but **"is any non-reciprocal
slot still EXCHANGED"** -- the audit now prints that count (0 in both meshes).

**Measured after the fix:** cubed_sphere_{smr,mhd_conv,uniform,rigidrot} BIT-IDENTICAL to
HEAD; hydro SMR under MPI at 2/3/6 ranks matches serial to 1e-14, no hang; the multilevel
MHD cubed sphere with no level boundary runs 200 cycles clean (dp/p 1.2e-5, dt steady).

**WHAT IS LEFT.** With a REAL coarse/fine interface the telescoping is now round-off
(2e-25, both families, was 1.0e-9 = 14%) but the run still goes NON-FINITE AT CYCLE 43
(clean at 42, dp/p 1.7e-5 and growing slowly), and MPI still hangs at 2/3/6 ranks. That
is the FACE-CENTRED prolongation/restriction numerics at the radial interface -- the
"real numerics" piece that was always the stated remainder. Guard stays up.

**`mesh/cs_dev_allow_mhd_refinement = true`** now lifts the startup FATAL for diagnostics
(development only; documented at the guard). No more patching mesh.cpp each session.

Working inputs live in the session scratchpad: `lvl_mhd.athinput` (multilevel, NO level
boundary -- the control) and `ref_mhd.athinput` (real radial level boundary, 54 blocks);
the difference is only `<meshblock>/nx1` 16 vs 8.


## THE FC SIDE OF THE RADIAL INTERFACE -- ae7d4f2b + 9db0cbdc. RUN NOW SURVIVES.

After 23201e2c the flux telescoped at round-off yet the refined MHD run still went
non-finite at CYCLE 25. **It was not the prolongation arithmetic. It was ghost cells that
nothing ever wrote.**

**Localise by splitting the interface FACE, not by looking at it as a whole:**

| coarse/fine interface EMF (exact = 0) | face INTERIOR | face EDGE RING |
|---|---|---|
| fine side | 5.3e-9 | **3.3e-2** |
| coarse side | 3.6e-10 | **2.4e-2** |

The edge ring is where the level boundary meets a tangential block boundary, and the EMF
there reads the DIAGONAL ghost block beyond the interface. Those cells were **EXACTLY
ZERO** -- 1050 per block -- and the error was **FLAT in resolution** (9.999e-3, 9.996e-3,
9.989e-3 at nx2 = 16/32/64) at |B| = 1e-2. *A flat error equal to the field's own
magnitude is a cell that is never written, not one interpolated badly.* The density in
the same cells was FLOORED to 1.2e-38.

**Two reasons, both fixed:**
1. `IsSkippedPanelDiagonal` drops cross-panel diagonals at another level; nothing else
   wrote them. `FillPanelCorners{CC,FC}` does NOT -- it fills the TANGENTIAL corner, and
   these regions have one tangential index ACTIVE. New
   **`FillSeamLevelDiagonals{CC,FC}`** fill them from the two flanking halos.
2. **`Prolongate{CC,FC}` prolongated those slots anyway, from a coarse buffer nothing had
   filled, writing ZEROS over everything.** They now skip a slot whose exchange was
   skipped, and the corner fills run AGAIN after prolongation (their stencil reads the
   radial ghost strip prolongation has only just written).
3. `flux_correct_fc.cpp` skipped only cube vertices, so it still SENT EMFs on the skipped
   diagonals -- harmless in serial, **that was the MPI hang**.

**Measured (iprob=8, b0c=1e-2, USER radial BCs):** edge-ring EMF 3.3e-2 -> 2.0e-10;
600 cycles with dt unchanged to 7 digits, max|dp/p| 2.4e-8, max|v| 4.2e-7, mass drift
-3.2e-9; **div B EXACTLY preserved** (1.59511e-06 at cycle 1 and at cycle 50); MPI 2/3/6
ranks run and agree with serial to 1e-14. All four cs regression tests bit-identical.

**USE `user` RADIAL BCs FOR ANY MHD GATE.** `reflect` is not exact for B_r and
manufactures an O(|B|) error at the wall (1.07e-2 radial-ghost error, and a 7.5e-6 mass
drift that is pure wall artifact -- with b0c=0 the same run drifts EXACTLY 0). It
confounded a whole round of measurements.

## THE DIAGONALS ARE NOW DONE PROPERLY -- 0626f1c9. The blanket skip is GONE.

**The search asymmetry was the CHILD SELECTION.** `GetLeaf` is indexed in the
NEIGHBOUR's axes, but the octant indices passed to it were built from the offsets in
OURS. A seam can send our x2 to their x3, so a coarse block enumerating the two fine
blocks along a shared edge **walked the wrong axis**: panel 5 listed panel 1's two blocks
stacked in x3 when its real neighbours were the two stacked in x2. Each side then named a
block that did not name it back (`reverse slot = -1`), which is why no `dest` fix could
help -- **do not chase `dest` for this class again; check the child selection first.**

New **`PanelEdgeMap`** (mesh.hpp) carries the three things a diagonal needs across a
seam: the offsets pointing back in THEIR axes, WHICH of their axes runs along the edge,
and whether that axis is REVERSED (which swaps the two subblock halves). Used in
meshblock.cpp for the child selection, the subblock index AND `dest`, in all three edge
families and the corners. Identity when the panels share a chart.

    x1x2 edges 16 non-reciprocal -> 0 ; x3x1 edges 32 -> 0 ; corners 80 -> 64

Every remaining non-reciprocal slot is a CUBE VERTEX (three panels, no fourth block --
geometry, not bookkeeping). So **`IsSkippedPanelDiagonal` is DELETED**, with all its call
sites, and so is `FillSeamLevelDiagonals{CC,FC}` -- the extrapolation only existed to
paper over the skipped exchange. `IsCubeVertexCorner` alone remains. Audited clean in
THREE meshes: 54 blocks, 108 blocks, and a 444-block THREE-level nested one.

Exchanging beats extrapolating, and it also repairs the FACE prolongation (its stencil
reads the coarse buffer's edge cells): ghosts on the face 1.2e-4 -> 7.5e-6, 2-D diagonal
8.2e-4 -> 8.0e-5.

**Verified:** all four cs tests BIT-IDENTICAL to 07a680d6; hydro SMR bit-identical and
still exact under MPI 2/3/6; refined MHD 400 cycles, dt unchanged to 7 digits, div B
preserved to 6 digits, mass drift -2.5e-9, MPI 2/3/6 agrees with serial to 5e-14.

## CONVERGENCE MEASURED, 41ab15cb -- REFINED MHD DOES NOT CONVERGE. Guard stays up.

Rebuilt from `inputs/tests/cubed_sphere_mhd_conv.athinput` VERBATIM (iprob=9, b0c=1,
tlim=1.0, nx1=8, user radial BCs), scaling only nx2=nx3; the refined twin adds
`<meshblock>/nx1 = 4` (2 root radial blocks) plus a level-1 `<refined_region1>` over
x1 in [1.0,1.5], so the level boundary is the sphere r = 1.5. Inputs are in the
scratchpad as `conv2/{ctl,ref}_{16,32,64}.athinput`.

| nx2 | 16 | 32 | 64 | ratios |
|---|---|---|---|---|
| L1(B) unrefined CONTROL | 5.278e-4 | 2.203e-4 | 1.264e-4 | 2.40, 1.74 |
| L1(B) REFINED | 3.984e-3 | 3.146e-3 | -- | **1.27** |
| Linf(B) refined | **0.60** | **1.27** | -- | GROWS |

**The refined case is 7.5x the control and does not converge, and its Linf is O(|B|=1).**
That is the answer: refined MHD is not usable, and the guard is correct.

**LOCALISED, and it is not where the last round looked.** The new Linf line reports the
block's level and its x1 neighbours' levels. EVERY refined run, at t = 0.05 already, puts
the worst cell on the **COARSE block at i = 0** -- its INNER radial face, i.e. the
coarse/fine interface -- with its -x1 neighbour one level finer, and at a tangential
block edge (j or k = 0). **So the defect is on the COARSE SIDE of the interface, in the
first coarse cell.** The fine-side ghosts are fine (that was ae7d4f2b + 0626f1c9).

**And it needs v != 0 to appear**: the static-field gate (iprob=8) is clean at 1e-9 with
the same mesh. So the next suspect is what the FLUX/EMF CORRECTION writes into the coarse
cell -- the length-weighted EMF restriction (68ee6663) and `AverageBoundaryFluxes` /
`ZeroFluxesAtBoundaryWithFiner` on a gnomonic grid -- not prolongation, not the ghosts.

**A TRAP IN MY OWN SWEEP, fix it before quoting the control again:** the control's ratios
are 2.40 then 1.74, i.e. it is NOT second order either, because the sweep scales only
nx2=nx3 and leaves **nx1 = 8 fixed**, so the RADIAL truncation stops falling and takes
over by nx2 = 64. Scale nx1 with the tangential resolution before judging any of this.

## STILL OPEN (the guard stays up)

1. **CONVERGENCE with refinement: MEASURED and FAILING, see above.** iprob=9 with a radial level
   boundary: L1(B) 7.35e-4 -> 5.09e-4 over nx2 16->32 (ratio 1.44). But the UNREFINED
   control on the same input is 2.99e-4 -> 1.61e-4 (ratio 1.86), i.e. **the control is
   not clean 2nd order either**, so the INPUT is the first suspect, not the scheme --
   rebuild it from `inputs/tests/cubed_sphere_mhd_conv.athinput` verbatim (one radial
   block, its exact bc_* settings) and only then judge the interface. Scratchpad has
   `conv_ref*.athinput` / `conv_unref_*.athinput`.
2. Exchanged ghosts are NOISIER than the old extrapolation on the smooth static field
   (spurious KE 2.2e-11 vs 5.1e-14 at t=1.5). Expected -- the exchange carries the coarse
   cell's real interpolation error while extrapolation never sees the coarse grid -- and
   the exchange is still the correct general answer. Worth re-checking on a test with
   structure.
3. 18 three-dimensional corner ghost cells at cube vertices are still exactly zero. Not
   read by an edge-centred EMF; the interface EMF is at round-off with them left alone.

## HISTORICAL: what the diagonals looked like BEFORE 0626f1c9

* SAME-LEVEL cross-panel diagonals: **YES**, exchanged, audited reciprocal.
* MIXED-LEVEL cross-panel diagonals: **NO.** They are skipped and filled by local
  quadratic extrapolation. The fill is 2nd order and converges at the SAME order as the
  prolongated face ghosts beside it, so it is not the limiting error -- but it is not the
  neighbour's data, so the EMF on those edges is not exactly single-valued.
* 3-D CORNERS: 54 cells still at zero. Left deliberately; an edge-centred EMF reads
  2-D diagonals only, and the interface EMF is at round-off with them untouched.

**Two concrete obstacles to doing it properly, both MEASURED, do not re-derive:**
1. `dest` for a MIXED-LEVEL edge is computed with plain `NeighborIndex` and **no panel
   map** (`meshblock.cpp`, the "neighbor at coarser level" branches). Adding it (plus a
   subblock-parity remap) moves `dest` onto the right EDGE FAMILY but the subblock parity
   is case-dependent: the naive rule and its inverse each fix some pairs and break
   others. Attempted and REVERTED -- inert today because those slots are skipped.
2. Half the failures are a **SEARCH asymmetry**: the partner does not list the block in
   ANY slot (`reverse slot = -1` in the audit). No `dest` fix can repair those; the
   cross-panel diagonal search at mixed level has to enumerate both fine blocks covering
   a coarse block's diagonal region on the other panel.

## REFINED MHD NOW CONVERGES -- the defect was the EDGE FLUX BUFFERS, not the interface

Found 2026-08-30. Two gaps, both of the same shape: **a cubed-sphere code path that was
written when the radial direction was never split, and is only REACHED once a radial
level boundary exists.** Neither is about prolongation, restriction or the flux
correction arithmetic -- all of which the last four sessions had already measured clean.

### The dominant one: x1x2 and x3x1 EDGE flux buffers had NO seam transform

`flux_correct_fc.cpp` gives the x1face, x2face and x3face branches the full panel-seam
treatment (reverse the along-seam index, flip the SIGN of the along-seam EMF, write it
into the slot the neighbour reads, `vout`). The **x1x2 edge** branch (carries x3e) and
the **x3x1 edge** branch (carries x2e) had *none of it* -- not even the index reversal.
Its own header comment said why: *"x1x2 and x3x1 edge buffers never exist here -- the
radial direction has physical boundaries at both ends, so a block has no x1-direction
neighbour."* **True until radial refinement, false after.** Those slots then appear, and
they land exactly on the ring where the level boundary meets a seam -- which is precisely
where every localisation had been pointing.

**MEASURED (iprob=9, nx1 scaled with nx2, refined twin = 2 radial root blocks + a level-1
region over x1 in [1.0,1.5]):**

| nx2 | 16 | 32 | 64 | ratio |
|---|---|---|---|---|
| L1(B) unrefined control | 5.278e-4 | 2.067e-4 | 9.386e-5 | 2.55, 2.20 |
| L1(B) refined, BEFORE | 3.984e-3 | 3.146e-3 | -- | **1.27** |
| L1(B) refined, AFTER | **2.931e-4** | **1.077e-4** | -- | **2.72** |
| Linf(B) refined, before -> after | 0.60 -> 3.2e-3 | 1.27 -> 1.6e-3 | | |

L1(B) falls 12x, Linf 189x, and the worst cell moves OFF the interface into a block
interior. The refined run is now BETTER than the control and converges at the control's
own rate (the control is not clean 2nd order either -- that is the separate, closed
radial-BC phase lag, [[cubed-sphere-mhd-convergence]]).

### The second one: cross-panel FINE->COARSE buffers were never packed at all

The `do_cs` branch of `PackAndSendCC` (and of `PackAndSendFC`) had an **empty**
`if (neighbor at coarser level)` case -- the one place the normal branch reads `ca`/`cb`.
An earlier session found this gap and correctly ruled it out as *unreachable*, because
mixed-level cross-panel diagonals were skipped. **0626f1c9 deleted that skip, which made
it reachable.** A never-written Kokkos buffer is ZERO, so the coarse block read a ghost
state of exactly zero across the seam.

That is why it hid from every gate: with v = 0 (iprob=8) zero IS the right momentum, and
hydro never reconstructs through a corner ghost -- only the CT corner EMF reads
diagonally. **Worth 10% of L1(B) (3.984e-3 -> 3.582e-3); real, but not the story.**

### The control that made this findable in ten minutes

**Run the SAME level boundary with HYDRO** (iprob=3, which is iprob=9's own hydro state,
and `problem/conv_errors=1` reports L1(v), L1(p) for it). Refined hydro came out at
1.729e-4 against a 3.396e-4 control -- BETTER, converging, and bit-identical before and
after both fixes. One run split the search space in half: everything shared with hydro
(geometry, CC flux correction, restriction, prolongation of u0, the seam momentum
transform) was exonerated, leaving only the FC/EMF side.

### Also fixed while in there
* the stale "edge buffers never exist" comment
* the along-seam resample in bvals_{cc,fc} now takes its cell COUNT and index origin from
  the coarse arrays when the source is `ca`/`cb` (cnx2/cnx3, cjs/cks), not always nx2/nx3

**Verified:** all SIX `inputs/tests/cubed_sphere_*.athinput` bit-identical to 41ab15cb
(compare the `.hst`, and run each in its OWN directory -- `cubed_sphere_smr` and
`cubed_sphere_rigidrot` share the basename `cs_rigidrot` and will interleave one `.hst`
if run concurrently, which looks exactly like a regression).

## THE "RANK DEPENDENCE" WAS THE GATE ITSELF -- RETRACTED, and FIXED in 2943fd02

I reported, in 74cbc8df's own commit message, that a refined cubed-sphere run is
rank-dependent. **That is WRONG and is retracted.** There is no rank dependence.

`CSTestConvErrors` and `CSTestConsSums` loop over `pmbp->nmb_thispack` -- the MeshBlocks
in THIS RANK's pack -- and printed the totals as if global. L1 divides by an `ncell`
counted the same way, so the number is the mean over whatever share of the mesh the rank
owns; Linf is that rank's max; and EVERY rank appended its own line to
`<basename>-cs-errs.dat`, so which survived was a race. One identical solution therefore
printed a different L1 at every rank count:

| | np=1 | np=2 | np=3 | np=6 |
|---|---|---|---|---|
| refined MHD L1(B), rank-local | 2.93e-4 | 3.37e-4 | 3.23e-4 | 2.59e-4 |
| refined hydro L1(v), rank-local | 1.73e-4 | 2.06e-4 | 1.40e-4 | -- |
| refined MHD L1(B), REDUCED | 2.93118652e-4 | same | same | same |

**THE MEASUREMENT THAT SETTLED IT: compare BINARY DUMPS, not norms.** `file_type = bin`
every 0.05 in time, then `cmp` the files. All 21 dumps of the full refined MHD run are
bitwise identical at 1, 2, 3 and 6 ranks, and all 21 of the refined hydro run likewise.
That is a direct statement about the SOLUTION and it cannot be confounded by a
diagnostic. A norm cannot make that statement; it is one number, and a wrong reduction
in it is indistinguishable from a wrong solution.

**THE LESSON, and it is the same one as `fmax(NaN,x)` and the one-sided-bucket gate:
validate the instrument before believing what it says about the code.** A cheap check
here would have been to run the UNREFINED case at 2 ranks -- it would have shown the same
spread, on a configuration already known to be bit-identical, and named the gate
immediately. Cost: several hours, a wrong finding in a commit message, and a wrong entry
in this file.

Fixed in **2943fd02**: both gates now `MPI_Allreduce` (SUM for the sums and counts,
MAXLOC for Linf plus an `MPI_Bcast` of the owning rank's formatted location line, since
block ids are rank-local) and report on rank 0 only. `CSTestSeamFluxCheck` and
`CSTestLevelFluxCheck` pair faces BY GEOMETRY across blocks and cannot be reduced without
communicating the pairs; both already say they are serial-only, and that stands.

**So MHD + refinement now passes MPI too**, at 1/2/3/6 ranks, bitwise.

## GPU: PASSES. Job 11240524 on apudev (2x MI300A, vipa1001), 2026-08-30.

`bench/cs_smr_gpu/` holds the inputs and `submit_cs_gpu.sh`. Build: the standard
[[viper-hip-build-recipe]] plus `-D PROBLEM=cs_test`, into `build_cs_gpu`. No warnings in
bvals or cs_test (the run was a FRESH build, so warnings were not hidden -- see
[[gpu-this-capture-device-lambda]]). `xnack+` confirmed in `rocminfo`, so this is the real
unified-memory APU path.

| run | GPU | CPU serial |
|---|---|---|
| refined MHD, 1 GPU, L1(B) | 2.931187e-04 | 2.93118652e-04 |
| refined MHD, 1 GPU, Linf(B) | 3.158684e-03 | 3.15868386e-03 |
| refined MHD, 2 GPUs | **identical to 1 GPU in every digit** | |
| unrefined MHD control | 5.277894e-04 | 5.277894e-04 |
| refined HYDRO, L1(v) | 1.728746e-04 | 1.728746e-04 |

Every figure agrees with CPU to all seven printed digits, and **1 vs 2 GPU is BITWISE
identical over all 5 binary dumps**, which is the same gate that settled the CPU MPI
question -- compare DUMPS, not norms.

**apudev caps a job at 15 minutes** (`QOSMaxWallDurationPerJobLimit`); 00:29:00 is
rejected outright at submit time.

## THE FULL EVIDENCE FOR MHD + REFINEMENT, in one place

Measured at 2943fd02, all on the iprob=9 rigid rotation with a radial level boundary:

* **converges**: L1(B) 2.931e-4 -> 1.077e-4 over nx2 16->32 (2.72) vs the unrefined
  control's 5.278e-4 -> 2.067e-4 (2.55); refined is BELOW the control at both
* **level-boundary flux telescopes** to round-off (1.7e-21 against a 2.9e-5 scale)
* **div B** preserved; **conservation** at round-off
* **MPI** bitwise identical at 1/2/3/6 ranks, 21/21 dumps, MHD and hydro
* **GPU** matches CPU to 7 digits, 1 vs 2 GPUs bitwise, 5/5 dumps
* all six `inputs/tests/cubed_sphere_*.athinput` bit-identical to 41ab15cb

Not measured: AMR (still a separate startup FATAL and a separate problem), a level
boundary lying ON a seam (still refused), fofc/viscosity/conduction (still refused).

## THE EDGE-BUFFER FIX IS NOT AN SMR FIX. It repairs any RADIALLY SPLIT cubed sphere.

Asked directly ("if I don't use SMR, none of this runs, right?") and MEASURED. The
x1x2/x3x1 edge flux buffers exist whenever a block has an x1-direction NEIGHBOUR. Static
refinement forces that, but so does plain `<meshblock>/nx1 < <mesh>/nx1` with
`refinement = none`. iprob=9, nx2=16, NO refinement anywhere, only the radial block count
changed:

| | 1 radial block | 2 radial blocks |
|---|---|---|
| L1(B) at 41ab15cb | 5.277894e-04 | **1.291796e-03** |
| Linf(B) at 41ab15cb | 3.457e-03 | **1.263e-01** |
| L1(B) at 74cbc8df | 5.277894e-04 | **5.277894e-04** |
| Linf(B) at 74cbc8df | 3.457e-03 | **3.457e-03** |

2.4x in L1 and 36x in Linf, and now the answer is DECOMPOSITION-INDEPENDENT to every
printed digit -- which is the real gate: a correct scheme must not care how the mesh is
cut into blocks.

**Why no regression test caught it: all six `cubed_sphere_*.athinput` use ONE radial
block** (`<meshblock>/nx1` == `<mesh>/nx1`). That is also why they stayed bit-identical
through this whole change. **Any new cs regression test should split x1**, and
"same answer under a different block decomposition" is the cheapest strong gate available
on this grid -- cheaper than a convergence sweep and sensitive to exactly the class of
bug that has dominated this thread.

The other pieces ARE refinement-only: the cross-panel fine->coarse pack needs a neighbour
at a different level, and the guard is refinement by definition.

## HOW TO PICK THIS UP (build dirs, inputs, the reproducer)

Build dirs in the repo root, all untracked and gitignored, all at 790d02c0:
`build_cs` (serial), `build_cs_mpi` (MPI), `build_cs_gpu` (HIP/MI300A). Each is
`cmake -S . -B <dir> -D PROBLEM=cs_test -D CMAKE_BUILD_TYPE=Release` plus, for MPI,
`-D Athena_ENABLE_MPI=ON` (needs `module load gcc/14 openmpi/5.0` -- openmpi is INVISIBLE
until a compiler is loaded), and for GPU the [[viper-hip-build-recipe]] flags.

GPU submission lives in `bench/cs_smr_gpu/submit_cs_gpu.sh`; **apudev caps a job at 15
minutes** and rejects 00:29:00 at submit time with `QOSMaxWallDurationPerJobLimit`.

The convergence sweep is generated, not stored: take
`inputs/tests/cubed_sphere_mhd_smr.athinput` (refined) and
`inputs/tests/cubed_sphere_mhd_conv.athinput` (control) and scale nx1 = nx2/2 with
nx2 = nx3, `<meshblock>/nx1` = nx1/2 for the refined twin. **Scale nx1 with nx2** -- a
fixed nx1 caps the CONTROL at ratio 1.74 and makes everything look broken.

When comparing the six cs regression tests, run each in its OWN directory:
`cubed_sphere_smr` and `cubed_sphere_rigidrot` share the basename `cs_rigidrot` and will
interleave one `.hst` if run concurrently, which looks exactly like a regression.
