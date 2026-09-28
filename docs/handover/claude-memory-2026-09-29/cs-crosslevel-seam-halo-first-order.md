---
name: cs-crosslevel-seam-halo-first-order
description: ROOT-CAUSED and ESSENTIALLY CLOSED by TWO fixes (a01ace75 + 180a9b3e) -- a level boundary no longer meaningfully degrades the cubed-sphere seam halo -- the along-seam resample was silently skipped on cross-level seam EDGE buffers because cs_seam was classified from buffer EXTENTS, which face-centred staggering breaks at a coarse/fine boundary; 7.8x better, but the doubly-ghost x2x3-edge/corner buffers still take no resample and are now the leading first-order term
metadata:
  type: project
---

**ROOT-CAUSED AND PARTLY FIXED, a01ace75.** Found by the GHOST SCAN added in 397b4ad3.
Read THE ANSWER first; the long ruled-out list below is kept because it is what made the
final measurement obvious, and because the remaining term needs it.

## THE CORNER EXTRAPOLATION IS NOW AT ITS PRACTICAL LIMIT (e5e8983a)

`FillPanelCornersFC` averaged its two quadratic extrapolations EQUALLY although they do
not reach equally far (dj cells along j, dk along k). Weighting each by the OTHER's
distance, `w = dj/(dj+dk)` on the k-extrapolation, favours the shorter reach:

    tangent ghost halo   2.2562e-03 -> 1.4860e-03  (1.52x)  at BOTH levels
    corner, same level   1.9807e-03 -> 1.3047e-03  (1.52x)
    corner, refined      2.0911e-03 -> unchanged (its max is the DIAGONAL cell, dj = dk)

**QUADRATIC IS THE OPTIMUM**, bracketed both ways on the same-level corner: LINEAR
(amplification 3, 5) loses to truncation at 5.4389e-03; CUBIC (15, 49) loses to
amplification at 2.4572e-03; quadratic (7, 17) sits at 1.9807e-03.

**AN EARLIER MEASUREMENT IN THIS FILE IS NOW REVERSED, and that is the real lesson.**
Before the strips were fixed, linear beat quadratic 3x on the refined case because the
corner was AMPLIFICATION-dominated. With accurate strips it is TRUNCATION-dominated and
linear is decisively worse. **A tuning result is only valid for the state of the code it
was measured in -- re-measure it after any upstream accuracy change.**

FACE-CENTRED ONLY: the same weighting in `FillPanelCornersCC` left both hydro gates
identical and moved mhd_smr marginally worse in L1 / better in Linf, i.e. noise, so it
was not made.

**THE CEILING.** Beating ~1.3e-03 needs a REAL cube-vertex exchange -- sampling the third
panel that meets at the vertex, which does hold that data -- instead of extrapolating from
our own halo. Substantially larger work; NOT started, and the user should be asked first.

## THE THIRD PIECE: the DOUBLY-GHOST buffers, and what was masking them (e4be92c0)

Three coupled changes that only work together -- which is why the first two measured
INERT on their own and were once reverted:

1. **cs_seam now covers slots 40-55.** They are ghost in BOTH tangential directions so
   the slot index alone names no along-seam axis, but the NEIGHBOUR TABLE does: at a
   non-cube-vertex corner exactly ONE flanking face is a seam, and that one fixes the
   seam normal. (Both = a cube vertex, whose exchange is skipped anyway.)
2. **Resample stencil bounds from the SOURCE'S ACTIVE range, not the buffer's.** A no-op
   for every buffer that already resampled; for a doubly-ghost buffer the along-seam
   direction is only ng deep in the BUFFER and clamping there would extrapolate from two
   cells when the source block holds the data.
3. **`FillPanelCorners*` fires only on a TRUE cube vertex** (BOTH flanking faces seams).
   The loose "either face" test also caught ordinary corners that DO have a real diagonal
   neighbour and overwrote properly exchanged data with a one-sided extrapolation.
   **That is what masked (1).**

    before   dlev 0 seam 0  x2x3 EDGE / 3D CORNER   0.2077 cells
    after    every category, every level relation   0.0000 cells

**No `cs_seam == 0` category survives anywhere.** Justification is CORRECT PLACEMENT, not
a measured accuracy win: 4 of 5 dumps differ, but refined resistive is unchanged, ideal
MHD+SMR moves only marginally (L1(B) 2.924630e-04 -> 2.924587e-04) and hydro is identical.
Bitwise at 1/2/3/6 ranks and 1/2 GPUs.

**A MEASUREMENT ERROR THAT NEARLY COST THE FIX:** the first before/after dump comparison
said "0 of 5 differ" only because the `git stash` meant to build the baseline had SILENTLY
FAILED from the wrong working directory -- both runs used the same binary. Redone properly:
4 of 5. **Check that your baseline is actually a different binary.** Second time this
session a bad instrument nearly produced a wrong conclusion.

## THE SECOND AND LARGER HALF: coarse_b0's cube-vertex corner was filled by NOTHING (180a9b3e)

A cube-vertex corner's exchange is skipped as non-reciprocal and `FillPanelCorners*`
extrapolates it instead -- **but that fill only ever ran on the FINE array.** The COARSE
array's corner is skipped by the same predicate and written by nothing, and `ProlongateFC`
interpolates exactly that block into the fine ghosts at a level boundary.

**It is READ, not merely stale**: poisoning cb's cube-vertex corner with 1e30 moves the
prolongated fine EDGE halo by 24%; nothing blows up only because the prolongation limiter
clips it. Fix: call `FillPanelCornersFC(cb, true)` / `FillPanelCornersCC(ca, true)` with
the coarse index set, when multilevel.

    cb's own corner   2.09e-01 (= |B|, garbage)  ->  1.86e-03
    r x EDGE SEAM     2.3278e-03 -> 4.6508e-04  (5.0x)   order 2.03 -> 2.78
    r x CORNER SEAM   1.8152e-02 -> 2.0911e-03  (8.7x)   order 2.01 -> 3.57

Both now sit at the SAME-LEVEL halo's own accuracy (3.90e-04, 1.98e-03). **Ideal MHD +
refinement improves independently**: Linf(B) 3.158684e-03 -> 2.332732e-03 (nx2=16) and
1.568119e-03 -> 9.820341e-04 (nx2=32), L1 still converging at 2.72. Hydro SMR and
rigid-rotation BYTE-FOR-BYTE identical; bitwise at 1/2/3/6 ranks and 1/2 GPUs.

**HOW IT WAS FOUND, after several wrong turns:** the residual survived EXACT inputs (not
the inputs); the pack-kernel mapping dump put 24 of 30 buffer categories at EXACTLY 0.0000
cells (not placement); and the worst cell then turned out to sit on a FINE block whose
ghost is PROLONGATED -- a path nothing had examined. Probing `coarse_b0` directly showed
its active and edge ghosts healthy (7.6e-05, and 3.9e-04 converging at 3.89) and its
CORNER flat at |B|. **The lesson: when a halo cell is wrong, first establish WHICH ROUTINE
WRITES IT. Twice this session the answer was "nothing does."**

The cell-centred twin was made symmetric and is INERT (hydro SMR unchanged), kept only
because it is the analogue of a demonstrated defect in the paired path.

## THE FIRST HALF: cs_seam classified from buffer EXTENTS (a01ace75)

`cs_seam` -- which decides which way the along-seam resample runs -- was inferred from the
buffer EXTENTS: `nj == ng && nk > ng` means an x2-normal seam, etc. **That test is only
valid at a SAME-LEVEL boundary.** At a coarse/fine boundary a FACE-CENTRED buffer carries
one extra layer in the direction its component is STAGGERED in, so the ghost-direction
extent is ng+1 = 3 rather than ng, the equality fails, and `cs_seam` silently falls to 0 --
skipping the resample for exactly ONE of the two angular components on every cross-level
seam edge buffer. Measured, x3f on the x3x1 edge buffers: same level nk = 2 = ng and the
resample runs, cross level nk = 3 and it does not. **The radial component is not staggered
in either tangential direction, kept extent ng, and kept its resample -- which is exactly
why x1f was ten times less wrong all along, the clue that was in front of me the whole
time.** Fix: classify off the SLOT INDEX (0-7 x1 faces, 8-15 x2 faces, 16-23 x1x2 edges,
24-31 x3 faces, 32-39 x3x1 edges, 40-47 x2x3 edges, 48-55 corners), which is what the
extents were a proxy for and is level-independent. `bvals_cc.cpp` got the same treatment
so the two cannot drift apart.

**THE MEASUREMENT THAT SETTLED IT -- reuse this.** Instrument the pack kernel to compare
the PHYSICAL POSITION of the RESAMPLED source sample against the destination cell that
buffer offset actually lands on. The destination index is recoverable inside the kernel:
take the receiver's own index set for slot `dest` (`rbuf[dn].isame/icoar/ifine[v]`, chosen
by the receiver's level relation) and invert its `i-il + ni*(j-jl + nj*(k-kl))`.

    cs_seam 2 or 3 (resample runs)   0.0000 cells   at EVERY level relation
    cs_seam 0, same level            0.2077 cells
    cs_seam 0, cross level           1.8242 cells

A fixed number of CELLS is O(h) -- that IS the first order. Nothing else in the chain was
ever off by more than round-off.

    r x EDGE SEAM  1.8223e-02 -> 2.3278e-03 (nx2=16),  8.4520e-03 -> 1.1447e-03 (nx2=32)
    cells > 1e-2 in that category: 150 -> 0;  domain L1(B) 1.3902e-05 -> 1.3054e-05

**IT IS STILL FIRST ORDER.** After the fix the edge halo goes 2.3278e-03 -> 1.1447e-03
and the corner 1.8152e-02 -> 9.0263e-03 over nx2 16 -> 32, both ratio ~2.0. The magnitude
is 7.8x smaller and the ORDER is unchanged, so a second mechanism of the same order
remains. Do not read a01ace75 as closing this.

Same-level runs are BITWISE unchanged (the slot rule reproduces the extent rule there);
refined resistive still bitwise at 1/2/3/6 MPI ranks.

## SCOPE: FACE-CENTRED, CROSS-LEVEL, CROSS-PANEL.  Hydro is NOT affected.

Measured, not argued: instrument BOTH `bvals_cc.cpp` and `bvals_fc.cpp` to compute what
the OLD extent-based rule would have given and compare it with the slot-based one, per
buffer per step.

    MHD refined:    24 x  FC dlev=+1  old=0 new=3   nj=10 nk=3  ng=2
                    24 x  FC dlev=+1  old=0 new=2   nj=3  nk=10 ng=2
                    24 x  FC dlev=-1  old=0 new=3   nj=8  nk=3  ng=2
                    24 x  FC dlev=-1  old=0 new=2   nj=3  nk=8  ng=2
    HYDRO refined:  nothing at all

**Zero CC disagreements**, in a hydro run AND in an MHD run's own cell-centred exchange.
Cell-centred buffers are not staggered, so their ghost extent stays ng at cross-level and
the old rule was already right there. Four conditions are each NECESSARY: face-centred
data, a level boundary, a panel seam, and the cubed sphere. **The production
deep-hot-Jupiter runs are SPHERICAL POLAR, so none of this touches them.**

Caveat kept honest: the hydro face-halo gate on a refined mesh is BETTER than its
unrefined control (6.0e-05..2.7e-04 vs 2.4e-04..6.3e-04, just finer cells), but that gate
reports per-panel FACE halos, not cross-level seam EDGES, so it supports rather than
proves the hydro case. And the RESIDUAL term is not proven as narrow: the doubly-ghost
buffers take cs_seam = 0 by design in bvals_cc too, so the cell-centred corner halo gets
the same plain copy; whether that is first order for CC has NOT been measured. The
standing argument that a dimensionally split PLM+HLLC sweep never reconstructs through the
diagonal block is an argument, not a measurement.

## TWO ATTEMPTS AT THE REMAINING TERM THAT ARE INERT -- do not repeat them

Both were built, fired, and measured to change NOTHING on any gate (halo Linf per category,
domain L1/Linf, heating ratio, all identical to every printed digit), so both were REVERTED
rather than committed. They are not obviously wrong -- they are unmeasurable, which in this
part of the code is the same reason not to ship them.

1. **Extend `cs_seam` to the doubly-ghost slots 40-55.** For a NON-cube-vertex corner only
   one of the two flanking faces is a seam, so the neighbour table does name an along-seam
   axis even though the slot alone does not: read the flanking face panels (slots 8/12 and
   24/28, as `IsCubeVertexCorner` does) and set cs_seam 2 or 3 from whichever is the seam.
   VERIFIED TO FIRE (slots 40,42,44,46,48+, both outcomes, dozens of times per step) and
   still changed no number.
2. **Take the resample stencil bounds from the SOURCE'S ACTIVE range instead of the buffer
   range** (`blo = ks_`, `bhi = ke_ + (staggered?1:0) - 2`). Confirmed a no-op for every
   buffer that already resampled, which is the point -- it exists so a doubly-ghost buffer,
   whose along-seam direction is only ng deep in the BUFFER, can still interpolate from the
   source block's real data instead of extrapolating from two cells.
3. **Tighten `FillPanelCorners*` to fire only on an ACTUAL cube vertex** (BOTH flanking
   faces seams, matching `IsCubeVertexCorner`) instead of when EITHER is. The loose test
   overwrites ordinary corners that have a real diagonal neighbour, discarding exchanged
   data -- which is why (1) looked masked. Also inert here.

The likely reason all three are inert: in these meshes the corner cells that any gate
actually measures are TRUE cube vertices, filled by the `FillPanelCorners*` extrapolation
and never by a buffer at all.

## THE CORNER EXTRAPOLATION: measured, and it is an ACCURACY-vs-AMPLIFICATION trade

`FillPanelCornersFC` extrapolates the cube-vertex corner QUADRATICALLY from the two
flanking strips, weights ((d+1)(d+2)/2, -d(d+2), d(d+1)/2). Those amplify strip error by
|3|+|-3|+|1| = 7 at d = 1 and 6+8+3 = 17 at d = 2. Swapping to LINEAR (weights d+1, -d;
amplification 3 and 5) was measured:

| r x CORNER SEAM Linf | quadratic | linear |
|---|---|---|
| refined (cross-level strips) | 1.8152e-02 | **6.1001e-03** |
| unrefined (same-level strips) | **1.9807e-03** | 5.4389e-03 |

3x better where the strips are noisy, 2.7x WORSE where they are accurate -- so the corner
error is AMPLIFIED STRIP ERROR at a level boundary and the extrapolation's own TRUNCATION
at same level. **Reverted: it would regress the certified same-level path.** The lesson is
that the corner is downstream, so the way to fix it is to fix the EDGE halo it reads, not
to retune the extrapolation. Do not re-propose the linear form on its own.

## THE CONCRETE ANOMALY TO START FROM

`r x CORNER SEAM` got WORSE when the edge halo it is built from got 7.8x BETTER:
1.3124e-02 -> 1.8152e-02 at nx2=16 across a01ace75. `FillPanelCornersFC` extrapolates that
corner quadratically from the two flanking strips, which at a RADIAL-GHOST layer are the
x1x2 / x3x1 edge halos. If the corner were simply inheriting strip error it would have
improved by the same 7.8x. It moved the other way, which means the old corner value was
benefiting from a CANCELLATION between the strip error and the extrapolation's own
truncation error. **So the remaining first-order term is the cube-vertex corner
EXTRAPOLATION itself, not the buffers** -- and it is now the only thing left above 1e-2.

## WHAT IS STILL OPEN, and is now the leading term

The **x2x3 EDGE and 3D CORNER buffers (slots 40-55)** are ghost in BOTH tangential
directions, have no single along-seam axis, and take NO resample at all -- by design.
Their mismatch is the 1.8242 cells above, so they are still O(h) at a level boundary.
`r x CORNER SEAM` therefore does not improve (1.3124e-02 -> 1.8152e-02 at nx2=16; the
corner extrapolation reads the edge halos, so changing those moves it), and the edge
category, though 7.8x smaller, still converges at first order. **Closing it needs a
TWO-dimensional resample for the doubly-ghost buffers.** Ideal MHD + refinement barely
moves either way (L1(B) 2.931187e-04 -> 2.931095e-04, Linf 3.158684e-03 -> 3.207754e-03)
because ideal stencils scarcely reach these cells.

## What it is

A halo region that is BOTH a panel seam and a level boundary -- the x1x2 / x3x1 EDGE
buffers, plus the cube-vertex corner that extrapolates from them -- is only FIRST order.
The three-way split is unambiguous (iprob=11, nx2=16, absolute Linf against the exact
static B, |B| ~ 0.2):

| | Linf |
|---|---|
| same-panel, level boundary | 2.7e-04 |
| same-level, panel seam | 3.9e-04 |
| **panel seam AND level boundary** | **1.8e-02** |

First order: 1.82e-02 -> 8.62e-03 over nx2 16 -> 32 (ratio 2.1). Present in FULL after
ONE step at tiny CFL, so it is halo FILLING, not evolution. Flat along the seam, uniform
over all six panels, both edge directions, and it hits the two ANGULAR components ~10x
harder than the radial one (1.36e-02 / 1.82e-02 against 1.78e-03).

**NOT a resistivity bug.** An `eta = 1e-10` run reproduces it identically. It is a
property of the FIELD halo, so the ideal-MHD refinement path (supported since cf0eb77c)
has carried it all along -- ideal stencils simply never read those cells; the gnomonic
resistive curl does, via BcovXi/BcovEta at (i-1, j-1). This is why resistivity + SMR was
still LIFTED (c4735f18): the domain norms converge and beat the unrefined control, and
refusing resistivity for a defect the supported ideal path shares is not coherent.

## RULED OUT -- do not re-run these

* **Buffer classification.** Printed it from inside the kernel: every cross-level
  cross-panel edge buffer gets `cs_seam` = 2 or 3 correctly (nj/nk = 2 vs 8/10, ng = 2).
* **The angle the seam transform is evaluated at.** Scanning a deliberate shift of the
  angle feeding BOTH the basis rotation and the resample puts the null at ZERO.
* **The resample magnitude.** Scanning a multiplier on the along-seam offset gives a flat
  minimum over 0.75-1.25 in the bulk count and 0.5-1.5 in Linf. It is right at 1.0.
* **Reciprocity and buffer sizing** (75a225d1, now a standing gate): 0 non-reciprocal and
  0 size mismatches over 1056 slots refined and 204 unrefined. So NOT a wrong slot and
  NOT a wrong stride.
* **Ghost width** (nghost=3 vs 2) -- ruled out earlier, see [[cubed-sphere-resistivity]].
* **Outer-layer-only.** The FIRST ghost layer is bad too (8.4e-03), so it is not an
  end-of-stencil artefact in the outermost layer.
* **One direction only.** Both the F->C (restriction) and C->F (prolongation) directions
  show it, so it is not specific to either operator.

## CONFIRMED engaged

Disabling the basis rotation on cross-level buffers makes it 17x WORSE (0.315); disabling
the along-seam resample makes it 1.6x worse (2.98e-02). So the machinery runs, at the
right angle, with the right offset -- and still leaves O(h) that same-level and
same-panel transfers do not have.

## A LEAD THAT DID NOT HOLD UP -- and the method lesson

A POSITION PROBE (evaluate the exact solution at neighbouring indices and see which
matches the halo value best) gave a 330x better match at offset (0,+2,-1) for the worst
cell, which read as "right data, wrong place". **Do not trust that on its own.** Over 914
bad cells the best offsets SCATTER rather than agreeing, and with a smooth field and a
7x7x7 search some offset matches by luck. The reciprocity/size audit then ruled out
wrong-slot and wrong-stride outright. **A brute-force best-match search over a smooth
field is not evidence; check how many candidate offsets there are and whether the winners
agree before believing one.**

## THE CO-LOCATION PATH IS EXONERATED, and the bug is now sharply cornered

Measured directly (temporary probe, numbers kept here because the probe is gone):

* **The restricted array is fine.** On the FINE blocks -- whose `coarse_b0` is what the
  `cs_coar` seam path actually reads -- the restricted primary b.x2f is 7.6e-05 from the
  analytic value, and the co-located SECONDARY that `srcval` builds is 3.2e-04 / 5.7e-05 /
  5.1e-05 over nx2 16/32/64, i.e. SECOND order. The clamped range ends are no worse than
  the interior (2.9e-04 vs 3.2e-04), so the one-sided extrapolation is not the problem.
* **The resolution argument does not apply.** Refinement only ADDS finer blocks; the
  coarse block's own cells are the same size as in the unrefined control, and the fine
  block restricts to exactly that resolution. So the inputs are the same quality as the
  same-level case, and the expected halo accuracy is the same 3.9e-04. **It is 1.8e-02,
  47x worse.**
* **The cross-level seam buffer IS the writer.** SENTINEL TEST, which is the technique to
  reuse: add 1000.0 to the value `srcval` returns when `nghbr.lev != mblev`, then have the
  gate skip any halo cell whose error exceeds 100. Run with `time/nlim=0` so no evolution
  can smear the tag. The first-ghost-layer error of the category drops from 8.4e-03 to
  **8.3e-05, a factor of 100** -- every bad cell in that layer is written by that path.

So: a buffer whose inputs are second-order accurate, whose reciprocity and strides are
verified consistent (including the comp-1 <-> comp-2 swap: sender ni/nj/nk = receiver
ni/nk/nj, checked by hand on both the working and the failing case), and whose rotation
angle and resample offset both null at their current values, still emits an O(h) error.
**The defect is inside the cross-level branch of the seam transform itself, not in what
feeds it and not in where it lands.**

## THE REPLICATION TEST: inputs and evaluation point BOTH exonerated IN SITU

Rather than replicate `srcval` on the host (which risks reimplementation bugs), feed the
REAL code exact inputs and see what survives. Inside `srcval`, for cross-panel buffers,
overwrite `bxi` and `bet` with the ANALYTIC iprob=11 values at the sample's own (xi,eta,r)
-- built from `cubed_sphere::PanelToCart` + `PanelTangents`, which are already available
in bvals_fc.cpp -- and leave the rotation, the resample, the index map and the placement
untouched. **VALIDATE THE INSTRUMENT ON THE WORKING CASE FIRST**, which it passes: the
same-level seam halo improves 5x (r x EDGE SEAM 3.90e-04 -> 7.46e-05, tangent ghost
2.26e-03 -> 1.09e-03), so exact inputs demonstrably reach the right place.

    cross-level r x EDGE SEAM:  1.8223e-02  ->  1.8179e-02     (0.24 %)

**With EXACT inputs the cross-level halo is unchanged.** The inputs are conclusively
exonerated, in situ, not by a side probe.

Then scan the point at which the analytic field AND the source tangent basis are
evaluated, one index direction at a time (radial / xi / eta, in cells), still with exact
inputs. **All three null cleanly at zero** -- xi: 2.82e-02 / 2.32e-02 / **1.82e-02** /
2.64e-02 / 3.44e-02 over -1,-0.5,0,+0.5,+1, and eta the same shape. So the source sample
sits at the physical point the code believes. The resample multiplier and an added base
index shift also null at their current values under the same clean conditions.

## THE SIGNATURE TO EXPLAIN: a symmetric V along the seam

With exact inputs, the residual on the coarse block's x3x1 edge ghosts, bucketed by the
along-seam index j (0..15 active):

    j    0     1     2     3     4     5     6     7     8 ...
    max 1.82  1.72  1.56  1.34  1.08  0.79  0.71  0.66  0.66  (x 1e-2, mirror-symmetric)

Maximum at BOTH panel corners, minimum at the seam midline, symmetric. **The error tracks
the magnitude of the along-seam OFFSET** -- the very quantity the resample exists to
cancel, which is zero at the midline and largest at the corners -- on top of a ~6.6e-03
floor that the midline still carries. Since that offset is a FIXED NUMBER OF CELLS
(resolution-independent), an offset that is not fully applied costs cells x h = O(h),
which is exactly the observed first order.

The natural reading is that the resample cannot reach the sample it needs at cross-level,
because a fine sender covers only HALF the seam and the 3-point stencil clamps at its
range end -- which would also explain why scaling the multiplier does nothing (already
clamped). **But that reading does not survive the j-distribution**: clamping would hurt
most where a sender's range ENDS, which at cross-level includes the MIDDLE of the coarse
block's range, and the middle is where the error is SMALLEST. So the mechanism is still
not identified. Do not adopt that explanation without testing it.

## Where to start next

The transform is applied to the right data at the right angle with the right offset, and
lands in the right place, yet is O(h) wrong only when the buffer crosses BOTH a seam and a
level. Everything that FEEDS the transform is now excluded, and so is the point it is evaluated
at. What is left is the correspondence between a SOURCE index and the DESTINATION cell it
serves -- and note that no UNIFORM shift improves it, so if that correspondence is wrong it
is wrong by an amount that VARIES from cell to cell. The next measurement is therefore not
another scan: instrument the pack kernel to emit, for one cross-level seam buffer, the
source index it reads and the destination index that buffer offset lands on (invert the
receiver`s `i-il + ni*(j-jl + nj*(k-kl))` on the host), and compare the two physical
positions cell by cell. That is the one link in the chain never directly observed.

Reproduce: `athena -i inputs/tests/cubed_sphere_resist_smr.athinput time/nlim=1
time/cfl_number=0.003` and read the `GHOST SCAN r x EDGE SEAM` line.
