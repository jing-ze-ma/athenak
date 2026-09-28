---
name: cs-resistive-seam-order
description: WHY the resistive current density is only 1st order on a panel seam and ~0.65 at a cube vertex, and why NO interpolation upgrade can fix it -- the halo cannot be more accurate than the O(h^2) solution it is built from. Records the 4-point experiment that was tried and REVERTED
metadata:
  type: project
---

Measured 2026-08-31 on `cubed_sphere_resist` (eta = 1e-3 so the hyperbolic CFL binds and
nx2 = 64 is affordable -- see [[cs-mhd-validation]]).

## THE CHAIN, established by measurement not argument

```
  resistive EMF (= J) L1 order    panel INTERIOR 1.91, 1.82     seam ring 1.39, 1.20
  tangent-ghost HALO error        first layer    2.00, 2.11     all layers 1.78, 1.76
```

The seam halo is **second order in VALUE**. The current density DIFFERENCES it over h, so
it can only be first order: `O(h^2)/h = O(h)`. In the panel INTERIOR the same counting
would give the same answer, and it does not, for one reason -- the interior truncation
error is a SMOOTH function of position, so the errors on opposite sides of the Stokes loop
CANCEL to O(h^3) and J comes out second order. At a seam the ghost values carry the
NEIGHBOUR panel's error field, which is a different function; the cancellation is lost.

## WHY NO INTERPOLATION UPGRADE CAN FIX IT

**The halo can never be more accurate than the data it is interpolated from.** The source
is the neighbour's discrete solution, itself O(h^2). Raising the interpolation order raises
nothing once the input error dominates. Verified, not assumed:

- **TRIED AND REVERTED:** the secondary-component interpolation in the seam transform
  (`x3f_at_xiface` / `x2f_at_etaface` in bvals_fc.cpp) is a 2-point average, O(h^2) in
  value. Replacing it with the 4-point centred midpoint form `(-1, 9, 9, -1)/16`, O(h^4),
  gave **-13 to -19% on the seam EMF error and NO change in ORDER** (1.46, 1.23 against
  1.39, 1.20) and no change at all in the tangent-ghost halo order (1.78 both). Reverted:
  a 15% constant that churns every recorded cubed-sphere number is not worth it, and it is
  not the fix. **Do not retry it expecting the order to move.**
- The along-seam resample is already a 3-point quadratic, O(h^3) in value.
- The MONOTONICITY LIMITER is not the cause: it is guarded (skipped at smooth extrema and
  whenever the resample extrapolates), and on smooth monotone data the quadratic already
  lies inside the node range, so it is a numerical no-op.

## THE PAPER ANALYSIS -- and then the MEASUREMENT THAT REFUTED IT

**READ THE REFUTATION AT THE BOTTOM BEFORE USING ANY OF THIS.** The argument below is
recorded because its corollaries about interpolation and partial circulation still stand
on their own (and one of them was confirmed by experiment), but its central premise --
that the discriminator is a JUMP in a smooth truncation-error field -- was tested directly
and is FALSE.

## THE PAPER ANALYSIS -- done, and it REFUTES the partial-circulation fix

Write the discrete solution on panel P as B + e_P with e_P = h^2 phi_P(x), phi_P SMOOTH
and set by the scheme contracted with P's OWN chart geometry. J = (1/A_e) * closed integral
of B.dl over the dual loop, A_e = O(h^2), each side carrying dl = O(h).

**INTERIOR.** Opposite sides of the loop are neighbouring points of the SAME panel:
`dl*[e_P(x+h) - e_P(x)] = h * h^2 * (h phi_P') = O(h^4)`, divided by A_e gives **O(h^2)**.
The cancellation IS the mechanism, and it needs phi smooth across the loop.

**SEAM.** One side is interior to A, the opposite side is a ghost carrying B's solution:
`dl*[e_A - e_B] = h * h^2 * (phi_A - phi_B) = h^3 * dphi`, divided by A_e gives
**O(h) * dphi**. First order, with the coefficient set by the JUMP in the TRUNCATION-ERROR
FIELD across the seam, dphi = phi_A - phi_B.

Three corollaries, and they close off everything cheap:

1. **Interpolation order is irrelevant.** The interpolation error enters only at O(h^4)
   after multiplying by dl; the O(h) term contains none of it. Predicted before it was
   re-checked, and it matches the 4-point experiment exactly: -15% constant, order 1.39 ->
   1.46 and 1.20 -> 1.23, i.e. unchanged.
2. **The PARTIAL-CIRCULATION exchange cannot work.** Whoever computes the loop, the side
   in B's territory carries phi_B and the side in A's carries phi_A. Rearranging who
   computes it, or exchanging partial sums instead of finished values, reproduces the same
   phi_A - phi_B. **This kills the design I had proposed -- do not build it.**
3. **Averaging the two panels' J cannot work either.** When B computes the same edge, its
   OWN side carries phi_B and its GHOST (from A) carries phi_A -- the SAME pairing, not the
   opposite one. There is no antisymmetry to exploit. (Worth knowing because the seam EMF
   is currently made single-valued by a DIRECT COPY, not an average -- see the note at the
   top of flux_correct_fc.cpp -- so averaging looks like a free thing to try. It is not.)

**WHAT dphi ACTUALLY IS.** The equiangular gnomonic grid IS mirror-symmetric about each
seam plane, so an exactly equivariant scheme would give phi_B = R phi_A and, on the seam
(fixed by R), dphi = 0 for components EVEN under the reflection and dphi = 2*phi_A for the
ODD one. Consistent with the per-component measurement: domain L1 order 1.94/1.86 for the
RADIAL edge EMF (the radial field crosses a seam as a SCALAR) against 1.67/1.53 and
1.80/1.68 for the two tangential ones.

**THE ONLY ROUTES LEFT**, none cheap:
  (a) raise the INTERIOR scheme to 3rd order so e = O(h^3), the jump is O(h^3) and J is
      O(h^2) -- what the high-order cubed-sphere literature does, a whole-scheme change;
  (b) make the discretization EQUIVARIANT under the seam reflection for the odd component
      so dphi = 0 by construction -- targeted but research-grade, no guarantee;
  (c) accept it: a 2-cell ring, and the evolved-field DOMAIN L1 still converges at 2.07.

## HOW MUCH IT MATTERS

The evolved-field DOMAIN L1 still converges at ~2.07, because a 2-cell seam ring is a
vanishing fraction of the volume. The damage is confined to Linf and to the seam ring
itself. The cube vertex is worse (L1 order 0.68, 0.64) and [[cs-wire-fill-wip]]'s corner
fill lowers its error 9-17% **without changing the order** -- same reason.


## THE REFUTATION -- measured, with the control that mattered

New gate in cs_test (`ERROR JUMP in B.rhat across a shared face`): take the error in
`b.x1f` -- B.rhat, a CHART-INDEPENDENT SCALAR, so the two sides are directly comparable --
in the first GHOST cell (which holds the neighbour's solution) and the first ACTIVE cell,
and difference them. The premise predicts O(h^3) where one smooth error field spans the
face and O(h^2) where a jump exists. **The control is the point:** the same difference at
an ordinary block boundary WITHIN one panel, where no jump can exist. Run with 4 MeshBlocks
per panel so both exist (nx2 = 16, 32, 64):

```
                     jump L1 / max                                order (L1/max)
  PANEL SEAM     2.026e-06/1.443e-05 4.470e-07/4.935e-06 1.764e-07/1.555e-06  2.18/1.55, 1.34/1.67
  same-panel CTRL 5.451e-07/7.206e-06 1.199e-07/3.114e-06 2.773e-08/1.420e-06  2.19/1.21, 2.11/1.13
```

**The same-panel control shows the SAME O(h^2) jump as the seam, and its interior current
is nevertheless SECOND order.** So an O(h^2) difference of the field error across a face
is NOT what makes the seam current first order, and the premise is dead. The reason the
naive count fails is that for a face-centred CT field the truncation error of the CURL is
not the difference of the truncation errors of B -- the discrete field is a 2-form and the
curl is its exact discrete exterior derivative.

**WHAT SURVIVES, empirically:** the seam jump is 4-6x LARGER than the control and its order
degrades (2.18 -> 1.34) while the control holds ~2.1. Across both, `J order ~ jump order -
0.8` (control 2.11 -> J 1.9; seam 1.34 -> J 1.2). So there IS a real seam-specific
low-order component in the halo -- but it is NOT identified, and route (b) ("make the seam
stencil reflection-equivariant") has **no demonstrated target to cancel**.

**THREE HYPOTHESES HAVE NOW FAILED** -- higher-order interpolation (measured: order
unchanged), the partial-circulation exchange (refuted on paper), and the jump mechanism
(refuted by this gate). That is the signal to stop guessing and to instrument the halo
error's h-dependence component by component before proposing a fourth.


## THE DECOMPOSITION -- done, and it LOCALISES the defect

Four measurements, each with its own control (eta = 1e-3, nx2 = 16/32/64):

```
1. HALO MACHINERY ALONE (time/nlim=0, so the interior IS the analytic solution and the
   halo error is the interpolation machinery with nothing propagated into it):
     seam jump L1  1.552e-06  2.187e-07  2.920e-08     order 2.83, 2.91   <- THIRD ORDER
     same-panel control jump  EXACTLY 0.0                                 <- instrument OK

2. FULL evolved seam jump      2.026e-06  4.470e-07  1.764e-07     order 2.18, 1.34

3. EVOLVED FIELD ERROR BY DISTANCE from the panel edge, in cells:
     d0 order 1.30, 1.24     d1 1.37, 0.74     d2 1.89, 1.44
     d3 2.04, 1.95           d4 1.96, 2.18     >=5 2.27, 2.17
   -> the damage is CONFINED TO THE FIRST TWO CELLS and does not propagate.

4. RESISTIVE EMF AFTER ONE STEP (time/nlim=1: the interior is still essentially exact, so
   this is J built from an exact interior plus the third-order halo):
     seam ring 2.223e-03 8.872e-04 3.871e-04   order 1.33, 1.20
     panel interior 8.406e-04 2.149e-04 5.438e-05  order 1.97, 1.98
```

**CONCLUSION: the halo interpolation is NOT the cause -- it is third order. The seam
current is ALREADY first order after a single step, from an exact interior.** So it is not
accumulation, not a feedback, and not the interpolation. The defect is in **how J is
FORMED at a seam-adjacent edge**: pass 1 of the gnomonic curl, when a side of its Stokes
loop lies in a ghost cell.

**THE REMAINING CANDIDATE, not yet tested:** the seam transform moves FACE-AVERAGED flux
densities between panels as if they were POINT VALUES. That discrepancy is O(h^2) in the
value -- harmless for the value, and it CANCELS in the interior curl because it is a
smooth field there -- but the two panels' faces have different orientations and shapes, so
across a seam it does not cancel, and differencing it over h costs exactly one order. The
fix would be the standard high-order-FV average-to-point correction (a transverse
Laplacian term) applied before the transform and inverted after. **This is a HYPOTHESIS.
Three have already failed; test it before building it** -- e.g. by adding the correction
only to the diagnostic and re-measuring the one-step seam EMF order.

**INSTRUMENT WARNING:** `time/nlim=0` leaves `efld_resist` NEVER FILLED, so the EMF gate
then reports the error of an empty array (0.49, order 0.01) and looks like a catastrophe.
Use `nlim=1` for anything that reads the resistive EMF. See [[validate-the-instrument]].


## PER-COMPONENT SPLIT -- the answer, and the face-average hypothesis SURVIVES

Extended the jump gate to all three field components, mid-seam only, at `time/nlim=0` so
it is the halo MACHINERY alone:

```
  component      nx2=16     nx2=32     nx2=64      order
  x1f (B.rhat) 1.552e-06  2.187e-07  2.920e-08   2.83, 2.91   <- THIRD order
  x2f          2.167e-05  6.010e-06  1.590e-06   1.85, 1.92   <- second
  x3f          2.689e-05  8.067e-06  2.191e-06   1.74, 1.88   <- second
  same-panel control, all three components: EXACTLY 0.0
```

**The TANGENTIAL components are second order and the RADIAL one is third**, and the
tangential errors are 14x (nx2=16) to 75x (nx2=64) larger. A second-order halo differenced
over h is a FIRST-ORDER current, which is exactly the measured seam EMF order of 1.2-1.33.
Chain complete.

**The partner interpolation is definitively NOT the cause.** Re-tested the 4-point form
with this clean per-component gate: x2f 1.84/1.92 (was 1.85/1.92), x3f 1.87/1.91 (was
1.74/1.88). Order unchanged. Reverted again. The earlier refutation was right.

**WHY THE TWO COMPONENTS DIFFER, which is the face-average hypothesis made precise.**
`b.x1f` is B.rhat on a face of r = const -- the SAME physical surface in both charts, and a
chart-independent SCALAR, so it crosses a seam as a plain number and reaches third order.
`b.x2f` and `b.x3f` are flux densities over the TANGENTIAL faces, and the destination's
tangential face is a DIFFERENT physical surface from any source face -- different
orientation, different area. The transform rotates vector components at a POINT, but the
stored quantity is a face AVERAGE; treating one as the other is an O(h^2) error, and it
exists only for the components whose face is not shared.

Supported by four things: the order split (2 vs 3), the magnitude split (14-75x), the
elimination of the partner interpolation, and the shared-vs-unshared face argument. It is
CONSISTENT, not yet DEMONSTRATED -- that needs the average-to-point correction (subtract
the transverse Laplacian term before the transform, add it back after) and a measurement
showing the tangential order rise to ~3. **NOT BUILT.**


## SWAP vs NON-SWAP -- the SHARED-FACE RULE, and why the correction is not yet buildable

Split the per-component jump by seam type (nlim=0, halo machinery alone). Panels 3 and 5
are polar; an equatorial-to-polar seam is a SWAP seam:

```
  x1f (B.rhat)  NON-swap 2.80, 2.89     SWAP 2.84, 2.91
  x2f           NON-swap 3.13, 3.13     SWAP 1.83, 1.91
  x3f           NON-swap 1.77, 1.89     SWAP 1.68, 1.86
```

x2f on a NON-swap seam is **third order and 20x smaller** (1.567e-06 against 3.172e-05).

**THE RULE THIS ESTABLISHES: a component's halo is THIRD order exactly when its face IS
the shared seam surface, and SECOND order otherwise.** An equatorial-to-equatorial seam is
a xi-face, so there x2f's own face IS the seam -> third order; x3f's eta-face extends in
the perpendicular direction and is NOT shared -> second. b.x1f's radial face is r = const,
shared in every chart -> third everywhere. On a SWAP seam the shared face is x3f's from one
side and x2f's from the other, and this gate lumps both sides, so each component reads as a
mixture and both come out second. **So per seam, TWO of the three components are third
order and exactly ONE -- the along-seam tangential flux density -- is second.** That single
component is what makes the seam current first order.

**WHY THE CORRECTION WAS NOT BUILT.** The average-to-point argument says to subtract
(1/24)*(transverse Laplacian) before the transform and add the destination's back after.
But the two in-face directions are (radial, transverse-tangential) on both sides; the
RADIAL part is shared and cancels, and the code's own comment records that the two charts
**share the seam-NORMAL coordinate exactly**, so the remaining parts very nearly cancel
too. If they cancel, the correction is a no-op and building it is a fifth speculative
attempt. **The cancellation structure has to be worked out before writing the code** --
that is a paper question, not a coding one, and it is exactly where this stopped.


## THE PROFILE ALONG THE SEAM -- mechanism CONFIRMED, correction NOT closed

Binned the non-swap-seam jump by u = |along-seam angle|/(pi/4); u0 is the seam MIDPOINT,
u3 the CUBE VERTEX (nlim=0, nx2=64, machinery alone):

```
  x1f (B.rhat, shared face)  6.34e-09  1.65e-08  3.36e-08  5.55e-08
  x2f (SHARED face here)     9.59e-09  3.16e-08  2.74e-08  8.29e-09
  x3f (NON-shared face)      9.74e-07  2.56e-06  4.93e-06  7.30e-06   <- 7.5x, MONOTONE
```

The non-shared-face component grows **monotonically 7.5x from the seam midpoint to the cube
vertex**, and is 100-700x larger than the two shared-face components. That is exactly the
predicted signature, and the SAME geometric quantity explains why the cube vertex is the
worst region of all (order 0.65 there against ~1.2 on the seam).

**THE PAPER RESULT.** For the non-swap +x/+y seam at P ~ (1,1,t), the two charts' xi
directions are  e_A ~ (-1, 1+t^2, -t)  and  e_B ~ (1+t^2, -1, -t).  They are ANTIPARALLEL
at t = 0 and non-parallel for every t != 0 (equal component ratios would need
(1+t^2)^2 = 1).  Expressed in the source basis {e_xi, e_eta}, the destination's transverse
direction has cosines (-1, 0) at the midpoint and (-1, -1) at the vertex -- from pure xi to
a 50/50 mix.  So the misalignment vanishes at the midpoint and is maximal at the vertex,
matching the measured profile.

**WHERE IT IS STUCK, and this is why nothing was built.** Two sub-results point opposite
ways and I could not reconcile them:
  - The SPACING part of the correction provably CANCELS. Both charts are equiangular with
    the same angular width, and |e_xi| = r(1+X^2)sqrt(1+Y^2)/delta^2 evaluates to the SAME
    number in both charts at a seam point (source +x: X=1, Y=t; destination +y: X=1, Y=t).
    So (ds_dst/ds_src) = 1 exactly and the simple ratio-form correction is a NO-OP.
  - The DIRECTION part looks nonzero by the geometry above -- but in INDEX space the
    destination's transverse direction and the source's jj direction traverse the same
    sequence of physical depths (the two charts share the seam-normal coordinate exactly),
    which argues the second differences coincide and the correction cancels after all.

**RESOLVED by writing the ghost FACE out explicitly.** Both were partly right. The
seam-NORMAL coordinate maps EXACTLY (xi_A + xi_B = pi/2), so the ghost layers coincide with
the source's active layers and nothing is wrong in that direction -- that is why the
spacing part cancels. The faces are bounded in the OTHER direction, and there the
receiver's face is a surface of constant eta_A = atan(z/x) while the source's is constant
eta_B = atan(z/y): off the seam, DIFFERENT SURFACES, SHEARED relative to each other. In the
source chart the receiver's face has eta_B drifting with xi_B at rate

    sigma = d eta_B / d xi_B = sin(2 eta_B) / sin(2 xi_B)

which is ZERO at the seam midpoint and rises to 1 at the cube vertex -- the measured
profile. Averaging over a line tilted by sigma rather than a flat one gives, the linear
term killed by centring,

    <f>_recv - <f>_src = sigma * (dxi^2 / 12) * d2f/dxi deta + O(sigma^2 h^2)

second order, LINEAR in sigma, a MIXED derivative.


## THE SHEAR CORRECTION: BUILT, DERIVED, VALIDATED -- and it does NOT fix the order

Implemented in `srcval` (bvals_fc.cpp), behind `CS_SHEAR` and **DEFAULT OFF**; a default
build is BITWISE identical to HEAD (checked against a HEAD build of the same file, not
against a recorded baseline -- the first comparison used stale dumps and read as a
regression, [[validate-the-instrument]] yet again).

- **Sign settled empirically, and the signature is exactly right:** +1 DOUBLES the halo
  error, -1 HALVES it. That is what a correct term with the right magnitude looks like.
- **Halo (non-swap seam, the non-shared-face component):** 4.930e-05 -> 2.820e-05,
  1.444e-05 -> 7.228e-06, 3.885e-06 -> 1.891e-06. Order 1.77/1.89 -> **1.96/1.93**.
  A 2x error cut, but still SECOND order, not third.
- **Resistive EMF seam ring:** 1.33/1.20 -> 1.38/1.23, error down 9-15%. Panel interior
  unchanged at 1.97/1.98, so it costs nothing elsewhere.
- **STENCIL TRAP that cost an hour:** a seam buffer packs cells AT THE EDGE of the active
  range, so a CENTRED difference across the seam normal reaches outside the data and the
  guard never fires -- silently, with no effect at any amplitude. Found by poisoning the
  block with 1e30 and bisecting the conditions. Use a ONE-SIDED difference there.
- **Only fires on NON-SWAP seams.** On a swap seam `vv != v` and the component condition
  picks the wrong face; fixing that should roughly double the benefit. NOT DONE.
- The 4-point partner interpolation ON TOP makes it 2.8x WORSE. Do not combine them.

**VERDICT: the correction is real, derived and harmless, but the seam current stays FIRST
order.** A 2x cut in the one halo component moved the EMF by 10%, so the seam current is
dominated by some OTHER O(h^2) contribution that is not this halo component. That is the
next thing to find, and the per-component/position gates are the instrument for it.


## THE SEAM CURRENT'S SPATIAL STRUCTURE (2026-09-01, sixth instrument)

The EMF ring gate now bins by distance to the nearest CUBE VERTEX, in CELLS and then in
PHYSICAL fraction of the half-seam (nlim=1, exact interior, eta=1e-3, nx2=16/32/64).
**Bins in cells confound order with profile** -- at fixed cell distance the bin tracks a
shrinking physical neighbourhood of the vertex, and the orders came out 0.3-0.95, which is
an artifact of that. At FIXED PHYSICAL distance:

```
                       shear OFF          shear ON
  <1/8 (vertex)       1.85, 1.62         1.92, 1.71     (10x the mid-seam AMPLITUDE)
  1/8-1/4             1.69, 1.13         1.80, 1.14
  1/4-1/2             1.27, 1.15         1.29, 1.16
  >1/2 (mid-seam)     1.59, 1.39         1.64, 1.42
```

**Two conclusions.** (1) The vertex-domination hypothesis is WRONG TOO -- the seam current
is ~1.1-1.4 order at every physical distance from the vertex, i.e. the first-order error is
DISTRIBUTED along the whole seam, not anchored at the corner; the vertex has 10x the
amplitude but actually the best order. (2) The shear correction improves amplitudes ~10%
everywhere and order marginally; it is not the dominant term anywhere on the seam.

**STANDING PUZZLE, sharper than before:** at nlim=1 the interior is exact and the halo is
3rd order in two components and 2nd (halved) in one -- yet J is ~1.2-1.4 order EVERYWHERE
on the seam. The value-accuracy of the halo CANNOT be the explanation; something about how
the seam-adjacent Stokes loops consume it (or the geometry/metric factors they use in the
ghost region, or the partner-average inside BcovXi/BcovEta evaluated in ghost cells, or the
seam EMF direct-copy) loses an order. The next instrument should compare a seam-adjacent
loop's DISCRETE circulation, term by term, against the same loop evaluated with the exact
analytic field at exactly the stencil's sample points -- that separates "bad inputs" from
"inconsistent loop" per term, in one run.


## CLOSED: the TERM-BY-TERM LOOP INSTRUMENT names the mechanism (2026-09-01)

New gate in CSTestResistCheck: recompute every x1-edge Stokes loop TWICE, verbatim from
resistivity_gnomonic.cpp -- once with the discrete field (Jd), once with the EXACT analytic
field at the stencil's own sample points on the extended chart (Je), same geometry factors
-- against the truth Jx = 2*b0c*(zhat.rhat).  nlim=1, eta=1e-3, nx2=16/32/64:

```
  RING      total    2.224e-03  8.873e-04  3.871e-04   order 1.33, 1.20
            OPERATOR 9.780e-04  2.480e-04  6.269e-05   order 1.98, 1.98   <- INNOCENT
            input    1.580e-03  7.280e-04  3.486e-04   order 1.12, 1.06   <- GUILTY
  interior  operator 1.97/1.98; input 30-600x smaller
```

**The OPERATOR -- loop, ghost-region geometry, BcovXi/BcovEta -- is second order
EVERYWHERE, ring included.  The first-order seam current is entirely INPUT error pushed
through the loop.**  The reconciliation with the halo measurements: an O(h^2) VALUE error
that is rough cell-to-cell becomes O(h) in J after the loop divides by area (~ eps/h) --
so the one non-shared-face tangential halo component (O(h^2), the shear-affected one) is
exactly enough to explain a first-order current.  And the 10% response to the shear
correction is explained by SEAM COUNT: the correction fires only on NON-swap seams, which
are 4 of 12 -- the 8 equatorial-to-polar swap seams dominate the ring L1 and were untouched.

**THE PATH TO SECOND ORDER, now with no mystery left:**
  (a) extend the shear correction to SWAP seams (the component condition must use the
      swap: vv != v there) -- mechanical, should cut the input term's constant ~2x;
  (b) even where it fires it only HALVES the O(h^2) halo error (residual order ~1.95) --
      the neglected O(sigma^2 h^2) tilt term and the first-order one-sided mixed
      difference are the candidates; the clean endgame is a QUADRATURE over the receiver's
      tilted face (2-point Gauss along the tilted segment, each sample via the existing
      along-seam resample at a shifted normal position), which makes the ghost
      face-average O(h^3) and J genuinely second order;
  (c) the vertex still needs its own story (corner fill is an extrapolation).


## STEP 1 DONE (all 12 seams), and the VERDICT ON A HIGHER-ORDER INTERIOR

The shear correction now fires on every seam.  Two of my three attempts at the swap
bookkeeping were wrong, and the third was a bug in my own stencil:

- **THE REAL BUG:** the mixed difference was keyed on WHICH COMPONENT is corrected while
  the inward neighbour index was keyed on WHICH AXIS is normal.  The moment those disagree
  -- exactly the swap case -- it indexes a j-neighbour as a k-neighbour.  Key the stencil
  on the AXIS ROLES and choose the component separately.
- **THE COMPONENT RULE DOES NOT CHANGE ACROSS A SWAP SEAM, and neither does the sign.**
  Flipping it (the obvious guess, since the transform maps components through the swap)
  measured WORSE with either sign.  The swap seams only LOOKED untouched because of the
  stencil bug.

```
  halo x3f jump      all seams        non-swap        swap
    OFF            1.74, 1.88       1.77, 1.89     1.68, 1.86
    ON             1.82, 1.89       1.96, 1.93     1.68, 1.85
    amplitude      -35%             -43%           -22%
  ring EMF total   1.33/1.20 -> 1.43/1.27   (error -23%);  operator UNCHANGED 1.98/1.98
```

Fixing 12 of 12 seams instead of 4 of 12 gave ~2.3x the benefit (10% -> 23%), exactly as
the seam-count explanation predicted.  **Still first order.**

## WOULD A 4th-ORDER INTERIOR (ppm4) FIX THE SEAM?  NO -- and it is already MEASURED

**The experiment has effectively been run.** `time/nlim=1` starts from the ANALYTIC initial
condition, so the interior data is exact to machine precision -- better than any 4th-order
scheme can ever be, an infinite-order interior -- and the seam ring still came out
**1.33/1.20**.  An interior upgrade cannot beat an exact interior, so it cannot fix this.

Two structural reasons, both measured:
1. The OPERATOR term (loop, ghost geometry, BcovXi/BcovEta) is 1.98/1.98 and is the
   RESISTIVE CURL, not the reconstruction; ppm4 does not touch it.
2. The INPUT term is set by the SEAM INTERPOLATION -- partner average O(h^2), shear
   O(h^2), resample O(h^3) -- all independent of the interior scheme's order.
A higher-order interior only pays if the halo machinery is raised to match, which is the
package the high-order cubed-sphere literature buys (spectral element / high-order FV).
The cube VERTEX is likewise an extrapolation and equally indifferent to the interior order.

**RECOMMENDATION: STOP.** The domain L1 is 2.07, mass and energy are exact, the mechanism
is fully explained and instrumented, and what remains is a scheme-wide high-order rewrite,
not a patch.
