---
name: cs-seam-destag-refuted
description: The 2-point de-staggering average is NOT the binding O(h^2) term at the cs seam; the SWAP-seam ghost is, and neither CS_DESTAG4 nor CS_SHEAR raises the order
metadata:
  type: project
---


**SUPERSEDED 2026-09-02 -- the question in this file is ANSWERED.**  The srcval
hypothesis at the bottom was RIGHT: an oracle A/B put every 1.8-order category at ~3.0
and a clamped-window cubic realises it.  See [[cs-seam-colocation-fixed]] for the fix,
the numbers, and the reason a CENTRED wide stencil measures exactly the baseline.  Read
the rest of this file as the record of what was refuted along the way.

2026-09-01. Ran the halo-only seam-jump refinement study (iprob=11, `nlim=0`,
nx2=32/64/128, `bench/cs_seamjump`, jobs 11309330/11309384/11309431 on apudev, ~16 s
each).  No timestepping; the existing gate in `CSTestResistCheck` already prints the
ghost-minus-active error jump, so step 1 needed NO code change.

**CONFIRMED, and it is a clean split.**  Components that arrive as a directly stored
shared-face value have an O(h^3) jump: x1f 2.80/2.90, x2f NONswap 3.11/3.07.  Every
component that goes through the ghost de-staggering is O(h^2): x2f SWAP 1.79/1.90,
x3f NONswap 1.79/1.89, x3f SWAP 1.75/1.88.  An O(h^2) jump differenced across the
ghost/active interface gives O(h) in J -- that is the measured seam order of 1.00.
LOOP SPLIT reproduced the recorded 1.98 operator / 1.06 input exactly.

**REFUTED: raising the de-staggering to 4th order does not fix it.**  Added
`CS_DESTAG4` (default OFF, host-read + capture-by-value like `CS_SHEAR`) switching
`x3f_at_xiface`/`x2f_at_etaface` in `src/bvals/bvals_fc.cpp` from `0.5*(f(-1)+f(0))`
to `(-f(-2)+9f(-1)+9f(0)-f(1))/16`.  x3f NONswap amplitude halves (1.29e-5 -> 6.27e-6)
and the order creeps to 1.91, but **LOOP input moves only 1.06 -> 1.09**.  The swap-seam
components get slightly WORSE.  So the 2-point average is a real O(h^2) term but NOT the
binding one.

**CS_SHEAR does not fix it either.**  All four (destag4, shear) arms leave LOOP input at
1.03-1.09.  `shear=-1` is the best AMPLITUDE (5.73e-4 vs 7.34e-4 baseline, 22% better);
`shear=+1` is 36% WORSE than off, so the sign matters and +1 is wrong.
`destag4=1, shear=-1` drives x3f NONswap to order 1.98/1.99 -- the only arm that reaches
2 on any category -- yet its amplitude is worse than baseline and LOOP input is unmoved.

**The binding term is the SWAP seam.**  x3f SWAP sits at 1.75-1.88 in EVERY arm, immune
to both switches -- consistent with the shear correction firing only on non-swap seams.
Any fix that does not touch the swap-seam ghost construction cannot raise the seam order.

**The cube vertex is still unexplained.**  The jump's vertex bin converges at the SAME
order as the midpoint (1.86-1.94), only with a ~9x larger constant, so the halo jump does
not explain the measured resistive vertex order of 0.65.  Note the gate deliberately
SKIPS the `ng` cells at each seam end, i.e. it never measures the true cube-vertex cells.
The claim that the corner is "strictly downstream of the seam" is UNTESTED.

**Why:** this kills the cheap 30-line input-side fix and redirects the work to either the
tilted-face quadrature or, cheaper, making the resistive operator not read the ghost at
all (one-sided quadratic extrapolation from active faces + the existing cross-panel EMF
averaging) -- which is immune to whatever remains in the ghost.

**How to apply:** `CS_DESTAG4` is left in the working tree, default OFF and therefore
inert; delete it or keep it as an A/B handle.  Before quoting any shear number, check
`CS_SHEAR` -- it defaults to 0.0 at `bvals_fc.cpp:65-67`, so no production run has ever
had the correction on.  See [[cs-resistive-seam-order]], [[cs-mhd-validation]],
[[validate-the-instrument]].

## Plan C (do not difference the ghost) -- PROTOTYPED, and it needs the vertex

Implemented as `CS_PLANC` (default OFF, bitwise safe): a working copy `bseam` whose
panel-seam ghost layer is overwritten by a one-sided quadratic extrapolation from this
panel's own active faces (`resistivity_gnomonic.cpp`, plus `bseam`/`cs_seamf` members).
Jobs 11309470/473/483.

**FIRST RESULT WAS VACUOUS -- the instrument was dead twice over.**  `nlim=0` NEVER
COMPUTES THE EMF, so both arms were bit-identical; the tell was `x1e max=1.0000e+00`
EXACTLY, which is `|0-analytic|/|analytic|`.  And `LOOP SPLIT` is a HOST-SIDE
reimplementation of the curl inside `cs_test.cpp` reading `b2h`/`b3h` -- it never calls
`AddEMFGnomonicResist` and CANNOT see any change to the operator.  Use the EVOLVED FIELD
gate with real timestepping (`tlim=0.01`), not `efld_resist` (pre-exchange) and not
LOOP SPLIT.

**Measured (evolved field, tlim=0.01, nx=32/64/128):** baseline L1 order 2.03/2.02,
plan C 0.98/0.61.  Linf 1.49/1.66 -> 0.49/0.20.  Plan C as prototyped is WORSE.

**But the premise SURVIVES.**  div B is IDENTICAL between arms to every digit (seam
2.036e-05, interior 1.648e-05 at nx=32) and Ohmic heating is unchanged (ratio 1.012939 vs
1.012979), so the two panels still agree on the shared seam EMF -- the existing
cross-panel averaging does reconcile the one-sided values, exactly as plan C assumed.

**The damage is the CUBE VERTEX.**  Panel interior is barely touched (2.12e-4 -> 1.93e-4);
the 2-cell panel-edge ring goes 5.60e-4 -> 1.00e-3 at nx=32 and 2.15e-4 -> 6.96e-4 at
nx=64, and the worst-error location MOVES to (34,34) / (66,66), i.e. the cube-vertex
corner.  In the prototype both de-ghost kernels build corner cells from raw unfilled
ghost, and the corner never gets the 2-D quadrant extrapolation.  So the vertex is not a
follow-on to plan C, it is a PRECONDITION for it.

**How to apply:** do not judge plan C on these numbers -- fix the corner first (2-D
quadrant extrapolation for cells where two seam flags are set), then re-measure.

## Plan C + the 2-D corner fill: STILL FAILS, and the failure is a GROWING MODE

Job 11309493.  Added kernel C: the four corner cells get the TENSOR PRODUCT of the two
quadratic extrapolations from the 3x3 block of ACTIVE faces, and kernels A/B were
restricted to the active transverse range so nothing is written twice and no
extrapolation reads another extrapolation.  Arm 0 still reproduces baseline exactly.

It did NOT fix it, and at nx=128 it is far worse:

  global L1   baseline 2.26e-06 5.52e-07 1.36e-07  order 2.03 2.02
              plan C   4.64e-06 2.56e-06 2.60e-06  order 0.86 -0.02
  global Linf baseline 1.93e-04 6.87e-05 2.17e-05  order 1.49 1.66
              plan C   2.08e-04 1.94e-04 2.37e-03  order 0.10 -3.61

**The error GROWS with resolution** (Linf 1.94e-04 -> 2.37e-03 from nx=64 to 128 at fixed
tlim), which is an amplifying mode, not truncation -- finer grids take more cycles to
reach the same physical time.  The worst cell is the CUBE VERTEX in every plan C run and
sits exactly at j=k=je+1: (34,34) at nx=32, (66,66) at 64, (130,130) at 128.  By nx=128
it has polluted the panel interior too (order -0.50).

div B is STILL bitwise identical to baseline (2.0360e-05 / 2.6180e-06 / 3.3210e-07) and
Ohmic heating is unchanged, so this is NOT a CT or telescoping failure.  It is positive
feedback: the corner value is extrapolated from active data, enters the curl, updates the
field near the corner, and feeds the next extrapolation.  The cross-panel ghost -- rough
as it is -- is what STABILISES the corner, because it carries the neighbouring panel's
information.

**CONCLUSION: plan C in the "extrapolate from your own panel" form is not viable.**  Any
future attempt must keep real neighbour information at the corner (blend, damp, or use
the true 2-D cross-panel sample of [[cs-cube-vertex-real-exchange]]) rather than replace
it.  Do not re-propose one-sided extrapolation without a stability argument.

## HYPOTHESIS for the second O(h^2) term (proposed 2026-09-01, NOT YET TESTED)

`srcval` co-locates the partner component by TWO SUCCESSIVE 2-POINT AVERAGES ON
ORTHOGONAL AXES, and fix 1 only touched one of them.  For `vv==1`
(`bvals_fc.cpp:366-369`): `bet = 0.5*(x3f_at_xiface(kk,..) + x3f_at_xiface(kk+1,..))`.
The INNER `x3f_at_xiface` (`:329-337`) interpolates across xi -- that is the one fix 1
raised to 4-point.  The OUTER `0.5*(kk, kk+1)` interpolates across ETA and is STILL a
2-point average.  Same for `vv==2` at `:372`.  Two comparable O(h^2) terms on orthogonal
axes: killing one HALVES the amplitude and CANNOT move the order -- exactly the measured
signature.

Why the SWAP seams are immune to both switches: the partner's error reaches destination
component `v` weighted by the OFF-DIAGONAL entry of `TransformFieldToDstNormals`
(`coordinates/cubed_sphere.hpp:245-253`).  On a NON-swap seam the destination component
whose face IS the shared surface has a vanishing off-diagonal weight, so it inherits the
clean primary -- hence x2f NONswap at 3.11.  Across a SWAP seam the axes exchange, NEITHER
destination face coincides with a source face, both components get O(1) off-diagonal
weight, and both sit at ~1.8.  The shear switch is structurally single-component (`vwant`
is one value, `:433-435`), so it can correct at most one of the two.

**Proposed fix:** replace `:366-372` with a single 2-D quadratic (or Balsara
divergence-preserving) reconstruction of the tangential pair inside the source cell,
evaluated at the (xi,eta) the destination ghost-face centre maps to.  Symmetric in the two
axes, so it hits swap and non-swap alike, and it SUBSUMES the along-seam resample.
Vertex stability is unaffected: this keeps real neighbour data, it only changes the
interpolation operator.

**Cheapest decisive test -- an ORACLE A/B, not a fix.**  In `srcval` add an env-gated arm
that replaces ONLY the partner component with the ANALYTIC iprob=11 field projected on
`PanelNormals`, everything else (transform, resample, index map) unchanged.  Read the
existing ERROR JUMP gate (`cs_test.cpp:2341`, printed `:2504-2523`), which already splits
per component into NONswap/SWAP, at nx=32/64/128.  If x3f NONswap and both SWAP jump to
~3, the co-location is the whole binding term; if they stay ~1.8, the transform or the
resample is binding.  **`nlim=0` IS valid for this one** -- it reads exchanged ghosts of
b0 and needs no resistive EMF -- and it does not touch the dead LOOP SPLIT instrument.

## THE REPO ALREADY PREDICTED BOTH FAILURES -- and my control was DEGENERATE

`cs_test.cpp:2343-2350` (the rationale block above the ERROR JUMP gate) states the paper
analysis: writing the discrete solution on panel P as B + h^2*phi_P with phi_P smooth, an
interior edge gives h*[e(x+h)-e(x)] = O(h^4) and J is 2nd order, while a SEAM gives
h*[e_A - e_B] = h^3*(phi_A - phi_B) and J is O(h)*dphi.  It then says explicitly:
**"Every cheap fix (higher-order interpolation, a partial-circulation exchange, averaging
the two panels) reproduces the SAME combination phi_A - phi_B and cannot move it."**
Fix 1 is higher-order interpolation and fix 2 is averaging the two panels.  Both failed,
exactly as written.  A second model independently re-derived the same thing as a PARITY
argument: the two-panel EMF average cancels every halo error that is EVEN under the seam
reflection (de-staggering, resample) and keeps only the ODD part (the face tilt / shear),
which is why fix 1 could not move the EVOLVED order no matter how much amplitude it
removed.  **Any future proposal must say what it does to phi_A - phi_B, or it is already
refuted.**

**BUT THE DECIDING MEASUREMENT HAS NOT ACTUALLY BEEN MADE.**  The gate's own design says
the question is the ORDER of the x1f (= B.rhat, chart-independent, no basis transform in
the way) seam jump: **3 means NO jump and the mechanism is something else; 2 means
dphi != 0**.  It requires a same-panel CONTROL at an ordinary block boundary, which
"needs more than one MeshBlock per panel tangentially to be non-degenerate".  EVERY run in
this session used `meshblock/nx2=nx3=` the full panel, i.e. ONE BLOCK PER PANEL, so the
control was EMPTY -- the logs say `same-panel CTRL L1=0.0000e+00 over 0 faces` and I read
straight past it.  The measured x1f seam jump was 2.80/2.90, which LEANS toward 3 (no
jump, so NOT the kink), but with no control it decides nothing.

**START HERE NEXT SESSION (a ~1 minute job).**  Re-run the seam-jump study with several
MeshBlocks per panel, e.g. `mesh/nx2=nx3=N` with `meshblock/nx2=nx3=N/2`, at N=32/64/128,
`nlim=0` (valid: this gate reads exchanged ghosts of b0 and needs no EMF).  Read the x1f
SEAM order against the now-live same-panel CTRL order.  Seam ~2 and control ~3 => dphi!=0,
the kink is real, and only a reflection-equivariant stencil can help.  Seam ~3 and control
~3 => there is no jump, the mechanism is elsewhere, and the halo work is back on.

## THE DECIDING MEASUREMENT, MADE 2026-09-02: **THERE IS NO JUMP.  dphi ~ 0.**

Ran step 1.  `inputs/tests/cubed_sphere_resist.athinput`, `mesh/nx2=nx3=N`,
`meshblock/nx2=nx3=N/2` (2x2 blocks per panel, 24 blocks) vs `N` (1 block, control
degenerate), N=32/64.  Serial `build_cs`, ~2 min each, run directly on viper13 -- NO
sbatch needed.  Working dir `/viper/u2/jinma/ATHENAK/bench/cs_seamjump`.

**FIRST ARM WAS VACUOUS AGAIN, AND FOR A NEW REASON.  `nlim=0` MAKES THE CONTROL
STRUCTURALLY ZERO.**  With 2x2 blocks the control finally had 4032 faces but read
`L1=0.0000e+00` EXACTLY.  Not a dead gate -- at `nlim=0` the discrete solution IS the
analytic projection, so phi_P = 0 identically; a same-panel ghost is a straight COPY of
the neighbour's analytically-set cell, so `|e_ghost - e_active| = |0 - 0|`.  **A
truncation-error JUMP cannot exist before the solution is evolved.**  The previous note's
"`nlim=0` IS valid for this one" is therefore WRONG for the control: it is valid only for
measuring HALO error.  Use `time/tlim=0.01`.  (The gate was alive: d1 and d2 agree
digit-for-digit on position bins u1/u2/u3, differing only at u0 where d2 skips `ng` cells
either side of the interior block boundary.)

**EVOLVED (`tlim=0.01`), the answer is seam ~3 AND control ~3:**

  x1f (B.rhat)  seam 2.5921e-07 -> 3.4577e-08   order **2.91**
                ctrl 7.5648e-09 -> 1.0478e-09   order **2.85**
  x2f           seam 4.5441e-06 -> 1.2647e-06   order 1.85
                ctrl 2.3937e-08 -> 3.5550e-09   order 2.75
  x3f           seam 4.0251e-06 -> 7.4470e-07   order 2.44
                ctrl 1.8840e-08 -> 3.3562e-09   order 2.49

By the gate's OWN stated criterion (`cs_test.cpp:2360-2364`: "3 means there is no jump and
the mechanism is something else, 2 means dphi != 0"), x1f is **3 on BOTH sides**.  The
seam carries ~34x the control's amplitude but the SAME order.  **So the chart-jump
phi_A - phi_B is NOT the binding term, the paper argument at `cs_test.cpp:2343-2350` does
NOT apply, and the "every cheap fix is already refuted" veto is LIFTED.**

**This puts the halo work back on.**  The tangential components say where: x2f and x3f
carry 100-200x the control's amplitude AND a lower order than their own control
(1.85/2.44 vs 2.75/2.49), which is the signature of HALO INTERPOLATION error dominating,
not a chart kink.  That is exactly what the UNTESTED `srcval` two-orthogonal-2-point-
averages hypothesis above predicts, and it is now the live lead.

**CAVEAT: two grid points only.**  N=128 (both decompositions) was relaunched detached in
`/viper/u2/jinma/ATHENAK/bench/cs_seamjump/ev_n128_d{1,2}/log.txt` and should be READ
FIRST next session to confirm 2.9 is not a two-point accident.

**How to apply:** the bench dir is `/viper/u2/jinma/ATHENAK/bench`, one level ABOVE the
repo -- searching inside `athenak/` will not find it.  See
[[cs-resistive-seam-order]], [[validate-the-instrument]].

## THE N=128 POINT LANDED 2026-09-02: the 2.91/2.85 reading was NOT a two-point accident

`bench/cs_seamjump/ev_n128_d2/log.txt`, baseline binary, evolved, 2x2 blocks per panel:

  x1f (B.rhat)  seam 2.5921e-07 3.4577e-08 4.3933e-09   order **2.91 / 2.98**
                ctrl 7.5648e-09 1.0478e-09 1.5048e-10   order **2.85 / 2.80**

Three points, seam and control both ~3 and the seam if anything CLEANER than the control
in the second interval.  **There is no chart jump; dphi ~ 0 is CONFIRMED**, and the veto
at `cs_test.cpp:2343-2350` stays lifted.  The halo lead this opened was then followed to
its end -- see [[cs-seam-colocation-fixed]].
