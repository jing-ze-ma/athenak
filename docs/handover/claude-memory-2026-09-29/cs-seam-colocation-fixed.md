---
name: cs-seam-colocation-fixed
description: The cubed-sphere seam halo's binding O(h^2) term was the PARTNER-COMPONENT co-location in bvals_fc srcval; an oracle A/B proved it and a clamped-window cubic fixes it, 5-11x on the evolved resistive field
metadata:
  type: project
---

2026-09-02. **CLOSED the seam-halo order problem that [[cs-seam-destag-refuted]] left open.**
Bench dir `/viper/u2/jinma/ATHENAK/bench/cs_oracle` (ABOVE the repo).

## The ORACLE A/B decided it in one run

Env-gated arm in `srcval` (`src/bvals/bvals_fc.cpp`) replacing ONLY the partner
tangential component with the ANALYTIC iprob=11 field at the primary's own (xi,eta);
transform, resample, index map, shear all untouched.  `nlim=0`, nx=32/64/128, serial
`build_cs`, ~20 s each.  **Arm 0 reproduced the recorded baseline DIGIT FOR DIGIT**
(1.8978e-07 / 5.2782e-06 / 7.1993e-06 at nx=32 -- instrument validated), and x1f, which
has no partner, was bitwise identical between arms -- an internal control.

  seam-jump L1     baseline (order)          oracle (order)
  x2f seam    5.28e-6 1.52e-6 4.06e-7  1.80/1.90   3.43e-7 4.39e-8 5.21e-9  2.97/3.08
  x2f SWAP    7.82e-6 2.26e-6 6.07e-7  1.79/1.90   4.19e-7 5.54e-8 6.56e-9  2.92/3.08
  x3f NONswap 1.29e-5 3.74e-6 1.01e-6  1.79/1.89   1.71e-6 2.07e-7 2.64e-8  3.05/2.97
  x3f SWAP    4.35e-6 1.29e-6 3.51e-7  1.75/1.88   3.37e-7 3.52e-8 4.20e-9  3.26/3.07

**Every category at 1.75-1.90 goes to 2.92-3.26, including BOTH SWAP ones that were
immune to CS_DESTAG4 and CS_SHEAR.**  Amplitude down 9-78x.  So the co-location is the
WHOLE binding term, not one of several.

## The fix, and the trap that made the first attempt fail

`srcval` co-locates the partner by two successive interpolations on ORTHOGONAL axes.
**Raising both to a CENTRED 4-point formula measured 1.80/1.90 -- exactly the baseline,
for a strictly more expensive stencil.**  The reason is structural and is the thing to
remember: **a seam buffer packs the cells at the EDGE of the active range**, so in the
seam-NORMAL direction a centred wide stencil NEVER FITS and silently degrades to the
2-point average over the whole packed range.  (The shear correction's comment already
said this about its own stencil; I did not connect it.)

What works is a fourth-order interpolation whose WINDOW SLIDES INWARD instead of
degrading: `lag4`, the cubic through four consecutive samples with the window CLAMPED
into the active range, evaluated at the target offset.  It never reads the source
block's own ghosts, stays 4th order to the edge (one-sided, not low order), and
extrapolates with the same window at the two seam ends.  Falls back to the old 2-point
form when an axis has fewer than 4 active samples (nx or cnx < 4) -- same degeneracy
guard as [[cs-narrow-block-resample-degeneracy]].

**It lands ON the oracle**: seam-jump within 1.5-3% of the oracle in every category
(x2f 3.54e-7/4.51e-8/5.30e-9, x3f 8.05e-7/9.31e-8/1.16e-8), i.e. the halo now reaches
the ceiling the oracle set.

## What it buys on the EVOLVED resistive field (tlim=0.01, the live gate)

                    nx=32              nx=64          gain
  global   L1   2.26e-6 -> 7.68e-7   5.52e-7 -> 1.91e-7   2.9x
  global   Linf 1.93e-4 -> 2.90e-5   6.87e-5 -> 1.15e-5   6.0x
  INTERIOR L1   7.41e-7 -> 5.99e-7   2.87e-7 -> 1.58e-7   1.8x
  SEAM     L1   5.25e-6 -> 1.19e-6   2.03e-6 -> 3.92e-7   5.2x
  VERTEX   L1   2.74e-5 -> 2.41e-6   1.07e-5 -> 9.34e-7  11.5x

The **CUBE VERTEX gains the most (11.5x)** -- the region [[cs-seam-destag-refuted]] left
"still unexplained" and that [[cs-mhd-validation]] measured at order 0.65.  Seam is now
2.5x the interior (was 7.1x) and the vertex 5.9x (was 37x).  Ohmic heating is unchanged
to 4 digits (1.003166 vs 1.003170) and div B is untouched, so the physics gates hold.

**Why:** the halo work that [[cs-seam-destag-refuted]]'s "no chart jump, dphi ~ 0" reading
put back on the table is now done, and the paper veto at `cs_test.cpp:2343-2350` is
confirmed inapplicable -- a cheap higher-order interpolation DID move it, once it was the
right one.

**How to apply:** the change is inside `do_cs`, so every non-cubed-sphere grid is bitwise
unchanged by construction.  Read the ORDERS at fixed CFL with suspicion --
[[cs-mhd-validation]] records that the radial-BC clock makes the evolved gate read ~1st
order everywhere; the AMPLITUDE ratios above are the signal.  See
[[cs-resistive-seam-order]], [[validate-the-instrument]].

## It also helps ACROSS A LEVEL BOUNDARY (SMR), and the narrow-block case is inert

A/B against a CPU baseline binary built from the pre-change HEAD in a git WORKTREE at
`/viper/u2/jinma/ATHENAK/bench/base_wt` (the bench `athena` is a GPU build and cannot run
on the login node -- that is why a fresh baseline was needed; note the worktree has no
`kokkos`, symlink the main one in or cmake fails).

`cubed_sphere_resist_smr` (cross-level, evolved):

                base -> new        gain
  global L1   1.3024e-05 -> 1.1070e-05   1.18x
  global Linf 5.4221e-04 -> 3.5288e-04   1.54x
  SEAM   L1   1.4823e-05 -> 1.0865e-05   1.36x
  VERTEX L1   4.9385e-05 -> 1.3853e-05   3.57x
  INTERIOR    1.1245e-05 -> 1.1051e-05   (unchanged, as it must be)

**The seam is now BELOW the interior** (1.087e-05 vs 1.105e-05); it was 1.32x above.
Gains are smaller than on the unrefined grid because the cross-level path has its own
error sources, but they are all in the right direction.

`cubed_sphere_resist_smr_narrow` is essentially INERT (interior 1.0169e-04 -> 1.0179e-04,
Ohmic 1.079615 -> 1.079135, ghost-scan counts 622/324 -> 621/325), which is the expected
signature: its blocks are too narrow for a cubic window, so the fallback keeps the old
form.  Its large ghost errors (6.2e-02, hundreds of cells over 1e-2) are PRE-EXISTING and
unchanged -- do not read them as a regression from this change.

## THE ORDER DID *NOT* MOVE -- the gain is a CONSTANT.  Measured, not assumed.

Two corrections to what I started from, both measured (`/viper/u2/jinma/ATHENAK/bench/cs_order`,
12 arms: {base,new} x nx={32,64} x cfl={0.3,0.15,0.075}, run in PARALLEL -- sequential
would have taken ~2 h):

**1. This gate is dt-INDEPENDENT, so no dt->0 extrapolation is needed.**  Over a 4x change
in dt the L1 moves 0.01% at nx=32 (2.2575 -> 2.2578e-06) and 0.005% at nx=64
(5.5242 -> 5.5245e-07).  The "at fixed CFL it reads 1st order everywhere from the
radial-BC clock" caveat in [[cs-mhd-validation]] is about the RIGID-ROTATION test
(iprob=9, tlim=1.0) and does NOT apply to this static iprob=11 gate at tlim=0.01.  The
fixed-CFL orders ARE the spatial orders here.

**2. The baseline seam order is 1.37/1.67, NOT 1.0**, and it TRENDS UP with resolution
(nx=32/64/128, three points).  The recorded "seam order 1.00" comes from the LOOP-SPLIT
input instrument, not from this evolved gate -- they are different numbers and I had
been conflating them.

  region     BASE 32/64/128            order      NEW 32/64/128            order    gain@128
  global   2.257e-6 5.524e-7 1.360e-7  2.03/2.02  7.681e-7 1.910e-7 4.743e-8  2.01/2.01   2.9x
  INTERIOR 7.410e-7 2.865e-7 9.650e-8  1.37/1.57  5.991e-7 1.582e-7 4.201e-8  1.92/1.91   2.3x
  SEAM     5.245e-6 2.034e-6 6.400e-7  1.37/1.67  1.190e-6 3.917e-7 1.208e-7  1.60/1.70   5.3x
  VERTEX   2.738e-5 1.074e-5 3.774e-6  1.35/1.51  2.406e-6 9.339e-7 3.030e-7  1.37/1.62  12.5x

ALL THREE RESOLUTIONS ARE IN.  Global is a clean 2.01/2.01.  The REGION orders TREND
UPWARD with refinement (seam 1.60 -> 1.70, vertex 1.37 -> 1.62), which is the signature of
a fixed-CELL region creeping physically closer to the seam -- they approach the true rate
from below rather than sitting at a defect.  See [[cs-seam-order-limiter]] for why an
order must not be read off these regions at all.  Ohmic heating at nx=128: 1.000725 vs
the baseline's 1.000726.

**So the halo fix buys a 5-12x CONSTANT, not an order.**  Something else limits the
evolved seam and vertex order to ~1.5-1.7, and the halo is no longer it -- the halo now
converges at ~3, well clear of the operator.  That is the next thing to find, and this
table is the baseline to beat.  (INTERIOR going 1.37 -> 1.92 is a real and separate
improvement: seam-halo error had been leaking into the interior norm.)

**Do not repeat the mistake this corrects:** I nearly wrote up the 5-11x as an order
improvement.  Measure the order on the SAME instrument you quote the amplitude from.
