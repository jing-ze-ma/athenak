---
name: cs-seam-order-limiter
description: After the co-location fix NOTHING limits the cubed-sphere seam -- J at the seam sits on the interior's operator floor and the global L1 is 2.01 (2.08 eta-independent, 1.99 resistive); the radial floor and the '~1.8 residual' are both REFUTED
metadata:
  type: project
---

2026-09-02, continuing [[cs-seam-colocation-fixed]].  Bench dir
`/viper/u2/jinma/ATHENAK/bench/cs_limit` (ABOVE the repo).  Question asked: after the
halo fix, what still limits the evolved seam order (region L1 read 1.60)?

## ANSWER: nothing seam-specific.  The seam now matches the panel interior.

**The LOOP SPLIT instrument settles it** (`cs_test.cpp:2760`, run with `time/nlim=1`; it
separates the resistive OPERATOR's own consistency from what the INPUT errors contribute,
on a 2-cell seam ring).  nx=32/64/128, both binaries:

  on the seam RING          32         64        128      order
  base input          7.135e-04  3.437e-04  1.686e-04   1.05/1.03
  NEW  input          5.077e-05  1.169e-05  2.665e-06   **2.12/2.13**
  base total          8.673e-04  3.799e-04  1.767e-04   1.19/1.10
  NEW  total          2.522e-04  6.218e-05  1.535e-05   **2.02/2.02**
  OPERATOR (both)     2.480e-04  6.269e-05  1.578e-05   1.98/1.99
  interior total      2.143e-04  5.422e-05  1.364e-05   1.98/1.99

The new ring TOTAL has fallen onto the OPERATOR floor, and that floor equals the
interior.  **J at the seam is now the same quality as J in the interior**, and the
"ring by distance to the nearest VERTEX" profile went from a 9.5x gradient
(3.61e-3 -> 3.80e-4) to FLAT (3.71e-4 ... 2.48e-4).  The vertex concentration in the
current is gone too.

## The evolved-field error has no seam structure left

Use the `BY DISTANCE from the panel edge (cells)` profile the evolved gate already prints
(`cs_test.cpp:2082`) -- it is the honest instrument, not the REGION split.

  new, nx=32   d0=1.625e-06 d1=8.431e-07 d2=7.506e-07 d3=7.291e-07 d4=6.733e-07 d5=5.243e-07
  new, nx=64   d0=4.808e-07 d1=3.259e-07 d2=2.327e-07 d3=1.983e-07 d4=1.934e-07 d5=1.464e-07
  order          1.76         1.37         1.69         1.88         1.80         1.84

Baseline for contrast had a MID-DISTANCE PLATEAU converging at 0.28-0.53 (d2/d3/d4) --
halo error diffusing outward -- which is now gone.  And **the deep-interior floor is
IDENTICAL between arms** (d5 = 5.252e-07 base vs 5.243e-07 new at nx=32): the fix does
not touch it, because it is not a seam error.  The profile is now essentially FLAT from
d1 outward, with no concentration at the seam.  (The per-bin "orders" in the last row are
NOT trustworthy as orders -- see the retraction below -- but the SHAPE of the profile at
a given resolution is exactly what this instrument is for, and it says the seam excess is
gone.)

## REFUTED: the radial floor.  It contributes ~2%.

I suspected the study's **fixed `nx1 = 8`** (every run refines only nx2/nx3) was a
non-converging floor -- [[cs-hydro-validation]] warns that a tangential-only refinement
gives a FALSE first-order reading.  **Wrong here.**  Refining ONLY the radial direction
at fixed nx2=nx3=64 (`R_8/R_16/R_32`):

  nx1     global L1     d0        d5       SEAM
  8     1.9097e-07  4.808e-07 1.464e-07  3.9166e-07
  16    1.8719e-07  4.656e-07 1.443e-07  3.8157e-07    -2.0%
  32    1.8366e-07  4.573e-07 1.416e-07  3.7442e-07    -1.9%

**A 4x radial refinement moves the total by 3.8%, and it is SATURATING** (-2.0% then
-1.9% per doubling), so the radial share of the error is at most ~4%.  Every bin moves
by the same ~2-3%, seam and interior alike.  The error is tangential through and through.
`R_8` reproduced the nx1=8 run exactly, so the sweep was a valid instrument.

## A MEASUREMENT CORRECTION worth keeping

**The global L1 order (2.01) is FLATTERED, not honest.**  The bad-region fraction
SHRINKS as the grid refines -- the d5 bin holds 47% of faces at nx=32 but 71% at nx=64 --
so the domain average converges faster than any individual bin in it.  The code comment
at `cs_test.cpp:1959` warns about this effect for the seam being hidden; it also runs the
other way and inflates the global rate.  **Quote the per-distance bins, not the global.**

## ANSWERED: BOTH components of the error are 2nd order.  Nothing is limiting it.

The eta sweep finished at BOTH resolutions.  Fit L1 = A + B*eta from eta=0.5 and 1.0
(both diffusion-limited, so like-for-like):

  nx=32   A(eta-indep) 1.7814e-07 (23%)   resistive@eta=0.5 5.8993e-07 (77%)
  nx=64   A(eta-indep) 4.2260e-08 (22%)   resistive@eta=0.5 1.4871e-07 (78%)

  order of the eta-INDEPENDENT part  **2.08**
  order of the RESISTIVE part        **1.99**
  order of the total at eta=0.5      **2.01**

The split is STABLE with resolution (23% -> 22%), and **both physical components converge
at 2.0**.  There is no anomalous sub-second-order component in the global norm.

## THEREFORE: I OVER-CORRECTED EARLIER.  The global L1 is HONEST; the REGION and
## PER-DISTANCE "orders" are the artefact.

My earlier reading in this file -- "the global 2.01 is flattered, quote the per-distance
bins, the field really converges at ~1.8" -- is **WRONG and is retracted**.  The reason is
that EVERY binned region here is defined in CELLS, so it is NOT THE SAME PHYSICAL REGION
across resolutions:

* d0..d4 are fixed-cell shells that shrink physically as the grid refines;
* **d5 is "everything at distance >= 5 CELLS", which GROWS to swallow cells that are
  physically much closer to the seam** as nx doubles.  Since the error rises toward the
  seam, d5's average is inflated at high nx and its apparent order (1.84) is DEPRESSED.

The global L1 is a norm over a FIXED physical domain and is the one number whose order
means what it says: **2.01**, and 2.08 / 1.99 in its two physical parts.

**Use the region and per-distance splits to LOCALISE error at a given resolution -- that
is what they are for and they did their job here.  Do NOT read an order off them.**
The one exception is d0, "the cell adjacent to the seam", which is a genuine pointwise
quantity; its 1.76 is meaningful and is the only sub-2 number left standing.

## THE eta=0 ARM IS VACUOUS -- do not repeat it

`mhd/eta_ohm_const=0` does NOT give a clean "ideal" control here, for TWO reasons at once:

* **dt explodes.**  Without the resistive diffusion limit dt goes 3.485e-05 -> 8.749e-03,
  so the run reaches tlim=0.01 in **2 cycles** instead of ~290.  Nothing is comparable.
* **The gate disappears.**  `CSTestResistCheck`'s output is tied to the resistivity
  module, so with eta=0 the run prints NO `EVOLVED FIELD` / `REGION` / `BY DISTANCE`
  lines at all.  The log looks like a normal completed run.

**Corrected recipe for next time:** use a TINY but nonzero eta (e.g. 1e-6) so the module
and its gate stay alive, and PIN dt across every arm (a fixed `time/cfl_number` small
enough that the eta=1.0 arm is not diffusion-limited, or an explicit dt cap) so all arms
take the same number of steps.  Only then is L1(eta) a like-for-like fit.  The eta=0.5
and eta=1.0 arms above ARE comparable to each other (both diffusion-limited, dt 3.485e-05
and dt halved), which is what makes the 23%/77% split usable.

**How to apply:** do NOT go looking for another seam bug -- the seam is done, and this
file is the evidence.  See [[cs-seam-colocation-fixed]], [[validate-the-instrument]],
[[measure-impact-before-claiming]].
