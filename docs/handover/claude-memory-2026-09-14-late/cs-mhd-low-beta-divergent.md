---
name: cs-mhd-low-beta-divergent
description: THE REPRODUCER -- the cubed-sphere MHD scheme DIVERGES under refinement for plasma beta < 1, everywhere, and 24-36x worst at the cube vertex; a 3-second CPU test with no gravity, RT or stratification
metadata:
  type: project
---

**2026-09-02. This is the minimal reproducer of the cs+MHD blow-up.**

`cs_test` iprob=9 -- the rigid-rotation equilibrium carrying a UNIFORM (hence force-free,
curl B = 0, so J x B = 0 at ANY beta) Cartesian field -- has an exact solution at every
time and already reports L1 split by panel INTERIOR / SEAM / CUBE VERTEX. Every previous
cs MHD validation ran it at the default `p0 = 1, b0c = 1`, which is **beta = 2**. The
failing hot-Jupiter atmosphere sits at **beta ~ 1e-3**. Nobody had ever run it there.

Sweep beta by `b0c` at fixed `p0 = 1` (beta = 2 p0 / b0c^2), with `tlim` scaled as 1/b0c
so every beta does a comparable number of Alfven crossings. CONVERGENCE ORDER of L1(B):

    beta   nx        interior    seam       VERTEX
    2.0    16->32      +1.42      +1.28      +1.06
    2.0    32->64      +1.25      +1.20      +1.03     <- converges, the only tested regime
    0.1    16->32      -0.75      -1.24      -1.46
    0.1    32->64      -0.90      -0.85      -0.86     <- DIVERGES
    0.01   32->64      -0.80      -0.86      -0.88
    0.001  32->64      -0.77      -0.86      -0.88

**For beta < 1 the scheme is DIVERGENT: refining makes the error grow, roughly as
h^-0.85.** Vertex-to-interior error ratio at nx=64: 2.8x at beta=2, but **36x at beta=0.1**
and 24x at beta=1e-3.

## Why this explains everything

* **Refining makes it worse** -- matches the production run exactly (nx=32 dies at 0.1506
  rot, nx=16 at 0.1901). Same sign, and it rules out under-resolution.
* **Scales with |B|** -- beta falls as bbot^2, which is the bbot threshold between 1 and
  3 G ([[cs-mhd-dhj-blowup]]).
* **Vertex-localised** -- 24-36x worse there, which is why ignition is at a cube vertex.
* **MHD-only** -- hydro has no beta and is clean past 74 rotations.
* Consistent with dt-independence, round-off divB, and an exact halo: a divergent
  OPERATOR needs none of those to be broken.

## RETRACTED, then settled: it IS confined to the vertex

I first read the table above as "the interior diverges too (-0.8), so this is a low-beta
failure EVERYWHERE". **That was wrong, and it was contamination.** At the sample times
above the vertex is already O(1) destroyed (at beta = 0.1, L1(p) at the vertex is 33-50%
of p0 and Linf(p) ~ 170%), and a destroyed vertex radiates into the interior long before
tlim.

Re-measured at EARLY time (tlim/10 = 0.002236 at beta = 0.1, so the sound crossing covers
~0.003 of the domain and the vertex CANNOT have reached the interior):

    nx     L1(B) interior     L1(B) VERTEX
    16       1.8022e-05        2.1199e-04
    32       9.1670e-06        5.9993e-04
    64       6.7244e-06        2.0605e-03

**The interior CONVERGES (1.8e-5 -> 6.7e-6); the vertex DIVERGES (2.1e-4 -> 2.1e-3, 10x
over two refinements).** So it is a vertex defect in B itself, and the "everywhere"
reading is retracted. **Always separate a diverging region from its downstream
contamination by sampling EARLY**, or a local bug reads as a global one.

## A second mechanism, distinct from the bug: catastrophic cancellation

At beta << 1, `p = (gamma-1)(E - KE - ME)` is a small difference of large numbers -- at
beta = 0.01, ME/p ~ 100 -- so a ~1% error in B becomes a ~100% error in p. Measured at
beta = 0.01: L1(p) reaches 6.6 against p0 = 1 while L1(B)/b0c is only ~1%. This is
INHERENT to energy-based conservative MHD, not a cs bug, and it is the amplifier that
turns the vertex's B error into a blow-up. The interior's own p error at low beta is
muddier and I did not resolve whether it is truncation in E amplified by 1/beta or
something else -- do not quote it as a result.

## How to run it (3 seconds, CPU, no GPU needed)

    build_cs/src/athena -i cubed_sphere_mhd_conv.athinput -d . problem/iprob=9 \
      mesh/nx2=32 mesh/nx3=32 meshblock/nx2=32 meshblock/nx3=32 \
      problem/b0c=4.472 problem/p0=1.0 time/tlim=0.02236

Sweep script `sweep.sh` in the session scratchpad. Both instruments are COMMITTED as
**22e2f712** (the FV `mhd_divb` and `CSCornerHaloSmoothness`); neither changes an answer. **Any fix must be gated at beta <= 0.1,
not at the default beta = 2, which hides it completely.** Cf. [[validate-the-instrument]].


## LOCALISED: it is the PLM RECONSTRUCTION at the vertex, and it is beta-dependent

Four measurements, all `nlim=1` (one timestep, nx=32) so the error is the OPERATOR's, not
an accumulated instability. **After ONE step the vertex error is already 66x the
interior** (6.54e-4 against 9.90e-6) -- so this is debuggable at one step, essentially
free.

**1. The one-step vertex error is beta-dependent, with a threshold at beta ~ 1:**

    beta     200      2       0.1      0.01     0.001
    vtx/int  1.7      1.7     66.1     57.8     53.3

**2. The scalings differ.** From beta=2 to beta=0.1 the field grows 4.47x; the INTERIOR
error grows 9x (~ B^2, the expected magnetic scaling) but the VERTEX error grows 350x
(~ B^3.5). Something at the vertex carries a higher power of B than it should.

**3. It is the RECONSTRUCTION.** Same test, varying only `mhd/reconstruct`:

    variant        beta=0.1 vtx/int     beta=2 vtx/int
    donor cell           5.9                  5.1      <- beta-INDEPENDENT
    plm                 66.1                  1.7      <- beta threshold appears

With piecewise-constant reconstruction the vertex/interior ratio does not care about beta
at all. With PLM it explodes only once magnetic pressure dominates. **So the beta-dependent
vertex anomaly is in the PLM reconstruction path, not in the fluxes, the EMF or CT.**
(`wenoz` fails to run in this setup -- more ghost zones needed -- not pursued.)

**4. Two vertex suspects EXONERATED by direct A/B, both to 4+ digits:**
   * the cube-vertex corner FILL (`mesh/cs_vertex_fill` true vs false): 2.1199e-04 vs
     2.1190e-04, and bit-identical at nx=32 and 64;
   * the three-way vertex EMF average (`AveragePanelCornerEMF`, skipped via a temporary
     getenv switch, since REVERTED): bit-identical at every resolution. Verify such a
     switch is really compiled in (`strings` the binary) before trusting a null result --
     identical-to-all-digits is exactly what a no-op edit also produces.

## RETRACTED: "it is the FIELD reconstruction, not the fluid one"

**The donor-cell substitution method is INVALID for this question, and every conclusion
drawn from it is withdrawn.** DC does not REDUCE the low-beta vertex error; it raises the
high-beta baseline, which is what makes both the ratio and the growth factor fall.
Absolute vertex L1(B) at beta = 0.1:

    PLM all          6.54e-04
    B.rhat -> DC     8.34e-04     <- WORSE
    B.e_xi -> DC     1.02e-03     <- WORSE
    FLUID  -> DC     6.42e-04     <- same

DC's own vertex error is comparable to the anomaly being hunted, so substituting it cannot
isolate the anomaly. Both the "ratio" metric (denominator rises) and the "growth from
beta=2 to beta=0.1" metric (numerator baseline rises) are confounded the same way.
**A substitution experiment only discriminates when the substitute is MORE accurate than
the effect you are chasing.** Cf. [[measure-impact-before-claiming]],
[[confounded-tests-rejected]].

The ONE result from that family that survives is a NULL: DC on B.f2 alone leaves the
answer identical to 5 digits (6.5439e-04 vs 6.5438e-04), and identical means B.f2's
reconstruction has no influence on the answer at all -- so **the B.f2 index-space kink,
real as it is (27x on a varying field), is NOT the cause.** A null from substitution is
safe; a positive is not.

## SUPERSEDED (kept for the trail): it is the FIELD reconstruction, not the fluid one

`mhd_fluxes.cpp` reconstructs the fluid and the field with two SEPARATE calls
(`PiecewiseLinearX*(... w0_ ...)` and `PiecewiseLinearX*(... b0_ ...)`). Patching them
behind two env switches (since REVERTED) and running one timestep at nx=32:

    variant                     beta=0.1 vtx/int    beta=2 vtx/int
    PLM both                        66.1                1.7        <- the bug
    FIELD donor cell, FLUID plm      6.0                4.8        <- beta dependence GONE
    FLUID donor cell, FIELD plm     44.4                2.1        <- beta dependence REMAINS

**Making only the FIELD reconstruction first-order removes the beta dependence entirely.**
So the defect is in PLM applied to `b0_` (the cell-centred field) near the cube vertex.

## What was CHECKED and is CORRECT (do not re-derive)

The per-sweep reconstruction FRAMES are right. bcc holds the orthonormal triple
(B.rhat, B.e_xi, B.f2) with f2 = (e_eta - c e_xi)/s. Since e_xi.f2 = 0 exactly:

  * x1 frame {rhat, e_xi, f2}      -> face-parallel slots 1,2 already correct;
  * x3 frame {f2, rhat, e_xi}      -> face-parallel slots 0,1 already correct;
  * x2 frame {nhat_xi, e_eta, rhat} -> needs B.e_eta = c B.e_xi + s B.f2, which is exactly
    what `GnomonicEquiangleFaceBX2` computes.

So x1 and x3 genuinely need no field rotation, and the "no rotation needed" comments in
those two sweeps are correct -- that asymmetry is NOT the bug. The trig pairs are right
too (x1 -> cell, x2 -> xi face, x3 -> eta face), matching hydro.

One measured non-lead: sin_cell in the CORNER GHOSTS degrades only mildly (0.88 at the
corner cell to ~0.82 two ghosts out at nx=32), so 1/s does not blow up there; a naive
"B.f2 = (...)/s diverges in the corner ghost" story does NOT hold up numerically.

## THE DEFECT: the corner-halo `bcc` third component B.f2 has an INDEX-SPACE KINK

New gate in `cs_test` CSTestGhostCheck under iprob=8 (`### CS CORNER-HALO bcc SMOOTHNESS
IN INDEX SPACE`). It needs no analytic solution: smoothness in index space is a property
of the SEQUENCE, so it compares |d2 f| = |f(n-1) - 2f(n) + f(n+1)| along the j index
ACROSS the active->ghost crossing at a cube vertex against the same stencil deep inside
the panel. nx=32, b0c=4.472:

    component   d2 at CORNER crossing   d2 deep INTERIOR   ratio
    B.rhat            2.6e-03               3.6e-03          0.7
    B.e_xi            6.8e-03               6.7e-03          1.0
    B.f2              7.8e-07               4.4e-16       ~1e9

**B.rhat and B.e_xi are perfectly smooth across the corner. B.f2 -- and only B.f2 -- has
a kink nine orders of magnitude above its own interior smoothness.**

**And B.f2 is the only component that involves the metric.** bcc stores
(B.rhat, B.e_xi, B.f2) with f2 = (e_eta - c e_xi)/s, so
B.f2 = (B.e_eta - c B.e_xi)/s carries sin_cell and cos_cell; the other two do not.
The two trig-free components are clean; the trig-carrying one is broken at the corner.
That is why only the FIELD reconstruction misbehaves (the fluid never uses f2), why it is
localised at the vertex (where the trig is most extreme), and why PLM matters -- PLM
limits second differences, which is exactly the quantity that jumps.

### It is UPSTREAM of bcc, in the halo FACE field, and it is SEAM-WIDE

Extending the same gate to the face fields and the trig, on the same index line:

    quantity   d2 at corner crossing   d2 deep interior
    b0.x2f          7.47e-03               6.76e-03      <- SMOOTH (ratio 1.1)
    b0.x3f          1.54e-05               4.44e-16      <- kinked
    sin_cell        2.63e-04               1.45e-06      <- ratio 181
    cos_cell        1.15e-03               1.45e-06      <- ratio 790

**The kink is already in the halo FACE field b0.x3f, before bcc is formed** -- so bcc.f2
inherits it rather than creating it, and the TRIG ARRAYS ARE NOT THE CULPRIT (181-790x is
ordinary extra curvature near a corner, not a jump). The component that is kinked is the
one TANGENTIAL to the crossing; the normal one (b0.x2f) is smooth.

And it is **seam-wide, worst at the ends**. Same j crossing, three along-seam positions:

    k = ks (at a vertex)   d2 = 1.54e-05
    k = mid-seam           d2 = 1.88e-06
    k = ke (other vertex)  d2 = 2.82e-05

### RETRACTED: the "9 orders of magnitude" framing

The huge ratios are inflated by a DEGENERATE DENOMINATOR. For this uniform field on this
grid, b0.x3f is EXACTLY LINEAR in j in the panel interior, so d2_interior = 4.4e-16 is
essentially zero rather than a normal truncation scale. In ABSOLUTE terms the kink is
1.9e-06 mid-seam and 1.5-2.8e-05 at the vertices against |B| = 4.14 -- i.e. **4e-7 to
7e-6 relative**, which is SMALL, and smaller than a generic O(h^2) term. The structure
(which component, where, seam-wide, worse at the ends) is established; **the MAGNITUDE is
NOT, and this degenerate test cannot establish it.** Do not quote the ratios as a measure
of severity. Cf. [[measure-impact-before-claiming]].

### The gate on a VARYING field (iprob=11) -- the real numbers

The gate is now a standalone `CSCornerHaloSmoothness(pmbp)` called from BOTH
`CSTestGhostCheck` (iprob 8) and `CSTestResistCheck` (iprob 11). iprob=11's field
B = b0c*(-y,x,0) + tilt VARIES, so the interior d2 is a meaningful reference rather than
the degenerate ~0 of the uniform case:

    component   d2 corner    d2 interior   ratio
    B.rhat       5.8e-05      8.1e-05       0.7
    B.e_xi       1.6e-04      1.5e-04       1.1
    B.f2         5.6e-05      2.1e-06      27.3   <- on 4 of 6 panels (0.8 on 3 and 5,
                                                     whose interior reference is 40x larger)
    b0.x2f       6.7e-04      1.5e-04       4.5
    b0.x3f       3.7e-05      4.1e-06       9.0

**B.f2 remains the unique standout at 27x while the two trig-free components sit at ~1x.**
So the STRUCTURE is confirmed on a non-degenerate problem. The magnitude is 27x, not the
1e9 the uniform test appeared to show -- and in absolute terms 5.6e-05 against |B| ~ 1.

**HONEST STATUS: the chain is not closed.** A 27x elevated second difference in one bcc
component is consistent with, but has NOT been shown to produce, the 66x flux error at
beta = 0.1. What IS established end to end: (a) the failure is a cs bug, (b) it is in the
FIELD reconstruction (field->donor-cell removes the beta dependence, fluid->donor-cell
does not), (c) B.f2 is the only bcc component with an anomalous index-space second
difference at a seam/corner, and (d) B.f2 is the only one built through the metric
(c and s). Closing the chain needs the amplitude of the PLM slope error at the vertex
compared against the observed flux error -- do not assert causation before that.

**Superseded suspect**: the sin_cell / cos_cell used when `bcc` is formed in the CORNER GHOST cells
(`GnomonicEquiangleRaiseVelMHD`). The halo FACE fields are exact there (verified), so if
the trig used to convert them to the orthonormal triple is not the true chart
continuation into the corner ghosts, B.f2 alone picks up a kink. Check what
sin_cell/cos_cell hold in the (both-indices-ghost) corner block, since they are 3D arrays
(m,k,j) with no radial extent and are filled by CoordGnomonicEquiangle.

CAVEAT on magnitude: in the UNIFORM-field test B.f2 is nearly constant, so the absolute
kink (7.8e-07 against |B| ~ 4.1) is tiny -- the test is degenerate in amplitude. What it
establishes is the STRUCTURE (which component, where), not the size; in a varying field
the kink should scale with the field's own variation.

## WHAT STILL STANDS, AND WHAT WOULD ACTUALLY WORK

Solid, because measured with PLM throughout and no scheme substitution:
 * beta > 1 converges, beta <= 0.1 does not, threshold near beta = 1;
 * sampled EARLY, the panel interior CONVERGES while the cube vertex DIVERGES;
 * after ONE timestep the vertex error is already 66x the interior;
 * the vertex error grows ~B^3.5 where the interior grows ~B^2;
 * B.f2 has a 27x index-space kink -- and is NOT the cause (clean null above).

NOT established: which reconstruction, or which component, carries it.

**The method that would actually answer it** is a DIRECT measurement, not a substitution:
compare the RECONSTRUCTED face states (wl/wr, bl/br) against the ANALYTIC face values at
the vertex, per component, on iprob=9 where the exact solution is known everywhere. That
measures the reconstruction's own error where it happens, with no substitute whose
accuracy has to be assumed. Do that before touching any operator.

## THE DIRECT MEASUREMENT: the HALO IS FINE, so the defect is in the OPERATOR

Done at last (gate `CSCornerHaloSmoothness`, "DIRECT" section): compare the STORED bcc in
the ghosts against the ANALYTIC value -- for a uniform Cartesian field, B_cart projected
on the chart-continued basis (rhat, e_xi, f2) at the ghost's own (xi,eta).

**Only B.rhat can be resolved this way.** The reference is a CELL-CENTRE projection while
bcc is a FACE AVERAGE, so the two differ at O(h^2); the active interior already shows
1.9e-04 in B.e_xi and 2.1e-04 in B.f2, which swamps any halo signal in those two. But
B.rhat is EXACTLY 0 in the active interior, so for that component the reference is exact.

    B.rhat error in the ghosts (b0c = 1, |B| = 1)
    nx      vertex ghost    mid-seam ghost    vtx/mid
    16       2.05e-04         1.21e-06          170
    32       3.01e-05         9.81e-08          307
    64       4.09e-06         6.82e-09          600

**It CONVERGES: order ~2.8 at the vertex, ~3.8 mid-seam.** Larger at the vertex, but
convergent -- so this is not a non-convergent halo defect either. (It also scales exactly
linearly with |B|, 3.0e-05 -> 1.35e-04 for a 4.47x field, so it carries no beta anomaly.)

**The deduction that matters:** the halo is convergent, yet the SOLUTION at the vertex
DIVERGES with refinement at low beta. So the defect is NOT in the halo inputs at all --
it is in the interior OPERATOR applied at the skew vertex cells, fed with correct halo
data. Every remaining suspect is on the operator side.

## THE BUG, CHARACTERISED: the vertex operator is INCONSISTENT below beta ~ 0.5

The decisive measurement is the PER-STEP error (nlim=1) against resolution. That separates
the operator's own local truncation error from any accumulation, and it is unambiguous.

**Per-step L1(B) convergence order at the CUBE VERTEX, against beta:**

    beta      nx16        nx32        nx64        vtx order   interior order
    200     8.93e-07    1.55e-07    3.30e-08      +2.23          +2.88
    2       1.01e-05    1.87e-06    4.58e-07      +2.03          +2.60
    0.5     2.55e-05    4.56e-06    1.17e-06      +1.96          +2.60
    0.1     6.35e-04    6.54e-04    6.71e-04      **-0.04**      +2.09
    0.01    3.55e-03    3.60e-03    3.66e-03      **-0.02**      +2.00

**There is a SHARP transition between beta = 0.5 and beta = 0.1.** Above it the vertex is
cleanly second order; below it the per-step error STOPS CONVERGING -- order zero, the local
truncation error no longer vanishing as h -> 0. **That is an INCONSISTENT operator**, which
is a far stronger statement than "inaccurate": no amount of refinement helps, and refining
makes the accumulated solution worse (the -2 order seen over many steps).
The panel INTERIOR stays second order at EVERY beta.

**All three variables carry it.** At beta = 0.01, one step, nx 16 -> 32 -> 64:
L1(v) 0.110 -> 0.128 -> 0.138, L1(p) 5.53 -> 5.98 -> 6.20, L1(B) 3.55e-3 -> 3.60e-3 ->
3.66e-3. The B error is only 2.5e-4 RELATIVE to |B| = 14.14 -- small, but non-convergent,
and the pressure is 620% wrong because at low beta p = (gamma-1)(E - KE - ME) amplifies it.

**The Riemann solver is NOT the cause.** llf, hlle and hlld all diverge at the same rate
over many steps (order -1.98, -1.99, -1.78 at beta = 0.1). hlld is merely a constant ~1.7x
worse in absolute terms at low beta while being the BEST at beta = 2 -- a factor, not the
mechanism. So it is upstream of the solver.

**The EMF-rotation asymmetry is NOT the cause either** (x1 rotates its EMFs, x2/x3 do not).
Verified correct: CT needs E on the EDGE directions; the x2 frame {nhat_xi, e_eta, rhat}
and x3 frame {nhat_eta, rhat, e_xi} have both face-parallel axes ON edge directions, while
x1's frame {rhat, e_xi, f2} has f2, which is not an edge -- hence
E.e_eta = c*e21 + s*e31 in x1 alone. Correct as written.

## THE FLUX ROTATION IS NOT IT -- all three verified algebraically

Derived the required index lowering independently and compared, using
f2.e_xi = 0, f2.e_eta = s, nhat_xi.e_xi = s, nhat_xi.e_eta = 0, nhat_eta.e_xi = 0,
nhat_eta.e_eta = s:

  * X1, frame {rhat, e_xi, f2}:      need m_2 = A, m_3 = c A + s B.  Code matches.
  * X2, frame {nhat_xi, e_eta, rhat}: need m_2 = s N + c E, m_3 = E.  Code gives
    IM2 <- s*IM2 + c*IM3, IM3 unchanged.  Matches.
  * X3, frame {nhat_eta, rhat, e_xi}: need m_2 = X, m_3 = s M + c X.  Code gives
    IM2 unchanged, IM3 <- s*IM3 + c*IM2.  Matches.

No O(1) term, and the 1/s factors are benign (sin_cell >= 0.88 everywhere).

**The GEOMETRIC SOURCE TERMS are also correct.** `SrcTermsGnomonicEquiangleImpl`:
src1 reduces exactly to z_ov_rE*(p + B_r^2/2 + rho|v_t|^2/2) as its comment claims, and
src2/src3 carry the right Maxwell partners (s B^eta)^2 = b3^2 and
(s B^xi)^2 = (s b2 - c b3)^2, given B^xi = b2 - c b3/s and B^eta = b3/s.

**The CORNER EMF is not it either.** The cs branch of `mhd_corner_e.cpp` deliberately
drops GS07 and uses a plain four-face average ("Until the cell-centred EMF is given its
own gnomonic form"), and spherical polar -- which works at low beta -- takes the GS07
branch, so this looked like THE difference. Forcing cs onto the GS07 branch via a
temporary switch (since reverted): per-step vertex order at beta = 0.1 goes from -0.04 to
**+0.00**. Still inconsistent, and larger in absolute terms (1.68e-3 vs 6.54e-4, expected
since e_cc is in the wrong frame on cs). **Changing the corner-EMF algorithm entirely does
not restore consistency**, so the defect is upstream of it.

## HOW TO GATE ANY FIX

Run iprob=9 at nlim=1 and check the per-step VERTEX order at beta = 0.1 and 0.01. It must
be ~+2, not ~0. **Do not gate at the default beta = 2, where the defect is invisible**, and
do not gate on a many-step run, where accumulation muddies the order.

## WHERE TO LOOK NEXT

PLM reconstructs PRIMITIVES, and at the vertex its stencil reaches into the corner halo.
The halo VALUES are exact (verified: transform 1e-16, resample 2e-15), but a slope taken
between a cell in chart A and a halo cell transformed from chart B need not be a
meaningful slope -- the values are right at their own locations while the sequence is not
smooth in the destination's index space. That would hurt B most exactly when B dominates
the flux, which is the beta threshold observed. Start at the PLM call for the tangential
field components in `mhd_fluxes.cpp` and how its stencil meets the cube-vertex halo.

## 2026-09-04 RE-MEASUREMENT: the per-step operator is SECOND ORDER EVERYWHERE at beta=0.01

Same recipe (iprob=9, p0=1, b0c=14.142, nlim=1), RAW HLLD (`cs_lowbeta_fallback=0`), nx 16/32/64,
region bins as cs_test now defines them (fixed PHYSICAL distance to the vertex):

    region     L1(v)                          L1(p)                       L1(B)
    interior   1.09e-3  2.13e-4  6.30e-5      1.04e-1 2.60e-2 6.49e-3     8.97e-5 1.37e-5 1.95e-6
    seam       1.64e-3  3.21e-4  7.58e-5      1.05e-1 2.61e-2 6.51e-3     1.54e-4 2.61e-5 4.68e-6
    VERTEX     2.73e-3  5.66e-4  1.34e-4      1.06e-1 2.61e-2 6.49e-3     2.26e-4 3.94e-5 9.82e-6

Ratios 4-5 per refinement in every region and variable, vertex included. **The "order zero at
the vertex" table earlier in this file is SUPERSEDED** -- it was already retracted in
[[cs-mhd-instability-characterized]] (fixed-CELL bins that shrink under refinement), and these
are the numbers on the corrected bins. The seam's per-step error is ~1.5x the interior's and
the vertex's ~2.5x, both converging at the same rate: a constant, not a mechanism.

`<mhd>/cs_seam_econsist` (the halo-energy fix, b9efbbfa) moves these by ~0.1 % (seam L1(v)
1.6405e-3 -> 1.6418e-3, marginally WORSE). It is NOT a consistency cure; there is no seam
consistency defect left at this level. NOTE the fallback defaults to 0.5 on cs, so a sweep
that does not set it to 0 measures HLLE, not HLLD -- that also converges at 2nd order.

**Consequence**: the cubed-sphere low-beta instability is a NONLINEAR instability of a
CONSISTENT HLLD+CT scheme on the non-orthogonal grid (the odd-even / carbuncle class, where
consistency does not protect you), and the HLLE fallback below beta 0.5 is a legitimate
stabiliser of it -- the same status as spherical polar's polar averaging. The production
configuration (corrected rotation + fallback) is defensible on those terms; what remains
unproven is that nx=64 recovers through its dt trough as nx=32 did.
