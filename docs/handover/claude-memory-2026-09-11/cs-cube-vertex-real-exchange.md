---
name: cs-cube-vertex-real-exchange
description: PROTOTYPED AND MEASURED, deliberately NOT built -- filling the cubed-sphere cube-vertex corner by SAMPLING the panel that contains each ghost cell (3.4-5.5x on two halo categories) instead of extrapolating; and the design correction that it needs NO new exchange, only a wider destination range on the existing face buffer
metadata:
  type: project
---

Prototyped 2026-08-31 behind `CS_VXFILL`, measured, then REVERTED. Read this before
proposing either "extrapolate better" or "build a new cube-vertex exchange".

## THE GEOMETRY, which I first got wrong twice

At a cube vertex only THREE panels meet. The corner ghost region of panel P, in P's own
extended gnomonic chart, is **covered by the TWO FLANKING panels A and B**, split by a
diagonal -- NOT by an inaccessible fourth panel, and not by some third panel needing a new
neighbour relation. Concretely for P = +z, A = +x, B = +y: beyond P's corner both x > z and
y > z, so the point is on +x where x is largest and on +y where y is largest.

**And the cells needed are A's / B's own ng x ng ACTIVE corner block** -- within ng cells
of the shared edge and within ng of their other edge. **The existing x2/x3 FACE buffer
already packs exactly those cells** (it takes a ng-deep strip over the FULL along-seam
range, whose end IS the corner block). So the data already crosses the wire.

## THEREFORE: NO NEW EXCHANGE IS NEEDED

The correct implementation is to **extend the face buffer's DESTINATION range** to cover
the corner block, and give those destination cells a 2-D map (target -> Cartesian ->
`FindPanel` -> `CartToPanel` -> fractional index -> biquadratic sample -> 
`TransformFieldToDstNormals`). It rides existing communication, so it is MPI-safe by
construction, and it is a fraction of the work a bespoke exchange would be. **Do not start
by designing a new buffer/slot/tag.**

## MEASURED PAYOFF (prototype, same-rank direct sampling)

|  | tangent | r x EDGE SEAM | r x CORNER SEAM |
|---|---|---|---|
| extrapolate, refined | 1.4860e-03 | 4.6508e-04 | 2.0911e-03 |
| extrapolate, control | 1.4860e-03 | 3.9044e-04 | 1.3047e-03 |
| **sample, refined** | **4.3188e-04** | 4.6508e-04 | **3.8037e-04** |
| **sample, control** | **4.3188e-04** | 3.9044e-04 | **3.8545e-04** |

3.4-5.5x, nothing regressed, and the cube-vertex corner lands at ORDINARY HALO QUALITY
(the same-panel r x EDGE halo is 2.9e-04). Domain gates barely move, which is why it was
not built: refined resistive Linf 5.6995e-04 -> 5.4140e-04 (5.0%), L1 1.3064e-05 ->
1.3012e-05, heating unchanged; ideal mhd_smr L1(B) 2.924547e-04 -> 2.924534e-04 with
Linf(B) marginally WORSE.

## WHY IT WAS NOT COMMITTED

The prototype reads the source block straight out of the MeshBlockPack, so under MPI a
corner whose source panel is off-rank finds nothing; falling back to extrapolation there
would make results RANK-DEPENDENT, which fails this project's standing bar. Building the
proper version means the destination-range change above -- real work in code where this
session alone produced two subtle regressions of my own making. **5% on one domain gate
does not justify that risk while no test demands it.**

**WHAT WOULD CHANGE THE ANSWER:** a problem with real STRUCTURE at a cube vertex. The
iprob=11 field is smooth, so an 8x-worse corner halo has little to bite on; a strong
gradient parked on a vertex would make the 3.4-5.5x matter far more than the 5%.

## A TRAP THE PROTOTYPE HIT, and it is the recurring one

The first version used the FINE index set (`nx2`, `js`) even when `FillPanelCornersFC`
runs on `coarse_b0`. It silently mis-sampled and knocked the prolongated fine EDGE halo
back 5x to its pre-180a9b3e value. **The coarse array is where these bugs hide** -- three
times now (180a9b3e, the narrow-block guard, this). Anything that indexes a grid inside
that routine must branch on `coarse`.


## THE BOOLEAN SUM IS MEASURED WORSE -- do not re-propose it (2026-08-31)

Tried and REVERTED. The corner is filled by two 1-D quadratic extrapolations, one from
each flanking halo, combined by a reach-weighted mean. The textbook completion is the
**Boolean sum** `ek + ej - E_j E_k(interior)`, exact for any `g(j) + h(k)`. It is WORSE:

```
                          evolved Linf   tangent ghost   r x CORNER SEAM
reach-weighted mean (HEAD)  5.6994e-04     1.6455e-03      2.1705e-03
Boolean sum, quadratic      5.4373e-04     3.3100e-03      5.3740e-03
Boolean sum, linear cross   4.7471e-04     9.3265e-03      1.0613e-02
```

**Why:** the cross term is a DOUBLE extrapolation, amplification 7x7 = 49 at the first
ghost cell -- the same amplification that already killed the cubic variant. The two 1-D
extrapolations read DIFFERENT halos, so their errors are largely independent and
averaging them CANCELS noise; the Boolean sum instead subtracts a hugely amplified
common term. Dropping the cross term to linear (amplification 9) makes it worse still.

**And note which gate can even SEE this.** `cubed_sphere_mhd_smr` is unchanged across all
three modes (2.925412e-04 identically in modes 0 and 2) -- **not** corroboration, but a
gate that is INSENSITIVE BY CONSTRUCTION: ideal MHD never reads the corner block, only
the resistive curl does. Do not quote it as a second opinion on any corner-fill change;
use `cubed_sphere_resist_smr` and the GHOST SCAN lines.

**And distrust the domain gate here.** The evolved-field Linf IMPROVES monotonically as
the corner halo DEGRADES monotonically -- a backwards anti-correlation, so the 17% "gain"
at linear-cross is cancellation against an existing error, not a better corner. A
principled fill must improve BOTH. Corner ghosts wrong by 1e-2 of |B| are disqualifying
for shock work regardless of what the smooth gate says.

**Conclusion: the local-data estimator family is EXHAUSTED** (linear, quadratic, cubic,
reach-weighted mean, Boolean sum x2). The corner ghost of P maps into panel A's ACTIVE
corner block, but P only ever STORES its x2 halo over its own active eta range, so the
data exists on A and P never receives it. Extrapolation is provably the best local
option. **Only the wire change closes this.**

## GPU COST OF THE WIRE CHANGE, for the record

The corner fill kernel is `nmb x 4 x nc1 x ng^2` (~320 cells/block against 4096 active),
so ANY arithmetic-only estimator is free on GPU -- that is why the Boolean sum was worth
testing and why a 6x6 least-squares solve per cell was rejected on occupancy grounds
without testing. The real exchange is the one that costs: widening every x2/x3 face
buffer by ng at each end grows the packed strip by ~2*ng/nx (about **25% on a 16-cell
block, on all six panels, not just the eight vertices**), and the extra destination cells
need atan/tan plus a biquadratic sample in the pack kernel, which is already part of the
~9.6% boundary cost. No atomics, no new sync, MPI-safe by construction.
