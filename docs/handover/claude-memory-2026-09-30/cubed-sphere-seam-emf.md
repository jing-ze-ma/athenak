---
name: cubed-sphere-seam-emf
description: Seam EMF consistency for cubed-sphere MHD -- FULLY FIXED, 42323a66 + 7b421818 + 2646f826; the seam field is exactly single-valued for 1, 4 and 16 MeshBlocks per panel, cube vertices and four-block seam edges included. Also records that div B is the WRONG gate for this
metadata:
  type: project
---

**CLOSED 2026-08-28 (sixteenth session), committed 42323a66 + 7b421818 + 2646f826.**
This closes [[cubed-sphere-for-hot-jupiter]] item 1 completely. **START HERE for cs+MHD.**

## The single most reusable lesson: div B is the WRONG gate

CT conserves `sum(area*B)` over a cell's six faces **identically, for ANY edge EMF
whatever** -- each edge enters two of that cell's faces with opposite sign, so it
telescopes regardless of value. **No EMF error can break div B, and measuring div B says
nothing about seam consistency.**

A div B growth to 6e-2 that looked like a seam defect was actually **the test's own radial
BC overwriting an ACTIVE face**: `CSTestRadialBC` wrote `x1f(ie+1)`, which belongs to the
last active cell. Generalisable: when a BC sets a face-centred NORMAL component, set the
ghost faces only.

## What the problem was

A seam face is ONE physical face stored as an ACTIVE face of BOTH panels
(`InitRecvIndices` excludes it from the halo), so each panel evolves its own copy. They
stay equal only if every bounding edge carries the same line integral on both.
`SendE`/`RecvE` (`flux_correct_fc.cpp`) had **no panel handling at all**.

**The EMF needs NO basis transform.** CT consumes `dxedge*E`, and a seam edge is the same
physical curve in both charts, so it is a plain scalar. Three things do change: the
along-seam INDEX may reverse; the along-seam EMF then flips SIGN; and across a swap seam
the value goes into the neighbour's OTHER slot (x2e <-> x3e). x1e is untouched.

## Measured, all three fixes in (nx2=nx3=32, iprob=8, `user` radial BCs)

max over all 3072 shared seam faces of |b_P - b_Q|, and the x1e spread over block-corner
radial edges, at 400 cycles:

| | 1 blk/panel | 4 blk/panel | 16 blk/panel |
|---|---|---|---|
| seam face, away from vertex | 1.3e-15 | 9.4e-16 | 7.8e-16 |
| seam face, near a vertex | 6.9e-16 | 7.8e-16 | 1.1e-15 |
| x1e, cube vertex | 2.7e-20 | 1.4e-20 | 1.4e-20 |
| x1e, four-block seam | -- | 1.4e-20 | 2.7e-20 |
| x1e, panel interior | -- | 0 | 0 |

Before the fixes: seam face 2.1e-3 (blind) / 1.9e-3 (vertex, 42323a66 only) /
6.7e-4 and 7.6e-4 (four-block seam, after 7b421818). Spurious KE stops growing secularly.

## The three fixes, and why each generic mechanism failed

1. **42323a66 -- seam FACE buffers.** Panel-blind index map made the copies worse. Added
   the along-seam reversal, the sign flip and the x2e<->x3e slot swap.
2. **7b421818 -- the cube vertex.** Only THREE panels meet there; the edge machinery pairs
   a corner with ONE diagonal, which cannot express a three-way relation and on this grid
   points at the **WRONG CORNER of the right panel** (panel 0's (-x2,-x3) edge at
   (0.577,-0.577,-0.577) is paired with panel 5's (+x2,+x3) edge at (0.577,+0.577,-0.577)).
   Fix: `SavePanelCornerEMF` + `AveragePanelCornerEMF` -- each face buffer spans its whole
   seam INCLUDING both end edges, so a block's two seam buffers already hold the other two
   panels' x1e on that edge. own + those two, over 3. Stash own BEFORE `SumBoundaryFluxes`.
3. **2646f826 -- the four-block seam edge (multi-block).** **The cause was fix 2's bogus
   pairing landing somewhere else.** 7b421818 concluded the wrong-corner buffer "contributes
   only noise" -- true with one block per panel, FALSE with several: the pairing still
   carries a destination slot, and that slot is a **legitimate four-block seam edge of the
   RECEIVING block**, so the vertex send COLLIDES with the correct sender and one value is
   lost. Fix: skip the x2x3 edge exchange entirely at a cube vertex, send and receive
   together (`IsCubeVertexCorner`, applied in the pack kernel, the MPI send loop,
   `SumBoundaryFluxes` -- which must also stop COUNTING the slot -- and `InitFluxRecv`).
   The vertex edge is overwritten by `AveragePanelCornerEMF` anyway, so skipping is free.

**7b421818's "dropping the diagonal is measurably WORSE (4.2e-4 -> 1.4e-3)" is now
obsolete, not wrong:** it was measured before `AveragePanelCornerEMF` existed.

## The prerequisite that unblocked everything: FillPanelCornersCC

Without it the gate reached 1e39 on **cycle 2**. The cell-centred twin of the FC corner
fill. Hydro never needed it, but **MHD does**: a panel-corner EMF comes from fluxes whose
reconstruction reads those cells, so O(1) garbage there becomes an O(1) EMF -- 1.3% of the
magnetic energy lost on cycle 1.

## The gates to reuse (all in `cs_test.cpp`, under `iprob = 8`)

* `CS SEAM FACE SINGLE-VALUED` -- finds shared faces **by GEOMETRY**, keying every
  panel-boundary face on its physical centre, so it cannot inherit a bug from the index
  map under test.
* `CS CUBE-VERTEX EMF SPREAD` -- same for the corner radial edges.
* `CS BLOCK-CORNER x1e SPREAD` (2646f826) -- **the one that localised the last bug.**
  Buckets EVERY block-corner radial edge geometrically into cube vertex / four-block seam
  / panel interior. The aggregate number could not tell "the seam is broken" from "corners
  are broken"; splitting it showed the panel interior was already exact, which ruled out
  the generic edge machinery in one measurement.

Bucket "near a corner" GEOMETRICALLY (a cube vertex has |x|=|y|=|z|), never by index --
with several blocks per panel a block's seam ends are mostly ordinary interior edges.
**Run with `user` radial BCs, never `reflect`.** Run with `problem/iprob=8`; the gates do
not print for other iprob.

Two more audits worth repeating if this ever regresses: dump the x2x3 neighbour table
(gid/panel/dest) and check each pairing **geometrically** (24 of 96 mismatched, all cube
vertices); and dump `nflx` for n=40..46 (uniformly 4, which ruled out the count).

## Cubed sphere + MPI: also FIXED this session, see [[cubed-sphere-mpi-hang]] (1c2e29d6)

Harness: session a76ff032 scratchpad -- `runs/`, plus `bcs`/`bnew` builds in session
28cde4ec. MPI build recipe in [[cubed-sphere-mpi-hang]].
