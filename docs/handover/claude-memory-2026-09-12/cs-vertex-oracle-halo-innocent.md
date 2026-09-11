---
name: cs-vertex-oracle-halo-innocent
description: ORACLE (2026-09-06) -- EXACT panel ghosts (edges, or edges+corners) change the cube-vertex one-step force residual by <1 %; the vertex's 1-in-500 is the CELL BALANCE, not the halo. Any ghost-stencil idea (incl. the through-vertex seam continuation along the third panel's diagonal) cannot touch it
metadata:
  type: project
---

cs_test iprob=8 (uniform Cartesian field at rest), one step, GS07 + WB source, raw HLLD,
`problem/exact_panel_ghosts` = 0 (code) / 2 (exact edge ghosts, code corners) / 1 (exact all).
|F| in units of (B^2/2)/r, per region (vertex = the 2x2 corner cells; scratchpad sp_src/vtx.py):

    nx/panel  ghosts       interior   edge rows  vertex 2x2  vertex max
    32        code         2.086e-3   3.630e-3   7.060e-3    3.222e-2
    32        exact edges  2.086e-3   3.644e-3   7.112e-3    3.204e-2
    32        exact all    2.086e-3   3.644e-3   7.111e-3    3.194e-2
    64        code         5.695e-4   1.227e-3   3.246e-3    1.465e-2
    64        exact all    5.695e-4   1.228e-3   3.228e-3    1.467e-2

So the halo is already "perfect" as far as the momentum balance can tell: the residual at
the edge rows (2x interior) and the vertex (6x, worst cell 25x) is the cell's own flux-
divergence-vs-source cancellation on a sheared basis. The user's geometric idea (a seam
continued through the vertex lands on the third panel's chart diagonal -> a centred
stencil) is exact geometry but has NO headroom here; DECIDED not to build. It could only
matter for an operator whose vertex error is halo-limited, and the resistive curl was
already shown 2nd order ([[cs-resistive-seam-order]]). What would move the vertex number is
the cell balance itself: a Cartesian-momentum update in the singular cells (not built).
See [[sp-geometric-source-residual]] for the sp twin.
