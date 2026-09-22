# The cubed-sphere CUBE VERTEX: what fills it, and what a real fill buys

Session of 2026-09-22.  Everything here is serial CPU, Release, `-D PROBLEM=cs_test`.
The two binaries compared are built from **pristine `HEAD` (0e11d533)** and from HEAD plus
the patch described below -- *not* from the working tree, because another session was
editing `src/hydro`, `src/mhd` and `src/reconstruct` concurrently and those edits change
the answers (see "A working-tree hazard" at the end).

```bash
./scan_vertex.sh   ./athena_new          # gate (b): the STATIC vertex-ghost error
./regions_vertex.sh ./athena_ref base     mesh/cs_vertex_fill_cc=false   # gate (c)
./regions_vertex.sh ./athena_new ccfill   mesh/cs_vertex_fill_cc=true
./regions_vertex.sh ./athena_ref fcextrap mesh/cs_vertex_fill=false mesh/cs_vertex_fill_cc=false
```

## What ALREADY existed at a cube vertex before this session

The brief assumed the vertex ghosts are filled by extrapolation and that a sampling fill
was never built.  That is **half true**, and the half matters:

| data | today, at a cube vertex | switch |
| --- | --- | --- |
| face-centred `b.x1f/x2f/x3f` | **SAMPLED** from the two flanking panels: the x2/x3 face and x1x2/x3x1 edge buffers are widened by `ng` along the seam, ownership is decided receiver-side on the panel diagonal, and `FillPanelCornersFC` stands down where that lands (`buffs_fc.cpp`, `bvals_fc.cpp`) | `<mesh>/cs_vertex_fill`, default **true** |
| cell-centred (`u0`, `w0`, radiation, ...) | **EXTRAPOLATED**: two one-sided quadratics from the flanking face halos, reach-weighted mean, sign-preserving fallback (`FillPanelCornersCC`) | none -- no sampling path existed |
| corner EMF `x1e` on the 3-valent radial edge | **THREE-WAY AVERAGE** of own + the two flanking seam flux buffers (`SavePanelCornerEMF` / `AveragePanelCornerEMF` in `flux_correct_fc.cpp`), *at the same level only* | none |

So (A) was already done for FC and missing for CC, and (B) was already done -- at one
level.  Both statements are now measured, not read off the comments.

## What was built

`<mesh>/cs_vertex_fill_cc`, default **false** (bitwise), the cell-centred twin of the
face-centred fill: `WidenForCubeVertex` in `buffs_cc.cpp` widens the same slots by the
same rule on send and receive; `RecvAndUnpackCC` gates doubly-ghost writes of slots 8-39
by the same receiver-side diagonal ownership test (cell centres, no staggering);
`FillPanelCornersCC` stops extrapolating exactly where the sampled data lands, including
the radial-ghost layers, which need both flanking x1x2/x3x1 edge buffers at the same
level.  No new slot, no new tag, nothing read from the sender -- MPI-safe by the same
argument as the face-centred fill, **but see "MPI" below: it could not be run.**

## Gate (a): defaults are bitwise

`cs_test` cs hydro blast (iprob 12), cs MHD blast, and the `cs_regions_loop` field loop,
5 cycles each, binary dumps compared variable by variable
(`cmp_bin.py`, which reads the payload: the raw files differ only by the extra echoed
input line):

```
cs_blast.hydro_w      max|diff| = 0.0
cs_mhd_blast.mhd_w_bcc max|diff| = 0.0
cs_reg_loop.mhd_w_bcc  max|diff| = 0.0
```

## Gate (b): the STATIC vertex-ghost error

`cs_test` iprob 15 (a spherical harmonic of the direction cosines, exact value in a ghost
known from its own chart-continued angles), `time/nlim = 0`, so one boundary exchange and
nothing else; `mesh/nx1=16 meshblock/nx1=8` so the radial-ghost ("corner") vertex bin
exists.  Max `|ghost - exact|` in the CUBE-VERTEX bins:

**Cell-centred (the new switch).  1x1 and 2x2 MeshBlocks per panel give IDENTICAL numbers.**

| n | blocks/panel | extrapolate (today) | sample (`cs_vertex_fill_cc=true`) | gain |
| --- | --- | --- | --- | --- |
| 32 | 1x1 and 2x2 | 2.9113e-05 | 2.5785e-06 | **11.3x** |
| 64 | 1x1 and 2x2 | 3.6619e-06 | 3.4288e-07 | **10.7x** |

Both converge (extrapolate 2.99 orders, sample 2.91), so the sampled fill is a constant
11x better, not a different order.  That beats the 3.4-5.5x the 2026-08-31 prototype
recorded for the face-centred side.

**Face-centred (the switch that already exists), same runs, `seam_halo_scan_fc`:**

| n | component | extrapolate | sample (default) |
| --- | --- | --- | --- |
| 32 | x1f / x2f / x3f | 2.7131e-05 / 2.4758e-04 / 2.8925e-04 | 6.8603e-05 / 6.6605e-05 / 3.8931e-05 |
| 64 | x1f / x2f / x3f | 8.1596e-06 / 6.0015e-05 / 4.8421e-05 | 1.6476e-05 / 6.5975e-05 / 3.5063e-05 |

**This is a finding, and not a happy one.**  At n = 32 the sampled fill is 3.7x / 7.4x
better on the two TANGENTIAL components (the recorded 3.4-5.5x), but 2.5x WORSE on the
radial one -- and it **does not converge**: x2f is 6.66e-05 at n = 32 and 6.60e-05 at
n = 64 (order 0.01), while the extrapolation it replaced converges at order 2.05
(2.48e-04 -> 6.00e-05).  By n = 64 the sampled fill is better on x3f only, and worse on
x1f and x2f.  A stagnant O(1e-4) face ghost at the vertex is a candidate explanation for
the field loop's vertex L1 order of 1.57 that nothing has moved; it is NOT explained by
the cell-centred side, which converges cleanly.

## Gate (c): the per-region dynamic gate

`tests_cs_regions` configuration exactly (plm, `cs_wellbalanced_src` off, nghost 3, one
block per panel, `conv_nband` 2/4/8 at n = 16/32/64, fixed physical time), three arms:
`base` = today's defaults, `ccfill` = `cs_vertex_fill_cc=true`, `fcextrap` =
`cs_vertex_fill=false`.

**`loop` (iprob 11), L1 of the evolved FACE field, and the measured order:**

| arm | INTERIOR n=16/32/64 | VERTEX n=16/32/64 | vertex order | vertex/interior at n=64 |
| --- | --- | --- | --- | --- |
| base (today) | 5.6015e-05 / 1.3758e-05 / 3.4411e-06 | 9.8103e-05 / 3.3125e-05 / 1.0542e-05 | 1.57 / 1.65 | 3.06 |
| ccfill | 5.6010e-05 / 1.3757e-05 / 3.4409e-06 | 1.0243e-04 / 3.3758e-05 / 1.0582e-05 | 1.60 / 1.67 | 3.08 |
| fcextrap | 5.6019e-05 / 1.3757e-05 / 3.4410e-06 | 9.7941e-05 / 3.3219e-05 / 1.0568e-05 | 1.56 / 1.65 | 3.07 |

The `base` column reproduces the 76359975 baseline to every digit (1.57/1.65, ratio 3.06).

**`strat` (iprob 13), L1 of the TANGENTIAL velocity of an atmosphere at rest:**

| arm | VERTEX/INTERIOR at n = 16 / 32 / 64 |
| --- | --- |
| base | 2.52 / 2.72 / 2.875 |
| ccfill | 2.53 / 2.72 / 2.874 |
| fcextrap | 2.56 / 2.72 / 2.876 |

**The verdict is negative and it is the important result of this session.**  Sampling the
cube-vertex ghosts instead of extrapolating them -- cell-centred, face-centred, or both --
moves the DYNAMIC vertex error by under 0.5 % at every resolution, leaves the vertex
convergence order at ~1.6 and the vertex/interior ratio at ~3.07, and the strat tangential
ratio at 2.87.  Turning the existing face-centred fill OFF changes nothing either.  So
**the cube-vertex ghost halo is not what limits the vertex error**: an 11x better
cell-centred ghost and a 3.7-7.4x better tangential face ghost buy nothing dynamically.
What remains on the list is the other half of the brief's background -- the flux-divergence
versus geometric-source cancellation, 1:600 in a panel interior and 1:260 at a vertex, and
the CT geometry of the 3-valent edge itself.

`max|div B|` is round-off in every arm and every resolution; the MHD blast centred on a
cube vertex runs 200 cycles cleanly with the fill on.  Cost, measured serially (loop,
n = 32, one run each): 14.65 s -> 14.72 s, i.e. **+0.5 %** for the cell-centred fill
(4.867e+05 -> 4.844e+05 zone-cycles/s); dt is unchanged.

## (B) The corner EMF at the 3-valent vertex

Measured with the existing `CS CUBE-VERTEX EMF SPREAD` diagnostic (`CSTestGhostCheck`,
reached with `problem/iprob=8 problem/conv_errors=0`), which keys the radial edge by
GEOMETRY and reports the spread among the panels that share it:

| mesh | spread | max&#124;e1&#124; | verdict |
| --- | --- | --- | --- |
| uniform (`cubed_sphere_mhd`, 20 cycles) | 1.36e-20 | 1.20e-04 | single-valued to **round-off** (1e-16 relative) -- the three-way average is correct and consistent, in every arm |
| SMR (`cubed_sphere_mhd_smr`, 20 cycles) | 3.22e-04 | 7.40e-05 | **multi-valued by 4x its own magnitude** |

So there is nothing to fix at one level -- `AveragePanelCornerEMF` already makes all three
panels average the same three numbers, and the shared seam faces stay equal to 8.9e-16.
The **level boundary is the real gap**, and it is the one the code documents and skips
(`AveragePanelCornerEMF` returns when a flanking neighbour is at a different level).
Closing it needs the finer neighbour's length-weighted EMF sum rescaled by this block's
own edge length and the coarser neighbour's value taken whole; that was NOT attempted
here.  `max|div B|` is unchanged by every switch in every arm (round-off on the uniform
mesh; the SMR mesh's 1e-5 is pre-existing and identical in all three arms).

## (d) Block decomposition and MPI

1x1 vs 2x2 MeshBlocks per panel give bit-identical vertex-ghost numbers in gate (b), for
both fills, at both resolutions.  **MPI is UNTESTED**: this machine has no MPI package
visible to CMake (`MPI package required but could not be found`), so no MPI binary could
be built.  The cell-centred fill inherits the face-centred one's MPI argument -- send and
receive widen by the identical rule so message sizes match, no new slot or tag is
created, and ownership is computed from the receiver's own geometry -- but that argument
is not a test.

## A working-tree hazard, unrelated to this work

With the concurrent session's in-flight `src/hydro` + `src/mhd` + `src/reconstruct` edits
compiled in, **every hydro-only cubed-sphere input segfaults in
`ProblemGenerator::UserProblem` at `time/nlim = 0`** (`cubed_sphere_blast`,
`cubed_sphere_rigidrot`; exit 139).  The same inputs run cleanly on pristine HEAD.  Worth
telling that session before it commits.
