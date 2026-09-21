# The cell-centred x1x2x3 CORNER buffers across a cubed-sphere panel seam

`tests_seam4/README.md` (commit `94c7165d`) gave the cell-centred along-seam resample to
the x2/x3 FACE buffers (slots 8-39) and to the x2x3 EDGE buffers (slots 40-47).
`tests_seam_fc/README.md` (commit `2268e633`) then measured `bvals_fc.cpp` and found that
it already covers slots 40-**55**, i.e. the x1x2x3 CORNERS as well, while `bvals_cc.cpp`
stopped at 48 -- so a cell-centred corner buffer crossing a panel seam was a PLAIN COPY.

This note measures that gap and closes it.  `x1` is RADIAL on the cubed sphere, so a
corner ghost is a radial ghost cell that is also tangentially ghost in x2 and x3; it
exists only when the mesh has more than one MeshBlock in x1, and it is read by multi-D
stencils near a radial block boundary at a seam and by prolongation.

## What changed

* `src/bvals/bvals_cc.cpp`, `MeshBoundaryValuesCC::PackAndSendCC`, the `do_cs` branch
  only: the slot test `n >= 40 && n < 48` becomes `n >= 40 && n < 56`, and the
  slot -> flanking-face decode is written for both layouts
  (edge `n = 40 + n1 + 2*(iy>0) + 4*(iz>0)`, corner
  `n = 48 + (ix>0) + 2*(iy>0) + 4*(iz>0)`, whose bit 0 is the RADIAL side and carries no
  tangential information).  Everything downstream -- the 3-point along-seam stencil, the
  monotonicity clamp with the `|u| > 1` extrapolation guard of `94c7165d`, the
  `seam_extent < 3` narrow-block guard, the clamp to the source's active range -- is
  the code the edges already run, unchanged and now simply reached by eight more slots.
  That is the whole source change: 2 hunks, ~30 lines, all inside `if (do_cs)`.
* `src/pgen/cs_test.cpp`, `CSTestSeamHaloScan`: a THIRD pass over the corner ghosts, with
  three bins (same-panel / SEAM / CUBE VERTEX).  Scan tooling only.

`src/bvals/buffs_cc.cpp` was NOT touched.  No non-cubed-sphere path is touched: the hunks
live inside `if (do_cs || do_pole)` and only act when `cs_seam != 0`, which no Cartesian,
spherical-polar or polar-boundary buffer ever sets.

## The measurement

The `iprob = 15` state is a spherical harmonic of the DIRECTION COSINES alone, so the
exact value in a ghost is known with no reference to the neighbour panel and does not
depend on radius -- which is why the same scan works verbatim in a radial ghost.
`time/nlim = 0`, so one boundary exchange is the only thing that has touched the ghosts.
`mesh/nx1=16 meshblock/nx1=8` gives the two radial blocks the corner slots need.

```bash
cmake -B build_cs_cc -D PROBLEM=cs_test -D CMAKE_BUILD_TYPE=Release
make -C build_cs_cc -j8
B=build_cs_cc/src/athena; I=inputs/tests/cubed_sphere_raddiff.athinput
for n in 32 64; do for mb in $n $((n/2)) $((n/4)); do
  $B -i $I -d out mesh/nx1=16 meshblock/nx1=8 mesh/nx2=$n mesh/nx3=$n \
     meshblock/nx2=$mb meshblock/nx3=$mb time/nlim=0 problem/seam_halo_scan=1 \
     output1/dt=1e9 | grep -E "GLOBAL MAX|max ="
  rm -rf out
done; done
```

### max |ghost - exact| in the CORNER bin (slots 48-55), before -> after

| bin | n | 1 x 1 | 2 x 2 | 4 x 4 |
| --- | --- | --- | --- | --- |
| corner, same panel | 32 | 0 (exact) | 0 (exact) | 0 (exact) |
| **corner, SEAM** | 32 | none exist | **1.407e-4 -> 7.986e-8** | **7.414e-4 -> 5.077e-6** |
| corner, CUBE VERTEX | 32 | 2.911e-5 (unchanged) | 2.911e-5 | 2.911e-5 |
| corner, same panel | 64 | 0 (exact) | 0 (exact) | 0 (exact) |
| **corner, SEAM** | 64 | none exist | **3.498e-5 -> 4.284e-9** | **3.538e-4 -> 5.639e-7** |
| corner, CUBE VERTEX | 64 | 3.662e-6 (unchanged) | 3.662e-6 | 3.662e-6 |

After the fix the corner-SEAM bin equals the x2x3 EDGE-seam bin to every printed digit in
every layout (7.986e-8 / 5.077e-6 at n = 32, 4.284e-9 / 5.639e-7 at n = 64), which is the
claim: corners now carry exactly the same seam halo accuracy as the edges next to them,
and the residual at 4 x 4 is the same `O(dx^3)` interior-block-end extrapolation
`tests_seam4/README.md` lists under "What is NOT fixed".

Controls that say the scan measures what it claims:

* the same-panel corner bin is EXACT (0) before and after -- a plain index copy on one
  chart, which is what it must be;
* every other bin (face same-panel / face SEAM / edge same-panel / edge seam / edge cube
  vertex) and the `SEAM-HALO GLOBAL MAX` are byte-for-byte the same before and after;
* at 1 x 1 blocks per panel the corner-SEAM bin is EMPTY (n = 0): every cross-panel corner
  is a CUBE VERTEX there, so `cs_seam` stays 0 and the branch cannot fire.

### the negative control: the resample is what does it

`CS_NORESAMP_CC=1` (the existing env hook at `bvals_cc.cpp:65`) forces `cs_seam = 0` on
every buffer.  With the fix built in, n = 32, 4 x 4 blocks:

| bin | resample ON | resample OFF |
| --- | --- | --- |
| edge, seam | 5.0766e-06 | 7.4137e-04 |
| **corner, SEAM** | **5.0766e-06** | **7.4137e-04** |
| corner, same panel | 0 | 0 |

The OFF column reproduces the pre-fix corner number (7.414e-4) exactly, so the new branch
and nothing else is what moved it, and the corner path now tracks the edge path bit for
bit in both directions.

### the CUBE VERTEX interaction

Two independent reasons the vertex fill is untouched.  (i) `PackAndSendCC` skips the
exchange altogether for a cube-vertex slot (`IsCubeVertexCorner`, top of the kernel), and
`FillPanelCornersCC` runs later and owns exactly those cells.  (ii) even if it were
reached, a cube-vertex corner has BOTH tangential flanks on a seam, so the `s2 && !s3` /
`s3 && !s2` decode leaves `cs_seam = 0`.  Measured: the corner CUBE VERTEX bin is
2.9113e-05 (n = 32) and 3.6619e-06 (n = 64) before AND after, bit for bit, in all three
layouts.  The two paths do not both write the same ghost cell.

## Gates

Reference binary = this tree with the two source files at `544cb36b`.  `-j8`, Release,
no MPI, one rank.

### (b) + (d) the cubed-sphere blasts: bitwise, and conservation

`inputs/tests/cubed_sphere_blast.athinput` (hydro, vertex-centred, 278 cycles to
t = 0.30) and `cubed_sphere_mhd_blast.athinput`, as shipped (1 radial block, 1 x 1 blocks
per panel) and with `mesh/nx1=16 meshblock/nx1=8` so the corner slots exist, at 1 x 1,
2 x 2 and 4 x 4 blocks per panel.  The comparison is the full-precision binned profile
dump (`cs_blast-blast.dat`) plus the code's own conserved sums.

| input | radial blocks | blocks/panel | dumps | mass, energy |
| --- | --- | --- | --- | --- |
| hydro blast | 1 | 1 x 1 | **BITWISE IDENTICAL** | identical |
| hydro blast | 2 | 1 x 1 | **BITWISE IDENTICAL** | identical |
| hydro blast | 2 | 2 x 2 | **BITWISE IDENTICAL** | identical |
| hydro blast | 2 | 4 x 4 | **BITWISE IDENTICAL** | identical |
| MHD blast | 1 | 1 x 1 | **BITWISE IDENTICAL** | identical |
| MHD blast | 2 | 1 x 1 | **BITWISE IDENTICAL** | identical |
| MHD blast | 2 | 2 x 2 | **BITWISE IDENTICAL** | identical |
| MHD blast | 2 | 4 x 4 | **BITWISE IDENTICAL** | identical |

Conservation, hydro blast, 2 radial blocks x 4 x 4 (`d0 = 1`, so the exact final mass is
the volume): volume 2.93215314335059070e+01, mass 2.93215314335055410e+01, i.e.
`|mass - volume|/volume = 1.2e-16`, identical to the last bit before and after.

**This is the honest result of the dynamic gates: nothing changed.**  Neither blast reads
a cell-centred x1x2x3 corner ghost -- the hydro update is a dimensionally split PLM sweep
and the CT corner EMF reads the x2x3 diagonal (slots 40-47, already resampled since
`94c7165d`), not the radial-tangential one.  The gates therefore show only that the change
is CONFINED: it cannot perturb anything that was already right.  The benefit measured here
is the static one above, for the consumers that do read these cells (multi-D stencils at
a radial block boundary on a seam, and prolongation).

### (c) the diffusion seam test at 1 x 1 / 2 x 2 / 4 x 4

`tests_adi/measure_cs.py` on `cs_test` `iprob = 15`, `n = 64` per panel edge,
`cfl = 0.0075`, `nlim = 80`, `mesh/nx1=16 meshblock/nx1=8`:

| solver | blocks/panel | amp | L1 | Linf | edge | interior |
| --- | --- | --- | --- | --- | --- | --- |
| explicit, before | 1 x 1 / 2 x 2 / 4 x 4 | 0.74768 | 0.2601 | 0.3406 | 0.3385 | 0.3406 |
| explicit, **after** | 1 x 1 / 2 x 2 / 4 x 4 | 0.74768 | 0.2601 | 0.3406 | 0.3385 | 0.3406 |

Every column is identical before and after and identical across the three layouts -- the
block-count independence `94c7165d` established is preserved, and the corner buffers add
nothing to this operator (it is purely angular and never reaches across a radial block
boundary).  `amp` is far from 1 only because this arm splits the radial extent over two
MeshBlocks, which is not what the `tests_adi` ladder was calibrated on; the comparison
that matters is before vs after, and it is exact.

The `adi` and `sts` arms CANNOT be run here: `conduction.cpp:530` fatals with
`rad_ang_solver = adi needs the whole x1 extent in one MeshBlock`, and with one radial
block there are no corner slots at all.  So those solvers are unreachable by this change
by construction.

### (c) the cubed-sphere MHD convergence test

`inputs/tests/cubed_sphere_mhd_conv.athinput` (`iprob = 9`, rigid rotation carrying a
uniform force-free field), `mesh/nx1=16 meshblock/nx1=8`, 4 x 4 blocks per panel
(`conv_{32,64}_{ref,after}.log`).  Note the printed `nx2=` is the MESHBLOCK size.

| n per panel | L1(B) | Linf(B) | L1(v) | L1(p) |
| --- | --- | --- | --- | --- |
| 32 | 2.222289e-04 | 1.147781e-03 | 4.642670e-04 | 1.482576e-04 |
| 64 | 1.084289e-04 | 5.563696e-04 | 2.248281e-04 | 7.370611e-05 |

`L1(B)` ratio 2.0495 -- second order, as the input's header requires.  **Before and after
are identical in every `###` diagnostic line at both resolutions**, region tables (panel
interior / panel seam / cube vertex), the L1(B) seam-distance and radial-distance
budgets, the best-fit uniform B and the flow rate included:
`diff <(grep '###' conv_N_ref.log) <(grep '###' conv_N_after.log)` is empty for
`N = 32` and `N = 64`.

### (e) the regression test

```bash
cd tst && python3 run_test_suite.py --test test_suite/rad/test_rad_cs_raddiff_cpu.py
```

**1 passed in 94.30 s.**

## Style

`src/bvals/bvals_cc.cpp` is cpplint-identical to `544cb36b` (9 reports, unchanged).
`src/pgen/cs_test.cpp` goes from 115 to 116 reports; the one addition is
`runtime/int` on `static_cast<long long>(knum[c])` in the new printf, which is the same
cast the 17 pre-existing `static_cast<long long>(cnum[c])` / `(cnum[c])` calls in this
same function already use -- deviating from them would be the odd thing.  No tabs, no
line over 90 columns, no `}}`, no trailing whitespace in any added line.

## What is NOT claimed

* No dynamic gate in this tree resolves a difference.  The corner ghosts are correct now
  where they were O(dx) wrong before, but nothing currently gated reads them.
* The 4 x 4 residual (5.077e-6 at n = 32 against the 4.103e-6 face background) is the
  interior-block-end extrapolation, inherited unchanged from the edge path.
* `bvals_fc.cpp` is untouched; it already covered slots 48-55.
* MPI was not built or tested here; the change is in the pack kernel and is
  rank-layout independent, but that is an argument, not a measurement.
