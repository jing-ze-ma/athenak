# The cubed-sphere seam halo of the FACE-CENTRED field, with more than one MeshBlock per panel

`tests_seam4/README.md` fixed two defects of the CELL-CENTRED along-seam resample
(`src/bvals/bvals_cc.cpp`, commit `94c7165d`) that switch on as soon as a cubed-sphere
panel is split into more than one MeshBlock, and closed with

> `bvals_fc.cpp` (MHD) has the SAME two defects in its own resample [...] It is NOT
> touched here.

**That statement is wrong.  It was written from a reading of the code, not a
measurement.  `bvals_fc.cpp` has NEITHER defect, and it also already covers the
CORNER slots (48-55) that `bvals_cc.cpp` still leaves as a plain copy.**  This note is
the measurement that says so, the negative controls that prove the guards are
load-bearing rather than vacuous, and the multi-block gates.  **No line of
`src/bvals/bvals_fc.cpp` is changed.**

## What was added

* `src/pgen/cs_test.cpp` -- `CSTestSeamHaloScanFC`, the face-centred twin of
  `CSTestSeamHaloScan`, enrolled for `iprob = 11` when
  `<problem>/seam_halo_scan_fc > 0`.
* `inputs/tests/cubed_sphere_resist.athinput` -- `seam_halo_scan_fc`,
  `seam_halo_scan_fc_i` and `faces_from_potential` listed at their defaults, so a
  command line can override them (a command line cannot ADD a parameter).

Nothing else.  `git diff` is two files.

## The measurement: a STATIC face-centred halo scan

`iprob = 11` is `B = b0c*(-y, x, 0) + B_unif`: a smooth, globally defined, divergence-free
CARTESIAN field.  Its exact face-normal projection in a GHOST face is therefore known with
no reference to the neighbour panel -- evaluate it at that face's OWN chart-continued
`(xi, eta)` and project on the chart-continued panel normal there -- which is exactly what
the seam field transform plus the along-seam resample are supposed to produce.  At
`time/nlim = 0` the only thing that has touched the ghosts is one boundary exchange.

Run it with `problem/faces_from_potential = false`: the active faces are then the analytic
projections to the last bit, so the ghost error measures the halo and nothing else.  (The
default `true` builds the faces as the discrete curl of a vector potential, which makes
`div B` round-off but puts an O(h^2) discretisation offset into the ACTIVE faces that
would floor the scan.  That default is what the dynamic gates below run with.)

```bash
cmake -B build_cs_fc -D PROBLEM=cs_test -D CMAKE_BUILD_TYPE=Release && make -C build_cs_fc -j8
B=build_cs_fc/src/athena; I=inputs/tests/cubed_sphere_resist.athinput
for mb in 32 16 8; do      # 1x1, 2x2, 4x4 blocks per panel at n = 32 per panel edge
  $B -i $I -d out mesh/nx1=8 mesh/nx2=32 mesh/nx3=32 meshblock/nx1=8 \
     meshblock/nx2=$mb meshblock/nx3=$mb time/nlim=0 \
     problem/seam_halo_scan_fc=1 problem/faces_from_potential=false
  rm -rf out
done
```

It reports, per COMPONENT (the three face components are staggered differently along the
seam, so a defect need not hit them alike):

* per block, per seam face, per ghost layer, the max `|ghost - exact|` and the
  ALONG-SEAM INDEX where it occurs -- which localises a defect to a block end;
* every tangential ghost face binned by what fills it, with the same five categories as
  the cell-centred scan.

`<problem>/seam_halo_scan_fc_i = -1` moves the scan from the middle ACTIVE radial shell to
a radial GHOST shell, which is where the `x1x2x3` CORNER buffers (slots 48-55) live.  Note
that the scan deliberately stops at face index `je+ng` / `ke+ng`: a component staggered
along a ghost direction has one more face than there are ghost cells, and that outermost
face is filled by the x2/x3 FACE buffers but not by the x1-EDGE ones, so including it
reports a stale value as a halo error.  The established face-centred ghost scan in
`CSTestResistCheck` stops at the same place.

### global max |ghost - exact| over every seam ghost face (`|B|` is O(0.1-0.25))

| n per panel edge | 1 x 1 | 2 x 2 | 4 x 4 |
| --- | --- | --- | --- |
| 32 | 3.1117e-5 | **3.1117e-5** (exact match) | **5.7403e-5** (1.85x) |
| 64 | 1.3440e-5 | **1.3440e-5** (exact match) | **1.3440e-5** (exact match) |
| 128 | 1.2755e-6 | **1.2755e-6** (exact match) | **1.2755e-6** (exact match) |

Compare the cell-centred table before `94c7165d`: 4.103e-6 / 4.671e-5 / 6.950e-4 at
n = 32, i.e. a 137x degradation at 4 x 4.  Here it is 1.85x at n = 32 and nothing at all
at n = 64 and n = 128.

### the bin table, n = 32, per component (HEAD; there is no "before")

max `|ghost - exact|` on the middle ACTIVE radial shell:

| category | component | 1 x 1 | 2 x 2 | 4 x 4 |
| --- | --- | --- | --- | --- |
| face ghost, same-panel neighbour | all three | -- | **0 (exact)** | **0 (exact)** |
| face ghost, SEAM | `b.x1f` | 1.2624e-5 | 1.2624e-5 | 1.2624e-5 |
| face ghost, SEAM | `b.x2f` | 2.4239e-5 | 2.4239e-5 | 5.4740e-5 |
| face ghost, SEAM | `b.x3f` | 3.1117e-5 | 3.1117e-5 | 5.7403e-5 |
| x2x3 edge, both flanks same panel | all three | -- | **0 (exact)** | **0 (exact)** |
| x2x3 edge, a SEAM flank | `b.x1f` | -- (none exist) | 8.7268e-8 | 9.4452e-6 |
| x2x3 edge, a SEAM flank | `b.x2f` | -- | 1.3745e-6 | 2.4239e-5 |
| x2x3 edge, a SEAM flank | `b.x3f` | -- | 3.1761e-6 | 4.2522e-6 |
| x2x3 edge at a CUBE VERTEX | `b.x1f`/`x2f`/`x3f` | 5.4814e-6 / 3.1019e-6 / 1.6991e-5 | identical | identical |

The same-panel ghosts are EXACT to the last bit, as they must be (a plain index copy on
one chart), which is the control that says the scan measures what it claims.

On a radial GHOST shell (`mesh/nx1=16 meshblock/nx1=8`, `seam_halo_scan_fc_i=-1`), where
the doubly-tangential ghosts are served by the CORNER slots 48-55, the global max is
1.7499e-5 / 1.7499e-5 / 5.8160e-5 for 1 x 1 / 2 x 2 / 4 x 4 -- the same picture.

### what the 4 x 4 residual at n = 32 is

It is the extrapolation at an interior block end, the one item `tests_seam4/README.md`
lists under "What is NOT fixed", and it converges at THIRD order, so it disappears with
resolution rather than persisting:

| bin, 4 x 4 | n = 32 | n = 64 | n = 128 | order |
| --- | --- | --- | --- | --- |
| x2x3 edge, seam flank, `b.x1f` | 9.4452e-6 | 1.1528e-6 | 1.4935e-7 | 3.04 / 2.95 |

At 1 x 1 the same bin does not exist and the face bin is flat.  Removing the residual
entirely needs the receiver-side restructure the cell-centred note describes, not a change
here.

## The negative controls: the guards are load-bearing, not vacuous

"The code already has the guard" is only worth anything if the guard actually fires.  Each
control is a one-line temporary edit to `src/bvals/bvals_fc.cpp`, measured on the same
static scan, then reverted (`bvals_fc.cpp` is byte-identical to HEAD at the end).

| control | what it breaks | n = 32, 1 x 1 | n = 32, 2 x 2 | n = 32, 4 x 4 |
| --- | --- | --- | --- | --- |
| -- | HEAD, global max | 3.1117e-5 | 3.1117e-5 | 5.7403e-5 |
| **A** | clamp applied to EXTRAPOLATING stencils too (`fpos` guard dropped) | 3.1117e-5 (unchanged) | 3.1117e-5 | **9.3883e-3 (164x)** |
| **B** | slots 40-55 back to a plain copy (`cs_seam` never set) | not run | 3.1117e-5 | 5.7403e-5 |
| **C** | corner slots 48-55 only, back to a plain copy | not run | 1.7499e-5 (radial-ghost shell) | 5.8160e-5 (radial-ghost shell) |

Control A was run at 1 x 1 as well and is unchanged there, which is the statement that the
guard is inert at one block per panel (the pull is inward at a panel end, so the stencil
never extrapolates).  B and C were not run at 1 x 1 because there is nothing to measure:
with one block per panel every cross-panel x2x3 edge IS a cube vertex and the exchange is
skipped, so the branch they disable never fires.

Controls B and C do not move the GLOBAL max, because the face bin dominates it; they move
the edge/corner bins, which is where their machinery lives:

| control | bin | 2 x 2 | 4 x 4 |
| --- | --- | --- | --- |
| HEAD | x2x3 edge, seam flank (`x3f` / `x2f`) | 3.1761e-6 | 2.4239e-5 |
| **B** | same bin | **2.7251e-3 (858x)** | **1.0691e-2 (441x)** |
| HEAD | radial-ghost corner, seam flank (`x3f` / `x2f`) | 3.2241e-6 | 2.2067e-5 |
| **C** | same bin | **2.7598e-3 (856x)** | **1.0792e-2 (489x)** |

So:

* **Defect 1 (the monotonicity clamp applied to an extrapolant) does NOT exist for B.**
  `seamval()` in `bvals_fc.cpp` guards it with `dd1*dd2 > 0.0 && fpos >= 0.0 &&
  fpos <= 2.0`, where `fpos = pos - bs` comes out of `cs_seam_stnc` -- a table built PER
  COMPONENT (`kd = coar*3 + vv`) and PER SEAM ORIENTATION (`sidx`), with `bhi` carrying
  the `+1` of the staggered component (`bhi = ke + (vv == 2 ? 1 : 0) - 2`).  Control A
  degrades all three components on both orientations, so the guard is reached for each of
  them.  The MONOTONE half of the guard is also still in place here (`bvals_cc.cpp` had to
  drop it for the flat-topped blast cap); the MHD blast gate below shows that costs
  nothing for B.
* **Defect 2 (the x2x3 edge buffers getting no resample) does NOT exist for B** -- and
  neither does its untreated cousin: `bvals_fc.cpp` sets `cs_seam` for `n >= 40 && n < 56`,
  i.e. the x2x3 EDGES (40-47) **and the CORNERS (48-55)**, from which flanking face is
  cross-panel, using the same decode as the cell-centred fix.  `bvals_cc.cpp` still stops
  at 48.
* `bs` is already clamped to the SOURCE's ACTIVE range (`ks_`/`ke_`), which is
  `94c7165d`'s third hunk.

## Gates

### (a) static scan -- above.  1 x 1 == 2 x 2 exactly at every resolution; 4 x 4 equal at n >= 64

### (b) 1 x 1 bitwise against the reference binary

Reference = the UNMODIFIED tree built in the same `build_cs_fc`, binary copied aside
before any edit.  `time/nlim = 20`, a `bin` dump every step, 1 x 1 blocks per panel,
`diff -r` over the dump directories:

| test | input | result |
| --- | --- | --- |
| hydro cs blast | `cubed_sphere_blast.athinput` | **BITWISE IDENTICAL**, 22 dumps |
| MHD cs blast | `cubed_sphere_mhd_blast.athinput` | **BITWISE IDENTICAL**, 22 dumps |

(Trivially so -- no solver source changed -- but it is the gate that says the `cs_test.cpp`
addition is inert when the flag is off.)

### (c) cs MHD dynamic tests, 1 x 1 vs 2 x 2 vs 4 x 4 blocks per panel

All at HEAD, i.e. these are simultaneously the "before" and the "after" table.  The
question they answer is how much the number of MeshBlocks a panel is cut into still
changes an MHD answer.  All use `faces_from_potential` at its default `true`, so `div B`
is round-off by construction and stays there.

#### `cubed_sphere_mhd_conv.athinput` (iprob 9, rigid rotation carrying a tilted uniform field), `tlim = 1.0`, n = 32 per panel edge

| blocks/panel | cycles | `L1(B)` | `Linf(B)` | `L1(v)` | `L1(p)` |
| --- | --- | --- | --- | --- | --- |
| 1 x 1 | 157 | 2.395101e-4 | 1.384850e-3 | 6.221145e-4 | 2.620363e-4 |
| 2 x 2 | 157 | 2.398218e-4 | 1.387510e-3 | 6.219153e-4 | 2.620955e-4 |
| 4 x 4 | 157 | 2.401527e-4 | 1.479317e-3 | 6.226583e-4 | 2.620344e-4 |

`L1(B)` agrees with 1 x 1 to 0.13 % (2 x 2) and 0.27 % (4 x 4); `L1(v)` and `L1(p)` to
0.03 % and 0.02 %.  The one column that moves is `Linf(B)` at 4 x 4, +6.8 %, which is the
same n = 32 interior-block-end extrapolation the static scan localises -- an `Linf` over
the whole mesh picks up exactly the one worst ghost.

The seam flux gate is unmoved by the layout, and the same-panel block faces telescope
exactly, which is the control that the extra internal boundaries themselves are clean:

| blocks/panel | seam `dM/dt` | seam `dE/dt` | NULL (same-panel block faces) sum / max / count |
| --- | --- | --- | --- |
| 1 x 1 | -1.08843e-15 | -2.85199e-15 | 0 / 0 over 0 |
| 2 x 2 | -1.08911e-15 | -2.84708e-15 | **0.00000e+00 / 0.00000e+00** over 3072 |
| 4 x 4 | -1.08899e-15 | -2.84853e-15 | 9.97466e-18 / 3.68629e-18 over 9216 |

#### `cubed_sphere_resist.athinput` (iprob 11, the resistive current-density test), `eta = 0.5`, `tlim = 0.2`, n = 32 per panel edge

| blocks/panel | **max `|sum A.B|` (div B), seam / interior** | Ohmic heating ratio | EMF `L1`, panel-edge ring / interior | static-B halo `Linf`, tangent ghost | `r x EDGE` seam | `r x CORNER` seam |
| --- | --- | --- | --- | --- | --- | --- |
| 1 x 1 | **2.470e-14 / 2.257e-14** | 1.002451 | 3.4336e-4 / 3.1007e-4 | 1.0109e-4 | 3.1250e-5 | 1.5776e-5 |
| 2 x 2 | **2.534e-14 / 1.962e-14** | 1.002451 | 3.1011e-4 / 3.1873e-4 | 1.0109e-4 | 3.1250e-5 | 1.7249e-5 |
| 4 x 4 | **2.836e-14 / 2.145e-14** | 1.002445 | 3.1621e-4 / 3.1495e-4 | 1.0111e-4 | 3.1217e-5 | 2.8883e-5 |

`div B` is at ROUND-OFF on every layout, on the seam as in the interior -- the CT/seam-EMF
path is untouched by anything here, as it must be.  The Ohmic heating ratio agrees to 6
digits (1 x 1 and 2 x 2 bitwise-equal in the printed digits, 4 x 4 in the 6th), the evolved
static-field halo `Linf` to 5 digits, and the EMF `L1` in the panel-edge ring is if anything
BETTER at 2 x 2 and 4 x 4 than at 1 x 1 (3.10e-4 / 3.16e-4 against 3.43e-4), because the
ring is a smaller fraction of the mesh.  The only column that grows with the number of
blocks is the CORNER-seam ghost `Linf` (1.58e-5 -> 1.72e-5 -> 2.89e-5, 1.8x over 4 x 4), the
same interior-block-end extrapolation as everywhere else in this note.

#### `cubed_sphere_mhd_blast.athinput` (iprob 12, a CLOSED MHD blast through a seam and a cube vertex), `tlim = 0.3`, contrast 100

| blocks/panel | cycles | total energy | min p | symmetry `Linf` (all / within 2 cells of a seam) |
| --- | --- | --- | --- | --- |
| 1 x 1 | 278 | 1.08771595925670e+2 | 1.483322 | 2.4637 / 2.1934 |
| 2 x 2 | 278 | 1.08771595317750e+2 | 1.483322 | 2.4637 / 2.4118 |
| 4 x 4 | 279 | 1.08771596424242e+2 | 1.483322 | 2.4637 / 2.4512 |

No NaN on any layout; the cycle counts agree to one (which means the `dt` history agrees),
the total energy to 8 digits and `min p` to all 7 printed digits.  This is also the gate
on the MONOTONE half of the clamp guard, which `bvals_fc.cpp` keeps and `bvals_cc.cpp`
had to drop: for the face-centred field it costs nothing here.

## Do the CORNER slots (48-55) matter for MHD?

**Yes, and `bvals_fc.cpp` already handles them.**  They are the ghosts that are ghost in
all three directions, so they only exist where a radial MeshBlock boundary meets a
tangential one; the static scan sees them with `mesh/nx1=16 meshblock/nx1=8` and
`seam_halo_scan_fc_i=-1`.  Control C, which puts ONLY slots 48-55 back to a plain copy,
degrades that bin from 3.2e-6 to 2.8e-3 at 2 x 2 and from 2.2e-5 to 1.1e-2 at 4 x 4 --
the same 400-900x the x2x3 edges would cost.  So the branch is doing real work and
removing it would be a real defect.

`bvals_cc.cpp` still stops at slot 48, i.e. the CELL-CENTRED corners across a seam are
still a plain copy.  `tests_seam4/README.md` left that open on the grounds that no gate
there reads them.  The measurement above says the face-centred equivalent is worth
400-900x, so the cell-centred one is the obvious next thing to measure -- it is NOT part
of this note, and the fix would be a two-line mirror of the `bvals_fc.cpp` decode.

## What was NOT done, and why

* **No change to `src/bvals/bvals_fc.cpp`.**  The scan found neither of the two defects it
  was sent to find, and the only residual it does find (the interior-block-end
  extrapolation at n = 32) is third order and is explicitly the item the cell-centred note
  left open as needing a receiver-side restructure.  Changing the file to chase a third-
  order residual that vanishes by n = 64 would risk the 1 x 1 bitwise property for nothing.
* The `tst/test_suite` cubed-sphere tests: `grep -rl "cubed_sphere\|cs_test" tst/test_suite`
  gives `rad/test_rad_cs_raddiff_cpu.py` and `rad/test_rad_cs_implicit_ang_cpu.py` (the two
  `fofc_cs` tests need a DEFAULT build).  Only the raddiff one was run here, and note that
  by the time the gates finished ANOTHER session had uncommitted edits in
  `src/{CMakeLists.txt, eos, mhd, outputs, pgen/pgen.cpp, utils/two_stream_rt.hpp}` in this
  same working tree, so a test-suite rebuild measures that tree and not this change alone.
  Read the result with that caveat.

| test file | result |
| --- | --- |
| `rad/test_rad_cs_raddiff_cpu.py` | **1 passed in 120 s** |
| `rad/test_rad_cs_implicit_ang_cpu.py` | not run (267 s in the cell-centred note, plus a build) |
| `hydro/test_hydro_fofc_cs_cpu.py`, `mhd/test_mhd_fofc_cs_cpu.py` | not run (need a DEFAULT build) |

The raddiff build directory was pre-populated with `cmake` + `make -j8` before pytest ran,
because the test's own build step is `make -j $(os.cpu_count())`, which is more than this
login node's memory cap allows; it was deleted afterwards (inode quota).

## Files

* `src/pgen/cs_test.cpp` -- `CSTestSeamHaloScanFC`, enrolled for `iprob = 11` when
  `<problem>/seam_halo_scan_fc > 0`, in place of `CSTestResistCheck`.
* `inputs/tests/cubed_sphere_resist.athinput` -- the three parameters at their defaults.
* `src/bvals/bvals_fc.cpp` -- **unchanged**.  The guards this note verifies live at
  lines 186-200 (the stencil table: `blo`/`bhi` clamped to the SOURCE's active range, with
  the `+1` of the staggered component, and `fpos` stored), 446 (`cs_seam` for slots
  40-55), and 831 (`dd1*dd2 > 0.0 && fpos >= 0.0 && fpos <= 2.0`).
