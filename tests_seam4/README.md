# The cubed-sphere seam halo with more than 2 x 2 MeshBlocks per panel

`tests_adi/README.md` section (d) recorded that splitting a panel into 4 x 4 MeshBlocks
degrades the cubed-sphere seam halo -- the EXPLICIT transverse operator's max residual
went 0.0895 -> 0.2618 at n = 64 and 0.2833 -> 0.6372 at n = 128, independently of block
SIZE, of `nghost` and of the solver (both are back to the 1 x 1 value after this note) -- and left it as "presumably the along-seam
resample".  It is the along-seam resample, in TWO separate places.  This is the
measurement, the cause and the fix.  The production cubed-sphere grid is 4 x 4 blocks per
panel, so it was on the path of the production.

## The measurement: a STATIC halo scan

`cs_test` `iprob = 15` (`inputs/tests/cubed_sphere_raddiff.athinput`) starts from a
spherical harmonic of the DIRECTION COSINES alone, so the exact value in a GHOST cell is
known with no reference to the neighbour panel: evaluate the same harmonic at the ghost's
OWN chart-continued `(xi, eta)`.  The gnomonic formulae are pure `tan()` and continue
analytically past `|xi| = pi/4`, so this is the value a perfect halo would produce, and
the along-seam resample is what is supposed to produce it.

`<problem>/seam_halo_scan = 1` enrols `CSTestSeamHaloScan` as `pgen_final_func`
(`src/pgen/cs_test.cpp`).  Run at `time/nlim = 0`, so the ONLY thing that has touched the
ghosts is one boundary exchange: no operator, no time integration, nothing to average
over.  It prints two things:

* per MeshBlock, per SEAM FACE and per ghost LAYER, the max `|ghost - exact|` and the
  ALONG-SEAM INDEX where it occurs (`= 2` adds the full along-seam profiles);
* a scan of EVERY tangential ghost cell, binned by what fills it: face ghost from a
  same-panel neighbour, face ghost from a seam, x2x3 edge ghost with both flanks on the
  same panel, x2x3 edge ghost with a seam flank, and the cube-vertex fill.  That second
  bin table is what separates the two defects.

```bash
cmake -S .. -B ../build_cpu_cs_test -D PROBLEM=cs_test -D CMAKE_BUILD_TYPE=Release
make -C ../build_cpu_cs_test -j
B=../build_cpu_cs_test/src/athena
I=../inputs/tests/cubed_sphere_raddiff.athinput
for mb in 64 32 16; do    # 1x1, 2x2, 4x4 blocks per panel at n = 64 per panel edge
  $B -i $I -d out_$mb mesh/nx2=64 mesh/nx3=64 meshblock/nx2=$mb meshblock/nx3=$mb \
     time/nlim=0 problem/seam_halo_scan=1 output1/dt=1e9 | grep -E "GLOBAL MAX|max ="
  rm -rf out_$mb
done
```

### global max |ghost - exact| over every seam ghost

| n per panel edge | blocks per panel | before | after | 1 x 1 reference |
| --- | --- | --- | --- | --- |
| 32 | 1 x 1 | 4.103e-6 | **4.103e-6** | -- (bitwise unchanged) |
| 32 | 2 x 2 | 4.671e-5 | **4.103e-6** | 4.103e-6 (exact match) |
| 32 | 4 x 4 | 6.950e-4 | **5.077e-6** | 4.103e-6 (1.24x) |
| 24 | 3 x 3 | 7.096e-4 | **7.638e-6** | 7.638e-6 (exact match) |
| 48 | 4 x 4 | 4.575e-4 | **2.426e-6** | 2.426e-6 (exact match) |
| 64 | 4 x 4 | 3.409e-4 | **5.639e-7** | 5.575e-7 (1.01x) |

1 x 4 and 4 x 1 cannot be run: `Mesh` requires the same number of MeshBlocks along x2 and
x3 on the cubed sphere.

### the bin table, n = 64, before -> after

max `|ghost - exact|`:

| category | 1 x 1 | 2 x 2 | 4 x 4 |
| --- | --- | --- | --- |
| face ghost, same-panel neighbour | -- | 0 (exact) | 0 (exact) |
| face ghost, SEAM | 5.58e-7 -> 5.58e-7 | 5.58e-7 -> 5.58e-7 | 5.64e-7 -> 5.64e-7 |
| x2x3 edge ghost, both flanks same panel | -- | 0 (exact) | 0 (exact) |
| **x2x3 edge ghost, a SEAM flank** | -- (none exist) | **3.50e-5 -> 4.28e-9** | **3.54e-4 -> 5.64e-7** |
| x2x3 edge at a CUBE VERTEX | 3.66e-6 (unchanged) | 3.66e-6 | 3.66e-6 |

(the face-ghost column is "after the first fix" in both entries -- before that fix it was
4.10e-6 / 4.67e-5 / 6.95e-4 at n = 32; the 5.6e-7 numbers here are n = 64.)

The same-panel ghosts are EXACT to the last bit, as they must be (a plain copy of a cell
on the same chart), which is the control that says the scan is measuring what it claims.

### localisation of defect 1, 4 x 4 at n = 32 (block = 8 cells)

Before the fix, the whole face-ghost error is in the LAST one (ghost layer 0) or TWO
(layer 1) cells of each block's along-seam segment, at ONE end of the segment, and WHICH
end is set by the block's position in the panel.  `-x2` face, layer 0, panels 0 and 1
(`s` is the along-seam index within the block, 0..7):

| gid | panel | lx3 | max err | at s |
| --- | --- | --- | --- | --- |
| 0 | 0 | 0 | 7.137e-5 | **7** |
| 2 | 0 | 1 | 6.119e-6 | **7** |
| 8 | 0 | 2 | 6.207e-6 | **0** |
| 10 | 0 | 3 | 9.521e-5 | **0** |
| 16 | 1 | 0 | 1.953e-4 | **7** |
| 18 | 1 | 1 | 6.872e-6 | **7** |
| 24 | 1 | 2 | 6.075e-6 | **0** |
| 26 | 1 | 3 | 2.036e-5 | **0** |

The blocks on the NEGATIVE half of the panel (`lx3` = 0, 1) fail at their HIGH end, the
blocks on the POSITIVE half at their LOW end: both fail at the end that faces the PANEL
CENTRE.  The along-seam profile of one block (gid 64, `-x3`, layer 1) is flat at the
5e-7 background over s = 0..5, then 1.007e-4 at s = 6 and 6.950e-4 at s = 7 -- exactly
the two cells whose resample stencil has run off the end of the block.  At 1 x 1 the same
scan is flat at 5e-7 with its (small) maxima at the CUBE VERTICES, which is the ordinary
end-of-seam behaviour.

## Root cause

Both defects are in `src/bvals/bvals_cc.cpp`, in the along-seam resample of the `do_cs`
branch of `MeshBoundaryValuesCC::PackAndSendCC`.  The mechanism they share is the shape of
the resample map.  Writing the source cell's seam-parallel angle as `a` and its
seam-normal angle as `n`, the ghost needs the source field at `atan(tan a tan|n|)`, which
for a cell in ghost layer `g` is

```
    a + delta,   delta = -(g + 1/2) dx sin 2a + O(dx^2)
```

i.e. a pull toward the PANEL CENTRE of up to `(g + 1/2)` cells, ZERO at the panel centre
`a = 0` and at the panel edges `a = +-pi/4`, MAXIMAL at `a = +-pi/8`.

### defect 1: the monotonicity clamp on an EXTRAPOLATING stencil

```c++
            int b = static_cast<int>(floor(pos + 0.5)) - 1;
            b = (b < blo) ? blo : ((b > bhi) ? bhi : b);      // blo = kl, bhi = ku - 2
...
            const Real lo0 = fmin(sv0, fmin(sv1, sv2));       // in seamval()
            const Real hi0 = fmax(sv0, fmax(sv1, sv2));
            return fmin(hi0, fmax(lo0, qq));                  // UNCONDITIONAL
```

`blo`/`bhi` are the SENDING MESHBLOCK's own range, not the panel's, so at an INTERIOR
block end the stencil runs off the block, `b` is clamped and the quadratic EXTRAPOLATES.
That end exists only when a panel is split more finely than 2 x 2: at 1 x 1 the only
block ends are the panel ends, where the pull is inward and the clamp never fires, and at
2 x 2 the one interior end sits exactly on the panel midline `a = 0`, where the pull is
zero.  At 4 x 4 the interior ends sit at `a = 0` and `a = +-pi/8`, the worst place there
is.  That is why the defect switches on at 4 x 4 and does not care about block SIZE.

Extrapolating a little is harmless (`O(dx^3)`).  Clamping an EXTRAPOLANT to the range of
its own three nodes is not: it replaces it by the nearest NODE, an `O(dx)` value.  That
is the 137x.

The face-centred twin `src/bvals/bvals_fc.cpp` already guards its clamp with "monotone
stencil AND interpolating" and carries the comment "It must also be skipped when the
resample is EXTRAPOLATING ... clamping there cost 3.5x on the smooth seam halo".  The
cell-centred copy had that guard removed -- correctly for the MONOTONE half (a blast cap
has a flat top, so a stencil astride its edge is not monotone and does need the clamp;
`iprob = 12` at contrast 100 went NaN with the monotone guard in place) -- but the
INTERPOLATING half went with it.

### defect 2: the x2x3 EDGE ghosts of a seam got no resample at all

```c++
          if (n >= 8 && n < 24) {
            cs_seam = 2;
          } else if (n >= 24 && n < 40) {
            cs_seam = 3;
          }
          // ... slots 40-47 (x2x3 edges) fall through with cs_seam = 0 = plain copy
```

with the comment "one that is ghost in BOTH tangential directions -- the x2x3 edges and
the corners -- has no single along-seam axis and stays a plain copy".  It does have one:
EXACTLY ONE of the two flanking faces of such a buffer is a panel seam, because if both
were, the buffer would be at a CUBE VERTEX and the exchange is skipped entirely
(`IsCubeVertexCorner`).  The seam normal is that face's axis and the resample runs along
the other, exactly as for the face buffer beside it.

With ONE block per panel every cross-panel x2x3 edge IS a cube vertex, so nothing was ever
wrong; the moment a panel is split, these become ordinary diagonal ghosts with a real
donor, served by a plain copy that is off by `delta` -- `O(dx)`, and again largest where
the block ends sit at `a = +-pi/8`.  **These are the ghosts the metric cross term of the
transverse operator reads**, which is why fixing only the face halo left the EXPLICIT
operator's ladder row exactly where it was (0.2618 -> 0.2617) while the implicit ones
improved 7-9x.

## The fix

`src/bvals/bvals_cc.cpp`, three hunks:

1. carry `cs_interp = (u >= -1) && (u <= 1)` out of the stencil construction and return
   the raw quadratic from `seamval()` when it is false.  The MONOTONE half of the `fc`
   guard is deliberately NOT restored.
2. set `cs_seam` for the x2x3 edge slots (40-47) from which flanking FACE neighbour is
   cross-panel.
3. clamp `b` to the SOURCE'S ACTIVE range instead of the BUFFER's.  For a face or x1-edge
   buffer the two coincide -- the along-seam extent of such a buffer IS the active range
   -- so this is bitwise for everything that resampled before; an x2x3 edge buffer is only
   `nghost` deep along the seam, and clamping to that would extrapolate from two cells
   while the source block holds the data.  The stencil reads the array, not the buffer, so
   any active cell is available.  `bvals_fc.cpp` says the same thing about its own table.

`bvals_fc.cpp` (the face-centred halo, i.e. MHD) has the SAME two defects: its `b` is
already clamped to the source's active range and its clamp is already guarded, but its
stencil still cannot leave the donor MeshBlock, and its `cs_seam` selection skips the
x2x3 edge buffers in the same way.  It is NOT touched here -- the gates in this note are
hydro/radiation ones and an MHD change needs its own measurement.

## What is NOT fixed

The residual 1.24x at n = 32, 4 x 4 (5.077e-6 against 4.103e-6) is the extrapolation
itself at an interior block end.  It is `O(dx^3)` and it vanishes at the other resolutions
tested.  Removing it entirely would need donor cells from the along-seam NEIGHBOUR
MeshBlock of the SAME panel, which the sender does not own at pack time (its own ghosts
there are a stage stale).  The clean way is to resample on the RECEIVER, after the whole
halo including the x2x3 edge buffers has been unpacked, where the ghost strip has
`nghost` extra cells at each end and the stencil can always interpolate.  That is a
restructure and is not done here.

## Gates

### bitwise, 1 block per panel

`cubed_sphere_raddiff.athinput`, `nlim = 20`, `cfl = 0.03`, a `hydro_w` dump every step,
head binary (`a8222bd3`) vs new binary, compared with `tests_adi/cmpbin.py`:

| layout | solver | result |
| --- | --- | --- |
| 1 x 1 | explicit | **BITWISE IDENTICAL**, 22 dumps |
| 1 x 1 | `adi` | **BITWISE IDENTICAL**, 22 dumps |
| 2 x 2 | explicit | differs (max `eint` 1.2e-7 = one float32 ULP of the dump) |
| 2 x 2 | `adi` | differs (max `eint` 2.3e-6) |

2 x 2 MUST differ: the caller's brief assumed it was clean, and the static scan says it
is not -- 4.671e-5 against the 4.103e-6 of 1 x 1, and 3.50e-5 in its x2x3 edge ghosts.
The fix is what removes that.  1 x 1 is bitwise for both defects by construction: the
clamp never fires there (the pull is inward at a panel end) and every cross-panel x2x3
edge is a cube vertex, whose exchange is skipped.

Everything OFF a panel seam is untouched by construction: all three hunks live inside
`if (do_cs || do_pole)` and only act when `cs_seam != 0`, which no Cartesian,
spherical-polar or polar-boundary buffer ever sets.

### the diffusion ladder (`tests_adi/measure_cs.py`, `l = 2`, `cfl = 0.0075`, `nlim = 80`)

`Linf` is the max residual of the discrete angular source against the exact
`-l(l+1) kappa (T - T0)/r^2`; `edge` is that max restricted to the outer ring of every
MeshBlock, which is the panel edge only at one block per panel.  "before" is this tree at
`a8222bd3`, "after" is the same tree with the three hunks.

#### n = 64 per panel edge, `Linf` (max residual) and `L1`

| solver | blocks/panel | before `Linf` | after `Linf` | before `L1` | after `L1` |
| --- | --- | --- | --- | --- | --- |
| explicit | 1 x 1 | 0.08951 | 0.08951 | 0.09571 | 0.09571 |
| explicit | 2 x 2 | 0.08951 | 0.08951 | 0.09570 | 0.09571 |
| explicit | **4 x 4** | **0.2618** | **0.08951** | 0.09701 | **0.09571** |
| `adi` | 1 x 1 | 0.05323 | 0.05323 | 0.02546 | 0.02546 |
| `adi` | 2 x 2 | 0.05323 | 0.05323 | 0.02547 | 0.02546 |
| `adi` | **4 x 4** | **0.4806** | **0.05323** | 0.03069 | **0.02546** |
| `sts` | 1 x 1 | 0.03504 | 0.03504 | 0.02443 | 0.02443 |
| `sts` | 2 x 2 | 0.03504 | 0.03504 | 0.02444 | 0.02443 |
| `sts` | **4 x 4** | **0.4871** | **0.03505** | 0.03006 | **0.02443** |

**At 2 x 2 and 4 x 4 every column now equals the 1 x 1 row to all printed digits**, for
the explicit operator and for `adi`; `sts` at 4 x 4 lands on 0.03505 against the 1 x 1
0.03504, i.e. the last digit.  That is the claim
the whole note is about: the number of MeshBlocks a panel is cut into no longer changes
the answer.

The EXPLICIT row is the one that identifies defect 2.  With only defect 1 fixed it stayed
at 0.2617 against 0.2618 -- the explicit operator reads the x2x3 DIAGONAL ghost for its
metric cross term, and that ghost was still a plain copy -- while `adi` improved to
0.06676 and `sts` to 0.05234.  Fixing the x2x3 edge buffers takes all three to the 1 x 1
values exactly.

#### n = 128 per panel edge

| solver | blocks/panel | before `Linf` | after `Linf` | before `L1` | after `L1` |
| --- | --- | --- | --- | --- | --- |
| explicit | 1 x 1 | 0.2833 | 0.2833 | 0.3719 | 0.3719 |
| explicit | **4 x 4** | **0.6372** | **0.2833** | 0.3724 | **0.3719** |
| `adi` | 1 x 1 | 0.1568 | 0.1568 | 0.02924 | 0.02924 |
| `adi` | **4 x 4** | **0.3337** | **0.1568** | 0.03095 | **0.02924** |

Again exact agreement with 1 x 1 in every column, `amp`, `edge` and `interior` included.


### the cubed-sphere regression tests

```bash
cd tst
python3 -m pytest test_suite/rad/test_rad_cs_raddiff_cpu.py \
                  test_suite/rad/test_rad_cs_implicit_ang_cpu.py -q
# the two FOFC cubed-sphere tests need a DEFAULT build (pgen shock_tube), reachable
# without run_test_suite.py's rebuild by putting the binary where testutils expects it:
ln -sf <a default build>/src/athena athena
python3 -m pytest test_suite/hydro/test_hydro_fofc_cs_cpu.py \
                  test_suite/mhd/test_mhd_fofc_cs_cpu.py -q
```

| test file | result |
| --- | --- |
| `rad/test_rad_cs_raddiff_cpu.py` | **passed** |
| `rad/test_rad_cs_implicit_ang_cpu.py` | **passed** (2 passed in 267 s) |
| `hydro/test_hydro_fofc_cs_cpu.py` + `mhd/test_mhd_fofc_cs_cpu.py` | **7 passed in 85 s** |

They are the only tests in `tst/test_suite` that touch the cubed sphere
(`grep -rl "cubed_sphere\|cs_test" tst/test_suite` plus the two `fofc_cs` inputs).

### the cubed-sphere BLAST, which is what the monotonicity clamp is there for

`inputs/tests/cubed_sphere_blast.athinput` (`iprob = 12`, a VERTEX-centred blast), at
pressure contrast 100 and 1000, with 1 x 1 and 4 x 4 MeshBlocks per panel:

| contrast `problem/blast_p` | blocks/panel | head | new |
| --- | --- | --- | --- |
| 100 | 1 x 1 | 278 cycles to `tlim`, no NaN | 278 cycles, no NaN |
| 100 | 4 x 4 | 278 cycles, no NaN | 278 cycles, no NaN |
| 1000 | 1 x 1 | 680 cycles, no NaN | 680 cycles, no NaN |
| 1000 | 4 x 4 | 680 cycles, no NaN | 680 cycles, no NaN |

The cycle counts are equal, which means the `dt` history is equal, which is a stronger
statement than "it did not blow up".

The `iprob = 12` NaN that motivated the unconditional clamp is a MONOTONICITY failure, not
an extrapolation one, and the monotone half of the guard stays dropped, so the blast is
unaffected -- it runs to the same cycle count as the head binary in every arm.

## Files

* `src/bvals/bvals_cc.cpp` -- the fix (three hunks).
* `src/pgen/cs_test.cpp` -- `CSTestSeamHaloScan`, the static scan, enrolled for
  `iprob = 15` when `<problem>/seam_halo_scan > 0`.
* `inputs/tests/cubed_sphere_raddiff.athinput` -- `seam_halo_scan = 0` listed so a command
  line can override it.

## Left open

* `bvals_fc.cpp` (MHD) has the same two defects in its own resample and is untouched.
* the CORNER buffers (slots 48-55, the x1x2x3 corners) are still a plain copy across a
  seam, for the same reason the x2x3 edges were.  Those ghosts are outside the ACTIVE
  radial range as well, so no operator in this note's gates reads them and there is no
  measurement to justify the change; it is the obvious next one if a gate ever sees them.
* the residual extrapolation at an interior block end (see "What is NOT fixed").
