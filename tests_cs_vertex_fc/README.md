# Re-verification: the FC cube-vertex "x2f does not converge" claim

Session of 2026-09-22, serial CPU Release, `-D PROBLEM=cs_test`, HEAD 145edc00 (own build
dir `build_csvfx`, not the working tree that another session holds in `src/`).

```bash
./scan.sh ./athena_new
```

Runs `cs_test` iprob = 11 (`seam_halo_scan_fc`), one boundary exchange (`time/nlim=0`),
`mesh/nx2=nx3=n`, 1x1 and 2x2 MeshBlocks/panel, at n = 32/64/128, both
`mesh/cs_vertex_fill = false/true`, and **both** `problem/faces_from_potential = true`
(the input file `in_scan_fc.athinput`'s own default, which is what the original
`tests_cs_vertex/scan_vertex.sh` actually ran with -- it never overrides the parameter)
and `= false` (what `CSTestSeamHaloScanFC`'s own docstring says the scan *must* be run
with, to keep the active faces themselves exact to machine precision instead of the
discrete curl-of-vector-potential discretisation).

## Verdict: (a), a scan artefact -- not a defect in the sampled FC vertex fill

The original scan used the wrong input. `in_scan_fc.athinput` does not set
`faces_from_potential`, so it silently ran at the function's own documented default
(`true`), the ONE setting `CSTestSeamHaloScanFC`'s docstring explicitly says not to use.
Rerunning with `faces_from_potential=false` makes the "does not converge" result
disappear completely: the sampled fill converges at 3rd order and is 30-40x better than
extrapolation, at every resolution, on every component.

**`edge, CUBE VERTEX` bin, max`|ghost-exact|`, `problem/faces_from_potential=false`
(1x1 and 2x2 MeshBlocks/panel agree to 4+ digits for the sampled fill; extrapolation
differs by a few % between decompositions, expected -- its stencil reaches into
neighbouring same-panel blocks):**

| n | component | extrapolate (`cs_vertex_fill=false`) | sample (default, `=true`) | gain | sample order |
| --- | --- | --- | --- | --- | --- |
| 32  | x1f | 1.25e-16 | 8.67e-17 | -- (round-off both) | -- |
| 32  | x2f | 2.4562e-04 | 7.5569e-06 | **32.5x** | -- |
| 32  | x3f | 2.6398e-04 | 7.5569e-06 | **34.9x** | -- |
| 64  | x1f | 1.36e-17 | 5.20e-17 | round-off both | -- |
| 64  | x2f | 2.9254e-05 | 9.0006e-07 | **32.5x** | 3.07 |
| 64  | x3f | 3.5527e-05 | 9.0006e-07 | **39.5x** | 3.07 |
| 128 | x1f | 1.39e-17 | 2.78e-17 | round-off both | -- |
| 128 | x2f | 3.5157e-06 | 1.0946e-07 | **32.1x** | 3.04 |
| 128 | x3f | 3.8754e-06 | 1.0946e-07 | **35.9x** | 3.04 |

Extrapolation converges at essentially the same order (x2f: 3.07, 3.06 across the two
doublings), so sampling is a constant ~32-40x improvement, not a different regime -- the
mirror image of the CC result already in `tests_cs_vertex/README.md` (11x, one order,
same shape). `x1f` (the radial/normal component, a true scalar under the seam transform)
is round-off-identical under both fills at the vertex: there never was a defect there.

**The same bin, `problem/faces_from_potential=true` (what the original scan actually
ran), reproducing the flagged claim:**

| n | component | extrapolate | sample (default) |
| --- | --- | --- | --- |
| 32  | x1f | 2.2107e-05 | 5.5899e-05 |
| 32  | x2f | 2.7354e-04 | 2.2580e-04 |
| 32  | x3f | 2.7354e-04 | 1.2645e-04 |
| 64  | x1f | 6.6486e-06 | 1.3425e-05 |
| 64  | x2f | 2.1528e-04 | 2.2356e-04 |
| 64  | x3f | 1.2673e-04 | 1.1781e-04 |
| 128 | x1f | 1.8414e-06 | 3.2332e-06 |
| 128 | x2f | 2.2114e-04 | 2.2254e-04 |
| 128 | x3f | 1.1244e-04 | 1.1437e-04 |

Sampled x2f order n32->64->128: log2(2.2580e-4/2.2356e-4) = 0.014, log2(2.2356e-4/2.2254e-4)
= 0.0066 -- flat, reproducing the reported "order 0.01". But it is not specific to the
vertex or to sampling: the FIRST-PASS global seam max (general tangential face ghosts,
`FC SEAM-HALO GLOBAL MAX`) is *also* pinned at 2.97e-4 / 2.73e-4 / 2.72e-4 across all
three resolutions under `faces_from_potential=true`, in **every** arm including plain
extrapolation and `cs_vertex_fill=false`. That is the signature of an O(1),
resolution-independent floor in the *reference value*, not in the halo mechanics: the
active faces built by `CSTestBlastFaces` (`faces_from_potential=true`, the discrete curl
of a vector potential in the `A.rhat=0` gauge) do not converge to the scan's
`bexact()` -- the analytic face-normal projection of the continuous field -- at the cube
vertex the way they do on an ordinary panel-interior face. `x1f` also moves ~3x worse
with sampling under this setting (2.2e-5 to 5.6e-5 at n=32), exactly the "2x worse" the
brief flagged, even though `x1f` is round-off-exact under the correct setting -- direct
proof the floor sits in the reference/active-face discretisation, not in the fill.

## What this means for the earlier finding

`tests_cs_vertex/README.md`'s Gate (b) FC table is real *as measured*, but it measured
the wrong thing: `scan_vertex.sh` never overrides `problem/faces_from_potential`, so it
ran at the input file's default (`true`), the one setting `CSTestSeamHaloScanFC`'s own
header comment says will contaminate the scan with an O(h^2)-on-paper-but-not-here
offset from the vector-potential curl. Once the scan is run the way the function's
docstring requires, the sampled cube-vertex FC fill is clean: 3rd-order convergent and
30-40x better than extrapolation on both tangential components, with the radial
component exact to round-off in both arms. There is no defect to localise or fix in
`buffs_fc.cpp` / `bvals_fc.cpp`'s cube-vertex sampling path, and no implication for the
field-loop vertex order 1.6 from this mechanism -- that residual has to come from
somewhere else (the dynamic Gate (c) result in `tests_cs_vertex/README.md`, showing under
0.5% change in the *evolved* vertex error from any of these switches, already pointed the
same way).

## Fix (for the record, not applied here)

`tests_cs_vertex/scan_vertex.sh` should pass `problem/faces_from_potential=false`
explicitly (or `in_scan_fc.athinput` should set it), matching what
`CSTestSeamHaloScanFC`'s own docstring already requires. No source file needs to change;
this is a test-harness/input fix in `tests_cs_vertex/`, not `src/`.

## Caveat

Serial only (no MPI on this machine); 1x1 vs 2x2 MeshBlocks/panel agree for the sampled
fill (identical to 4+ digits) at every n tested, matching the earlier CC/whole-vertex
finding. `mesh/nx1=8 meshblock/nx1=8` (single radial block): only the tangential
`edge, CUBE VERTEX` bin was scanned, not the triply-ghost radial corner (slots 48-55),
which is outside the disputed claim.
