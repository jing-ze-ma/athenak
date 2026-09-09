# Well-balanced restart bug: read this before restarting or judging any WB run (2026-09-09, orion)

Fixed and pushed as **5c0b98e4** on `polar-average-perf`. Pull it, rebuild, and the fix is in
at zero cost. This note is for the viper side, which cannot see orion's run directories.

## The bug

`src/hydro/hydro_fluxes.cpp` (and identically `src/mhd/mhd_fluxes.cpp`) rebuilt the
well-balanced background cache `wbq0` only when

    stage == 1 && ncycle % wb_cache_every == 0

`Mesh::ncycle` is **restored from the restart file** (`src/mesh/build_tree.cpp`), while `wbq0`
is freshly allocated -- all zeros -- in the Hydro/MHD constructor. So a restart at a cycle
that is not a multiple of `wb_cache_every` ran up to `wb_cache_every - 1` cycles with a ZERO
background: the x1 reconstruction subtracted nothing and the scheme was silently not
well-balanced for those cycles. In a body that is hydrostatic to 1e-6 that is a free-fall
kick of order g * (wb_cache_every-1) * dt, spherically coherent, at every restart that does
not happen to land on a multiple of `wb_cache_every` (9 in 10 for the usual 10).

It existed since the WB scheme was committed (140dbf9e .. 8c384e46), so **every chained WB
run made before 5c0b98e4 -- the dhj hot-Jupiter productions included -- received this kick
at each restart.**

## What it did (measured, red_giant, 6x32x32x320, wb_cache_every = 10)

| test | KE_r first hst row after the restart (before: 2.055e41 erg) |
| --- | --- |
| unfixed binary, plain restart | 1.26e44 in the first row (x614), peak 5.6e44 at cycle ~10 |
| unfixed binary, `hydro/wb_cache_every=1` on the restart line | 2.0550e41 (-0.01 %) |
| fixed binary (5c0b98e4), default cadence | 2.0550e41, identical to 4 digits |
| the lidded, fully convective production (prod11), its own unfixed binary | 1-momentum x27 in the first row, KE_r +61 % by the second |

The unfixed run that was left alone: a coherent radial shell (-1.35 km/s horizontal mean at
r/R 0.38, 200x the local rms) climbed the density gradient, went transonic under the
photosphere, and the run NaN'd 3.6e5 s after the restart in every column at once with dt
still at its normal value -- **no dt collapse, no precursor**. If a dhj run "died of NaN
everywhere in one interval with no precursor" (cs_mhd_prod at rot 41.6 is recorded that way),
check whether it was within a rotation or so of a chain restart before blaming anything else.

## The fix

`wb_cache_built` (a bool member in `Hydro` and `MHD`, false at construction) forces
`BuildWBCache` on the first flux call after start OR restart; afterwards the cadence is
exactly what it was. **No extra cost**: 3.15e7 vs 3.10e7 zone-cycles/s for the unfixed
default on the same segment. Only the stopgap `wb_cache_every=1` costs anything (2.5 % on the
red giant; more on the hot Jupiter, where the cache is a larger fraction of the cycle) -- use
it only for a run whose binary you must not change mid-chain, e.g. by adding
`hydro/wb_cache_every=1` (and `mhd/wb_cache_every=1` for MHD) to the restart command line.
Both parameters already exist in every WB input, so the override is legal on a restart.

## How to check a past run

`rg.hydro.hst` / the dhj `.hst` column 8 (1-KE) and column 4 (1-mom) across each restart
time: a clean restart continues the previous trend to ~1 %; the bug shows as a jump in the
first or second row after the restart, decaying over ~1e2 cycles but leaving a radial pulse.
On the cubed sphere the hst volumes are Cartesian, so use ratios, not absolute values.
