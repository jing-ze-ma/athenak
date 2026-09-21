# A standing PER-REGION error gate for the cubed sphere

`run_regions.sh` + `regions.py` measure the cubed-sphere error **separately in the panel
interior, next to a panel seam, and next to a cube vertex**, at three resolutions, and
report the convergence order of each region. One command reproduces the whole table:

```bash
cmake -B build_cs_reg -D PROBLEM=cs_test -D CMAKE_BUILD_TYPE=Release
make -C build_cs_reg -j8
tests_cs_regions/run_regions.sh all      # gate | matrix | wide | km | long | table
```

Serial CPU, no MPI, no Slurm, at most two concurrent solvers. Logs land in `logs/`, the
tables in `tables/regions.txt`, disposable solver output in `out/` (delete it afterwards).
`run_regions.sh table` rebuilds the tables from logs that already exist.

## The three tests

All three are existing `cs_test` problems — no new solver and no new physics. Each has an
exact solution that is known at the final time, which is what makes a per-region *error*
meaningful rather than a per-region *difference*.

| name | iprob | input | what it is | what the error is |
| --- | --- | --- | --- | --- |
| `strat` | 13 | `inputs/tests/cs_regions_strat.athinput` | isothermal hydrostatic atmosphere **at rest** under a constant inward `g` (a user source term), carrying a uniform force-free Cartesian field; beta runs 2.0 at the base to 1.3e-2 at the top | the exact solution is the initial condition forever and `v == 0`, so **`|v|` IS the error** — no reference solution to argue about |
| `rot` | 9 | `inputs/tests/cs_regions_rot.athinput` | rigid rotation `v = Omega zhat x r` (the pressure supplies the centripetal force) carrying a uniform field that precesses with it | exact at **every** time (`B(t) = Rz(Omega t) B0`), so this measures the *evolved* solution; it is the advection-across-seams-and-over-a-vertex test |
| `loop` | 11 | `inputs/tests/cs_regions_loop.athinput` | purely azimuthal `B = b0c(-y,x,0)` seeded from its **vector potential** (`faces_from_potential`), balanced by `p0 - b0c^2 R^2`; `B.rhat = 0`, so reflecting radial walls close the budget | static exact solution; the error is the evolved **face** field, which is what CT evolves |

`rot` and `loop` are copies of the stock `cubed_sphere_mhd_conv` / `cubed_sphere_toroidal`
inputs with three knobs *declared* (`<mhd>/cs_wellbalanced_src`, `<problem>/conv_errors`,
`<problem>/conv_nband`), because AthenaK refuses a command-line override of a parameter
the input file never mentions. `cs_regions_strat` additionally declares the K&M
well-balancing knobs. Nothing else differs from the stock files except `nghost = 3`.

## How the regions are defined

Classified by **angle**, not by block index, so the answer does not depend on how many
MeshBlocks span a panel. A cell (or face) is *near an edge* when it lies within `NBAND`
cells of `|xi| = pi/4` or of `|eta| = pi/4`:

* **CUBE VERTEX** — near in **both** xi and eta (codim 2)
* **panel SEAM** — near in exactly one (codim 1)
* **panel INTERIOR** — neither

`NBAND` is `<problem>/conv_nband`. The gate **scales it with the resolution** (2 / 4 / 8
at n = 16 / 32 / 64, i.e. `NBAND32 = 4` cells at n = 32) so the band has a fixed
*physical* width. A fixed **cell** count is a region that shrinks as the grid refines, and
an L1 taken over a shrinking region is not a convergence rate — that is the trap that once
read the cs seam as first order. Set `NBAND32` in the environment to change it.

This is *why* a domain norm cannot answer the question: the seam is a codim-1 set and the
vertex a codim-2 one, so a region that is first order but a few cells wide holds a
vanishing share of the total L1, and the domain norm converges at the full rate while the
region where the error lives does not.

The cell-centred split (`CSTestConvErrors`) already existed in `src/pgen/cs_test.cpp`. The
only source change made for this gate is that the **face-field** split in
`CSTestResistCheck` now reads the same `<problem>/conv_nband` instead of a hard-coded 2,
so both tables move together when the band is scaled.

## What the gate reports

Per test, per configuration, per resolution, and for each of the three regions: L1 and
Linf of `v`, of the **tangential** part of `v` alone, of `p`, and of `B`; the measured
convergence order between successive resolutions; the VERTEX/INTERIOR and SEAM/INTERIOR
error ratios; `max|v|/c_s` for the at-rest atmosphere; and the mean `dt` and wall clock,
so a "fix" that costs the time step shows up in the same table.

**The tangential velocity is the instrument on `strat`.** The dominant error of a
stratified problem is the *radial* hydrostatic imbalance — a cell-centred gravity source
does not exactly cancel the finite-volume pressure gradient — and it is the same in every
region, so it buries any vertex signal in the combined norm. That imbalance drives
`v.rhat` only; a defect in the *angular* operators drives tangential flow, of which the
exact solution has none.

**Fixed physical time, not a fixed cycle count**, and the time is per test (`strat` 0.04,
`rot` 0.25, `loop` 0.25), chosen so the *coarsest* grid still takes about 20 steps. At
three or four steps the answer is the start-up transient and every region comes out near
first order whatever the scheme does. `strat` is Alfven-limited at the top of the
atmosphere (`rho = e^-5 d0`), which is why its time is so much shorter for the same step
count.

## What actually honours `<mhd>/reconstruct` on the cubed sphere

Verified in `src/mhd/mhd_fluxes.cpp`, not assumed:

* the **radial (x1)** sweep on the cubed sphere always takes `GridPiecewiseLinearX1`, the
  position-aware PLM, whatever `reconstruct` says — `x1v` is the volume centroid on this
  grid, so the plain uniform stencil is off-centre even without a radial stretch
  (`const bool str_r1_ = pmy_pack->pmesh->use_cubed_sphere;`, line ~132);
* the **angular (x2, x3)** sweeps do honour `reconstruct` (`dc | plm | ppm4 | ppmx |
  wenoz`), because the `use_spherical_polar` guard on their `Grid*` branches is false on
  the cubed sphere.

So on the cs, `reconstruct` is an **angular-only** knob. (On spherical polar all three
directions ignore it, which is the note this started from.) `ppm4` and `ppmx` call the
same `PiecewiseParabolicX*` and differ only in the extremum-preserving limiter flag.
`ppm4`/`ppmx`/`wenoz` need `nghost >= 3`, so every run in the gate uses `nghost = 3`,
including the `plm` baseline, to keep the halo width out of the comparison.

## The two well-balancing levers are not the same thing

* **`<mhd>/cs_wellbalanced_src`** (commit 88064f67) rebuilds the *geometric source term*
  as the face sum `(1/V) sum_f A_f (T_cell.nhat_f).e_i(x_f)` — the same sum the flux
  divergence forms — so the two cancel to round-off whenever the face states equal the
  cell state. It attacks the **angular** flux-divergence-versus-source cancellation, which
  is the term that degrades at a vertex. Its cached geometry (`Coordinates::wb_geom`) is
  *mesh* geometry, rebuilt automatically when the pack changes; there is **no staleness
  knob and no staleness risk**.
* **`<mhd>/wellbalance_dynamic` + `wb_x1` + `wb_option` + `wb_cache_every`** is the
  Kappeli & Mishra hydrostatic *reconstruction* in the **radial** direction. It is a
  different scheme attacking a different term, and it needs a gravitational potential from
  the problem generator. `wb_cache_every` belongs to **this** one, not to the cs source.
  `cs_test` did not supply a potential, so this arm needed one addition to
  `src/pgen/cs_test.cpp`: for iprob 13 the gravity is a constant inward `g`, so
  `Phi = g*(r - r0)` exactly, filled into `phicc0` / `phi0.x1f` at the same `x1v` / `xx1f`
  the reconstruction reads. The gravity source term is unchanged.

**Known gap, not fixed here:** `cs_wellbalanced_src` is read only from the `<mhd>` block
(`src/coordinates/coordinates.cpp` ~line 44); the `<hydro>` branch reads
`sp_wellbalanced_src` but not `cs_wellbalanced_src`. A **pure-hydro** cubed-sphere run
therefore cannot switch the well-balanced geometric source on at all. Every test here is
MHD, so the gate is unaffected, but the one-line asymmetry is real.

## Results

50 runs, zero failures, ~50 min wall on one login node at two concurrent solvers. Full
output in `tables/regions.txt`; the rows below are the ones that carry the conclusion.

### 1. The baseline gate: is the vertex a lower-order region?

`plm`, `cs_wellbalanced_src = false`, `nghost = 3`, fixed physical time.

**`loop` (iprob 11, face field vs the exact static one) — yes, unambiguously:**

| region | n=16 | n=32 | order | n=64 | order | VTX/INT |
| --- | --- | --- | --- | --- | --- | --- |
| INTERIOR | 5.601e-05 | 1.376e-05 | 2.03 | 3.441e-06 | 2.00 | 1.00 |
| SEAM | 7.611e-05 | 1.973e-05 | 1.95 | 5.246e-06 | 1.91 | 1.52 |
| **VERTEX** | 9.810e-05 | 3.312e-05 | **1.57** | 1.054e-05 | **1.65** | **3.06** |

The interior is clean second order, the seam is essentially second order (1.91), and the
vertex is ~1.6. The VERTEX/INTERIOR ratio **grows** with refinement — 1.75, 2.41, 3.06 —
which is the signature of a genuinely lower-order region rather than a larger constant.

**`strat` (iprob 13, spurious tangential velocity at rest):** every region is second order
in L1 (INT 1.93/2.09, SEAM 1.95/2.00, VTX 1.82/2.01), but the vertex carries **2.5-2.9x**
the interior error at every resolution and the ratio does not shrink. `Linf` is only ~1.1
everywhere, including the interior, so the worst cell is not a vertex effect.

**`rot` (iprob 9):** all three regions degrade **together** (L1(v) orders INT 1.60/1.30,
SEAM 1.65/1.41, VTX 1.75/1.48) and VTX/INT actually falls, 1.49 -> 1.35 -> 1.19. L1(p) on
the same runs is a clean 2.00 in all three regions. So `rot` shows **no** vertex-specific
defect; what it shows is that the velocity/field channel of the rigid-rotation equilibrium
is not second order at these settings (radial `user` BC and joint refinement are the
suspects). This is an open item, not a vertex finding.

### 2. (a) The well-balanced geometric source, `cs_wellbalanced_src`

**On its own it does not reduce the vertex error, and on `strat` it makes it slightly
worse.** Error ratios against the `plm`/off baseline at the same resolution:

| test, n | L1(v_t) INT | SEAM | VERTEX |
| --- | --- | --- | --- |
| strat, 16 | 1.123 | 0.987 | 1.084 |
| strat, 32 | 1.119 | 0.988 | 1.087 |
| strat, 64 | 1.141 | 0.970 | 1.070 |
| rot, 32 | 1.054 | 1.002 | 1.031 |
| loop, 32 (L1 face B) | 1.030 | 1.022 | 1.009 |

Consistent across all three resolutions: about +8 % at the vertex on `strat`, neutral on
`rot` and `loop`. It does help when combined with a high-order angular reconstruction
(next section), where it nearly halves the *interior* tangential error — but the vertex is
not where it pays. Orders are unchanged. `max|v|/c_s` on `strat` goes the wrong way by
~4 % (1.865e-02 vs 1.617e-02 at n=16).

The `cs_wellbalanced_src` geometry cache holds **mesh** geometry only and is rebuilt when
the pack changes, so **there is no staleness knob and no stale-cache failure mode** for
this lever. `wb_cache_every` belongs to the K&M scheme, measured in section 4.

### 3. (b) Reconstruction: this is the lever that works

At n=32 on `strat`, ratio to the `plm`/off baseline (lower is better):

| config | L1(v_t) INT | SEAM | **VERTEX** | Linf(v_t) INT | SEAM | **VERTEX** |
| --- | --- | --- | --- | --- | --- | --- |
| plm / wb ON | 1.119 | 0.988 | 1.087 | 0.993 | 0.970 | 1.041 |
| **ppm4** | 0.807 | 0.684 | **0.529** | 0.808 | 0.902 | **1.020** |
| **ppmx** | 0.537 | 0.398 | **0.315** | 0.088 | 0.107 | **0.146** |
| **wenoz** | 0.537 | 0.397 | **0.314** | 0.088 | 0.113 | **0.160** |
| ppmx + wb ON | 0.297 | 0.222 | **0.271** | 0.057 | 0.099 | **0.150** |
| **wenoz + wb ON** | 0.297 | 0.222 | **0.272** | 0.058 | 0.103 | **0.160** |

* `ppmx` and `wenoz` are the same to three digits on this problem and cut the **vertex**
  tangential L1 by **3.2x** and the vertex `Linf` by **6.3-6.8x** — the 6.5x of the old
  worst-cell note, reproduced here as a region norm rather than a single cell.
* `ppm4` (the non-extremum-preserving limiter) gets **half** the L1 benefit and **none** of
  the Linf benefit, so the gain is the limiter behaviour at the skewed cells, not the
  higher formal order.
* `cs_wellbalanced_src` on top of `wenoz` buys a further **1.8x in the interior and at the
  seam but only 1.16x at the vertex**; its benefit is not a vertex benefit.
* Orders are preserved: `strat` `wenoz`+wb gives 1.98 / 1.92 / 1.91 (INT/SEAM/VTX).

**Where it does not help.** On `rot` at n=32 `ppm4` is actively **worse** — vertex L1(v_t)
1.535x and vertex Linf(v_t) **2.367x** the `plm` baseline; `wenoz` there gives only
0.84x/0.90x. On `loop` no reconstruction helps the vertex L1 at all (`wenoz` 1.076x,
`ppmx` 1.109x, `ppm4` 1.097x); only the vertex `Linf` improves, and only under `wenoz`
(0.894x). The 3.2x on `strat` is specific to the at-rest stratified problem.

### 4. (c) The K&M radial well-balancing, and what a stale cache does

`wellbalance_dynamic` + `wb_x1` + `wb_option = isothermal` on `strat` (the potential
`Phi = g*(r - r0)` supplied by `cs_test`):

| n | L1(v_t) VERTEX, baseline | K&M | max\|v\|/c_s baseline | K&M |
| --- | --- | --- | --- | --- |
| 16 | 4.871e-04 | 4.972e-04 | 1.617e-02 | 1.545e-02 |
| 32 | 1.377e-04 | 1.367e-04 | 7.249e-03 | 7.064e-03 |

**No effect, in either direction, beyond a few per cent.** The likely reason is visible in
the setup rather than the numbers: the reconstruction is made hydrostatically consistent
but the gravity **source** is still the plain cell-centred `-rho g` of `CSTestGravSrc`, so
the two halves are still discretised differently and nothing cancels. Making this lever
work would need the matched source term — i.e. new solver code, which is outside this
task. Reported as measured, not as a recommendation.

`wb_cache_every` over a full radial sound crossing (n=16, t=1.20, 460 cycles):
`Linf(v)` = 2.407878e-02 (every stage), 2.407721e-02 (`=1`), 2.406320e-02 (`=10`) — a
**0.06 % spread**. On a problem whose background does not evolve, a stale cache costs
nothing measurable, which is what the scheme's own comment predicts. This gate cannot
detect the staleness failure that matters on an evolving stellar envelope.

### 5. Spurious velocity at rest, after ~1 radial sound crossing

`strat`, n=16, t = 1.20 (dr/c_s = 1.22), 460 cycles, `max|v|/c_s`:

| config | max\|v\|/c_s | vs baseline |
| --- | --- | --- |
| plm, wb OFF (baseline) | 5.637e-02 | 1.00 |
| plm, wb ON | 4.911e-02 | 0.87 |
| K&M radial (any cache setting) | 5.894-5.898e-02 | 1.05 |
| **wenoz, wb OFF** | **1.134e-02** | **0.20** |
| wenoz, wb ON | 1.849e-02 | 0.33 |

`wenoz` alone gives a **5x** lower spurious velocity after a sound crossing. Note the sign
flip on `cs_wellbalanced_src`: at short time it helped when combined with `wenoz`, at long
time it **hurts** (1.849e-02 vs 1.134e-02). Short-time error ratios do not extrapolate.

### 6. What the remaining vertex error is

`strat`, n=32, **one step**, with the RHS-split diagnostics that drop one half of the
tangential momentum balance (`L1(v_t)` per region):

| what is dropped | INTERIOR | SEAM | VERTEX |
| --- | --- | --- | --- |
| nothing (the actual error) | 2.17e-06 | 4.42e-06 | 7.66e-06 |
| geometric source (`cs_diag_no_coordsrc`) | 1.307e-03 | 1.342e-03 | 1.983e-03 |
| flux divergence (`cs_diag_no_divf`) | 1.320e-03 | 1.358e-03 | 2.010e-03 |
| magnetic part of the source only (`cs_diag_no_magsrc`) | 3.78e-04 | 7.48e-04 | 1.398e-03 |

The two halves are each ~1.3e-03 in the interior and ~2.0e-03 at a vertex, and what
survives is 2.2e-06 and 7.7e-06. **The cancellation is 1 part in 600 in a panel interior
and 1 part in 260 at a cube vertex** — 2.3x worse, which is exactly the 2.5-2.9x ratio the
evolved `strat` runs report. So the residual vertex error is *not* halo error and not a
boundary-exchange error: it is the incomplete discrete cancellation between the flux
divergence and the geometric source in the skewed cells, and the magnetic part of the
stress is a large share of what has to cancel (dropping it alone leaves 1.4e-03 at the
vertex). This is also why the well-balanced *source* cannot fix it on its own: it makes the
two halves cancel exactly only when the face states equal the cell state, and a
high-order angular reconstruction is what brings the face states closer to that.

### Cost and stability

No configuration changed the time step: the cycle count at fixed `tlim` is identical
(`strat`: 16 / 32 / 65 at n = 16 / 32 / 64 for every config). Nothing went NaN, hit a
floor-driven failure, or reported a FATAL in 50 runs. `wenoz` costs about **1.28x** the
wall clock of `plm` (252 s vs 198 s at n=64 on `strat`); `cs_wellbalanced_src` costs about
**5 %** (208 s vs 198 s), consistent with the 4 % the caching commit measured.

## Conclusion, honestly

1. **The cube vertex is a real lower-order region on the MHD vector-potential test**
   (`loop`: order 1.6 against 2.0 in the interior, ratio growing 1.75 -> 3.06 with
   refinement) and a real **2.5-2.9x** error hot spot on the at-rest stratified atmosphere,
   with the order preserved there. It is **not** visible on the rigid-rotation test, which
   degrades uniformly for an unrelated reason.
2. **(a) `cs_wellbalanced_src` does not reduce the vertex error.** It is neutral on `rot`
   and `loop` and about **8 % worse** at the vertex on `strat`, at every resolution, and it
   makes the long-time spurious velocity worse when combined with `wenoz`. It does help the
   panel interior and seam markedly when paired with a high-order angular reconstruction
   (1.8x), so it is not useless — it is simply not a vertex fix. Its cache carries no
   staleness risk.
3. **(b) Reconstruction is the lever that works, and only `ppmx`/`wenoz`.** On the at-rest
   atmosphere they cut the vertex tangential L1 **3.2x** and the vertex Linf **6.3x**, and
   the spurious velocity after a sound crossing **5x**, at 1.28x the cost and with no
   change to the time step or to stability. `ppm4` gets half the L1 and none of the Linf,
   and is **worse than `plm`** at the vertex on `rot` (Linf 2.37x) — so this is a limiter
   effect, and it is not uniform across problems: on `loop` no reconstruction improves the
   vertex L1 at all.
4. **The remaining vertex error is the flux-divergence / geometric-source cancellation**
   in the skewed cells: 1 part in 260 at a vertex against 1 part in 600 in the interior,
   with the magnetic part of the stress a large share of it. That is a cell-balance
   (consistency) defect, matching the earlier "oracle" finding that the halo is innocent.
5. **Nothing measured here gets worse in stability, positivity or time step.** The only
   regressions are accuracy ones, all named above.

## Reproducing and cleaning up

```bash
tests_cs_regions/run_regions.sh all      # ~50 min, 2 concurrent, serial CPU
tests_cs_regions/run_regions.sh table    # rebuild tables/regions.txt from logs/
rm -rf tests_cs_regions/out build_cs_reg # dumps and binaries -- inode quota
```

`logs/`, `tables/` and the two scripts are the artefacts worth keeping; `out/` is not.

