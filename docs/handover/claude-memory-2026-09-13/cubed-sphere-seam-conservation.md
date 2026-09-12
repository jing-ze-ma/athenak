---
name: cubed-sphere-seam-conservation
description: CLOSED 2026-08-29, 985faa22 + f8940ed9 + 806a9e6f. The seam was NOT conservative; new flux_seam_cc.cpp reconciles the shared-face flux and mass and energy are now exact to ROUND-OFF (1e5-1e8 gain). Also fixed history.cpp, which weighted every cell by the CARTESIAN volume on every grid
metadata:
  type: project
---

The user's question, asked 2026-08-29: at a seam (and at the 8 cube vertices) the halo
has to interpolate -- so what guarantees the conservation laws, and what guarantees
second order there?  **A "seam" is one of the 12 cube edges projected onto the sphere,
where two panels meet; a cube vertex is where three meet.**

## THE ANSWER

**Second order: the halo guarantees it** ([[cubed-sphere-seam-interp]] resample +
[[cubed-sphere-seam-basis]] tangent-basis transform; seam <= interior, measured).

**div B: exact by construction, and by a DIFFERENT mechanism.** CT conserves
sum(area*B) for ANY edge EMF, provided every edge carries ONE value; the seam EMF
exchange enforces that ([[cubed-sphere-seam-emf]]). Structural, not accuracy.

**Mass and total energy: nothing DID guarantee them -- now FIXED (806a9e6f).** A
finite-volume update telescopes only because both neighbours multiply the SAME stored
flux by the SAME stored area. Inside a panel, and between MeshBlocks of one panel, that
holds bit for bit (ghosts are copies). At a seam the ghosts are INTERPOLATED, the two
panels' Riemann solves disagree, and the difference is mass and energy from nothing.
There was no reconciliation: `SendFlux`/`RecvFlux` restrict fluxes only when `multilevel`,
and a uniform cubed-sphere grid is not multilevel. `bvals/flux_seam_cc.cpp` now does it.

**Momentum: the question does not apply.** On the gnomonic grid the momentum equations
carry geometric source terms (src1/src2), so momentum is not a telescoping conservation
law at all -- its error is the source-term discretisation, independent of the seam. Only
mass and total energy have "exactly conservative" as a meaningful target.

## RETRACTED: "mass drifts at second order, ratio 4.00"

The parked measurement (1.667e-5 / 4.167e-6 at nx2=16/32) is wrong three times over:

1. `.hst` prints `%12.5e`, and 24.0000-23.9996 quantises to exactly ratio 4. **Fix, no
   code change: `data_format = %26.17e` in the `<output>` block.**
2. It refined nx2/nx3 with **nx1 held at 8** -- the exact trap [[cubed-sphere-seam-interp]]
   already recorded. At fixed nx1 the drift PLATEAUS (ratios 1.63, 1.14) and reflecting
   radial BCs are 8x worse and dead flat: what dominates is radial truncation error.
3. **`.hst` cannot measure conservation on this grid at all.** The user caught this.
   `HistoryOutput::Load{Hydro,MHD}HistoryData` (`src/outputs/history.cpp:123,209,309`)
   weights every cell by `dx1*dx2*dx3`, the CARTESIAN coordinate volume, and NEVER by
   `Coordinates::volume` -- so its "mass" is not sum(rho dV), and its 1/2/3-mom adds up
   components on a basis that points a different way in every cell. **This is a live bug
   for cubed-sphere AND spherical-polar, not just for this test.**

## The new instrumentation, in `cs_test.cpp` (uncommitted, +217 lines)

* **`CSTestConsSums`** -- `### CS CONSERVED SUMS`, at `%.17e`: sum(volume*u) for mass and
  energy, and momentum rotated into the FIXED Cartesian frame first
  (P = u1*rhat + u2*e1hat + u3*e2hat). Get t=0 by re-running with `time/nlim=0`.
* **`CSTestSeamFluxCheck`** -- `### CS SEAM FLUX MISMATCH`: the cell-centred twin of
  `CS SEAM FACE SINGLE-VALUED`. Finds shared seam faces BY GEOMETRY (key = physical face
  centre), and sums `area*(F_P + F_Q)` for the MASS and ENERGY fluxes -- true scalars, so
  they must agree up to the outward normal's sign with no basis transform. That sum IS
  the spurious source term. Split by distance to a cube vertex.
* Both are called from `CSTestGhostCheck` AND `CSTestConvErrors`, so iprob 3/5/6/7/8/9.
* **The NULL control is what makes it a measurement**: same-panel block-internal faces go
  in a separate bucket and must be round-off. **Measured exactly 0.0 over 6144 faces**
  (nx2=32, 4 blocks/panel). Keying, sign convention and areas are all validated.
* Bucket the vertex neighbourhood with `0.25*pi*(x2max-x2min)/indcs.nx2`, NOT
  `0.5*pi/indcs.nx2`: `indcs.nx2` counts cells per MESHBLOCK, so the naive form doubles
  the vertex set when a panel is split. The older FC gate still has the naive form.

## Measured -- hydro iprob=3, JOINT refinement (nx1,nx2)=(8,16)/(16,32)/(32,64), t=0.25

Observed drift rate from `CS CONSERVED SUMS` vs the gate's seam source:

| | 16 | 32 | 64 | ratio |
|---|---|---|---|---|
| mass, observed | 3.933e-4 | 6.355e-5 | 1.024e-5 | 6.2 |
| mass, seam gate | 2.558e-4 | 4.516e-5 | 6.149e-6 | 5.7-7.3 |
| energy, observed | 1.649e-3 | 2.567e-4 | 4.002e-5 | 6.4 |
| energy, seam gate | 1.245e-3 | 2.003e-4 | 2.773e-5 | 6.2-7.2 |

**The seam accounts for obs/gate = 1.3-1.7 of the drift, stably across resolution**, and
both converge at the same rate. Drift over fixed time is O(h^2.6). The residual factor
~1.4 is the radial BC plus the fact that the gate is a snapshot of the final stage, not a
time average. Mass is CREATED, systematically -- a secular drift, not a random walk.

**Cube vertices ARE worse (item 3, answered).** Per shared face:

| | 16 | 32 | 64 |
|---|---|---|---|
| away from a vertex | -5.91e-8 | -5.85e-9 | -2.30e-10 |
| within 2.5 cells of one | -3.46e-7 | -1.39e-8 | -4.44e-10 |
| ratio | **5.85** | 2.37 | 1.93 |

They dominate the total source at coarse resolution (78% of it at nx2=16 from 37% of the
faces) but converge faster and settle at ~2x per face. Decomposition dependence is small:
4 blocks/panel at nx2=32 reproduces the vertex sum to 6 digits and the total to 2.7%.

## THE FIX -- 806a9e6f, `src/bvals/flux_seam_cc.cpp` (new, ~330 lines)

Each panel packs the OUTWARD flux through its copy of the seam face and both replace
theirs with `F_out <- 0.5*(F_out - F_out_neighbour)`. Both sides evaluate the same
expression, so the result is exactly antisymmetric with no ordering or single-owner rule.
Wired into the EXISTING `SendFlux`/`RecvFlux`/`InitBoundaryValuesAndUtils` task slots,
which were no-ops without `multilevel`. **Reusing the `flux` buffers is safe because
refinement and the cubed sphere are mutually exclusive by a startup FATAL.**

* **NO along-seam resample**, unlike the halo: ghost cell CENTRES do not coincide across
  a seam, but the seam FACES do, one-to-one. A signed index map is exact.
* **Momentum needs `TransformMomentum`, the scalars do not.** Mass, energy and the RADIAL
  momentum flux are scalars (rhat is common to both charts).
* **`iflux_same[0]` for CC had to be added in `buffs_cc.cpp`** -- it existed in the struct
  but only `buffs_fc.cpp` ever filled it.

**THE BUG THAT COST THE MOST: the buffer numbering is not what it looks like.**
x1 faces 0-7, x2 faces 8-15, then the **x1x2 EDGES are 16-23** and x3 faces do not start
until **24**. A `n < 24 is a face` test silently dropped every x3 seam and made things
WORSE. `flux_correct_cc.cpp`'s own `((n<16) || ((n>=24) && (n<32)))` says so if read
carefully. The fix classifies a buffer from its own stored INDEX RANGE (one cell thick
along its normal, full range along the other two), never from `n`.

## Measured, A/B against the same build with the exchange disabled

Hydro rigid rotation, joint refinement, **reflecting radial BCs** (exactly zero physical
flux, so the whole drift is scheme error):

| | nx2=16 | nx2=64 |
|---|---|---|
| mass, no fix | 2.92e-4 | 6.54e-6 |
| mass, with fix | **-2.6e-12** | **-9.7e-12** |
| energy, no fix | 1.35e-3 | 2.88e-5 |
| energy, with fix | **-3.8e-12** | **1.2e-11** |

Exactly conservative, at round-off; a gain of 1e5 to 1e8. The seam's spurious mass source
AND its spurious torque both fall to round-off, cube-vertex bucket exactly 0.
**Angular momentum is not exact and cannot be** -- the momentum equations carry geometric
source terms, so Lz is not a telescoping invariant -- but its convergence improves from
~2.5th to ~3rd order (ratios 5.4-6.5 -> 8.1-9.7). With `user` radial BCs the mass gain is
only 2-3x, because what remains there is genuine boundary flux, not a conservation defect.

Regression: MHD `L1(B)` unchanged within 2% at nx2=16/32/64. MPI 1/2/3 ranks give
bit-identical dt, global sums agreeing to 1e-13, no hang.

## A TRAP I FELL INTO -- read before touching the momentum

`u0(IM1..IM3)` holds the **COVARIANT** momentum `rho*v_i`; `w0(IVX..IVZ)` holds the
**CONTRAVARIANT** velocity `v^i`. This is stated at `coordinates.cpp:609`.
`cubed_sphere::TransformMomentum` is exact ONLY for covariant input -- fed contravariant
components it is wrong by ~46% of |V|, an error that VANISHES on the seam midline and is
worst at the panel corners. I built the first version of the torque diagnostic treating
`u0(IM*)` as contravariant, and it manufactured a fake "Lz plateau" that looked exactly
like a bug in the fix. Raise with g^{-1}, g = [[1,c],[c,1]], c = e1.e2, before using them
as coefficients of e1hat/e2hat. A 40-line standalone test against `cubed_sphere.hpp`
settled it in one run -- do that rather than reason about it.

## Still open

* The gates are **RANK-LOCAL**: under MPI a seam whose two copies are on different ranks
  is invisible, and a rank can print an exact zero because it holds no complete pair.
  Gate the seam on one rank; check MPI via the `.hst` global sums.
* `history.cpp`'s **momentum columns** are still local-basis components on curvilinear
  grids (documented, not fixed) -- summing them is meaningful only on a Cartesian grid.
* `bvals_cc.cpp` calls `pmy_pack->pmesh->GetPanelBoundary(...)` and reads
  `pmy_pack->pmesh->...` INSIDE a device lambda. My new code copies the same pattern for
  consistency, but it dereferences a host pointer on device and looks unsafe for a GPU
  cubed-sphere run. Pre-existing; worth checking before anyone runs cs on a GPU.

Harness: session 98b1b27f scratchpad, `h/` (`base.athinput`, `t0_*`/`tE_*`/`null32`,
`cons.py`, `an.py`); binary is the 9cc449f4 scratchpad's `build/`.
