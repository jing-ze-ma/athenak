# Cubed-sphere seam defect in the deep hot Jupiter runs (2026-09-26)

## Summary

The dt collapses of the WASP-121b 1x sponge arms (`nobot`, `notop`) sit on panel seams. They
come from the seam flux reconciliation in `src/bvals/flux_seam_cc.cpp`
(`RecvAndUnpackFluxSeamCC`, rt-integration lines 470-476). That code replaces each seam-face
flux by the plain average of the two panels' outward fluxes:

    flx.x2f(m,v,k,j,i) = 0.5*flx.x2f(m,v,k,j,i) - 0.5*sgn*recv;     // and x3f

Each panel computes its seam flux from ghosts that are resampled along the seam out of the
other panel's cells (`bvals_cc.cpp`, the 3-point clamped resample around line 441). If a
seam cell is much emptier than its along-seam neighbours, the other panel's ghost copy of it
is a blend with those denser neighbours. That panel's Riemann solve then moves mass
accordingly. The average applies half of that flux to the real cell, and the result is
**not positivity-preserving**.

What was measured in `nobot`:
- The cell holds rho = 1e-14.
- Its own solve says it gains 6e-8 g/cm^2/s.
- The neighbour panel's solve says it loses 3.3e-7.
- The average removes 1.3e-7, which is 30 % of the cell's mass per step.
- Within about 10 cycles this grows to 38x the cell's mass per step. That is the dt collapse.
- The legacy density floor (which keeps momentum) then turns the drained cell into the
  v ~ 1e11 cm/s, T ~ 1e11 K runaway.

The healthy state shows no seam anomaly: the defect is a local instability, not a bias. It
only switches on when a seam cell becomes much emptier than its along-seam neighbours
(night side, 1e-7..1e-6 bar).

The fix is on branch `cs-seam-0926`: `<mesh>/cs_seam_flux = average | upwind | positive`.
The default is `average`, which is bitwise the old behaviour. The fix decides per face, from
the mass flux of both panels, which panel is the donor, and uses that panel's own flux for
every variable. `positive` also zeroes the mass flux (by convex weights) when both panels
claim inflow. Both options stay exactly conservative.

## Setup

- WASP-121b 1x, C32 cubed sphere (6 x 2 x 2 blocks of 16x16), nx1 76, hydro (general EOS
  table), PLM + HLLC, rk2, cfl 0.3, ck RT (`deep_hot_jupiter_rt`).
- Runs: `/viper/ptmp2/jinma/sparc_w121x1/{base,nobot,notop}`.
- Production binary: sparc_0925/athena.gpu, commit f637b3be, md5 e3a8442e.
- Crash reruns:
  - `nobot`: from rst 12 (rotation 5.995, cycle 46381). FATAL at cycle 47063.
  - `notop`: from rst 11 (rotation 5.5, cycle 78496). FATAL at cycle 83539.
  - Both reproduce bit for bit.
  - Evidence of the earlier analysis: `/viper/ptmp2/jinma/sponge_crash_0926/RESUME.md`.
- Work dir: `/viper/ptmp2/jinma/seam_0926`. Its RESUME.md lists every file and job.

## Q1: is the seam anomalous in the healthy base arm?

**No, except within about 3 cells of a night-side cube vertex.**

Data: base rst 20..38 (rotations 10-19, 10 dumps), from `seamstat.py`.

Cells were classed by position on the panel face:

| class | cells |
| --- | --- |
| `vert` | cube-vertex cell |
| `edgeV` | face edge, within 3 cells of a vertex |
| `edgeM` | face edge, mid-seam |
| `edge1` | one cell inside a face edge |
| `blk` | internal block boundary, a copy halo (a control) |
| `int` | interior |

The classes were compared within bins of local pressure (decades) and day/night. The
metrics were:
- grid-scale roughness of log T and log rho normal to the seam, with the seam ghost as the
  across-seam neighbour;
- the in-panel one-sided curvature;
- Mach number, with a KE that accounts for the cubed-sphere geometry;
- horizontal velocity divergence x r dtheta / c_s;
- the fraction of cells at the T floor;
- the fraction of cells with M > 5.

Mid-seam versus interior, night side (medians):

| p [bar] | aT edgeM / int | Mach edgeM / int | T-floor % edgeM / int |
| --- | --- | --- | --- |
| 1e-7..1e-6 | 0.0069 / 0.0062 | 1.11 / 1.14 | 0.011 / 0.94 |
| 1e-6..1e-5 | 0.0022 / 0.0022 | 3.39 / 3.13 | 0 / 0.073 |
| 1e-5..1e-4 | 0.0011 / 0.0017 | 3.53 / 3.03 | 0 / 0.010 |

The p99 values are also comparable (0.198 vs 0.305 at 1e-7..1e-6 bar). **Mid-seam cells
are no rougher and hit the T floor less often than interior cells.**

Within 20 deg of the 4 night-side cube vertices (lon +-135, lat +-35), the vertex-adjacent
cells (`vert`, `edgeV`) are slightly anomalous. At 1e-6..1e-4 bar:
- median M is 4.7-5.2 against 4.1-4.4 in the interior of the same window;
- the fraction with M > 5 is 39-73 % against 20-30 %;
- the median log-T roughness is 1.6-4x the interior value.

Both crash sites lie 2-9 cells from such a vertex.

Seam-ghost consistency (`gcheck.py`, base rst 20..36, night, 1e-7..1e-3 bar, mid-seam) was
measured as |log q_ghost - the linear extrapolation of the two active cells|. Against the
same statistic with a real third cell:
- p: median 0.047 vs 0.030 at M 5-8;
- T: 0.020 vs 0.011;
- p99 is lower at the seam.

**The seam ghosts are as consistent as real cells.** The ghost exchange has no M^2-growing
error in the healthy state.

Seam flux mismatch in a healthy run (base rst 30 run for 0.4 rotations, `CS_SEAMDIAG`):
- the per-step maximum over all seam faces of phiM = 0.5 |F_own + F_nb| dt / (dx rho) was
  0.026 median and 0.042 maximum;
- no face exceeded 0.3.

## Q2: which operator?

Six arms were run per crash case, one GPU each on apudev, from the crash restarts:
- `diag` is the production commit with the read-only diagnostic.
- `cfl29` and `cfl31` perturb the trajectory, as controls.
- `CS_*` are environment switches built into the snapshot `seam_0926/src`.

Job 11982302 (`nobot`) and job 11982303 (`notop`):

| arm | nobot | notop |
| --- | --- | --- |
| diag (= production, bitwise) | FATAL 47063 | FATAL 83539, 20 collapse warnings |
| cfl 0.29 | FATAL 47248 | survives, 4 collapse warnings |
| cfl 0.31 | FATAL 47243 | FATAL 81855, 8 collapse warnings |
| RT frozen (`ck_impl_every=1e9`) | survives, collapse warning at 46610 | **FATAL 78720 (earlier)** |
| no seam flux averaging (`CS_FLUXSEAM_OFF=1`) | **clean** (0 warnings) | survives, 2 warnings |
| no along-seam resample (`CS_NORESAMP_CC=1`) | **clean** | **clean** |

Jobs 11982326 and 11982366 tested the prototypes on the same arms:

| arm | nobot | notop |
| --- | --- | --- |
| upwind-owner (proto `=2`) | **clean** | **clean** |
| upwind-owner + cfl 0.29 | **clean** | **clean** |
| positive (proto `=3`) | **clean** | **clean** |
| positive + cfl 0.31 | **clean** | **clean** |

"Clean" means it reached tlim with zero dt-collapse warnings and dt stayed at 6-14 s. The
crash arms had dt of 1e-4 s.

The diagnostic (CS_SEAMDIAG=0.3) prints every seam face whose phiM exceeds 0.3.
- **nobot**: every one of the 1968 printed faces is at the crash site (m=23, k=17, j=8..10,
  i=50..58, an x3 seam face). The per-step maximum phiM was 0.03-0.08 for 280 cycles, near
  base values. It then rose to 0.15, crossed 0.3 at cycle 46667, and reached **38.7 at the
  first dt collapse (46680)**. It stayed at 5-19 until the FATAL.
- **notop**: phiM reached 170-300. The faces are at the crash site (m=0, k=2, j=3..8) and
  at the other collapse sites (m=17/19, j=17).
- A typical face: `own=-6.05e-08` (this panel: inflow), `nb=+3.27e-07` (neighbour panel:
  outflow from this cell), rho = 1.1e-14.

Reading of the evidence:
- **(b) the seam flux reconciliation is the operator that fails.** The fault is its
  positivity, not its conservation.
- **(a) the resample supplies the disagreement.** Remove either one and the collapses go.
  Keeping both, with a donor-based choice in place of the mean, also removes them.
- **(c) the ck RT is excluded.** It is column-by-column: `two_stream_column_ck.hpp` never
  reads a j/k neighbour, and `rad_angular = false`. Freezing it does not remove the
  collapse; notop crashes sooner.
- **(d) the floors and EOS are not the origin.** The C2P and floor path is the same for
  every cell. The legacy dfloor (which keeps momentum) only turns the drained cell into the
  runaway. The floor switches `dfloor_keep_*` and `vceil` (tested separately in
  sponge_crash_0926) also stop the FATAL, but they treat the symptom.

### Third case: the base arm crashed too (FATAL at cycle 221708, rotation 27.04)

Every dt-setting cell of the base crash is on one seam: panel 0, face j = 31. The cells are
(gid, k, j, i) = (1,10,17,60), (1,12,17,59), (3,11,17,59) and (3,9,17,58), all with j = 17
in the ghost-padded index. (1,2,17,49) is the cube-vertex cell at face (j 31, k 0).

Rerun from base rst 53 (rotation 26.5), job 11982463:
- The production-equivalent binary with the diagnostic reproduces the FATAL at 221708
  bitwise.
- Every seam face with phiM > 0.3 is on that seam (m = 1, j = 17; k = 10, 17, 9, 7).
- With the branch binary, `average`, `upwind` and `positive` all survive with no collapse.
  The branch is rt-integration HEAD, not the production commit, so its trajectory differs
  and this case does not discriminate on the branch.

### Solver speed settings (ck cadence and tolerance) are not the cause

Job 11982468 reran the crashes with the production binary (e3a8442e) and ck solved every step
(`ck_impl_every=1`; the input has 4) with `ck_impl_tol=1e-9` (the input has 1e-8):
- **nobot still crashes at the same cells.** The first collapse is at 46668 (production
  46680) and the FATAL at 47753 (production 47063). The cells are (23,17,8..11,50..54).
- **notop** has 4 collapse warnings starting at the same cycle 78558, at the known seam
  cells (17/19, j 17), but it reaches tlim. The notop FATAL is sensitive to perturbation:
  the cfl 0.29 control also survived, with 4 warnings. So this arm behaves like a
  perturbation control.

The ck cadence and tolerance are cleared.

### Branch binary on the crash reruns

Job 11982451 ran rt-integration cc529150 + the key, with the key given via `-i key.athinput`
on restart:

| key | nobot | notop |
| --- | --- | --- |
| `average` | **FATAL 46764** (3 warnings) | **FATAL 80495** (7 warnings) |
| `positive` | clean, tlim | clean, tlim |
| `upwind` | clean, tlim | clean, tlim |

### Smooth test

`cs_test` iprob 3, rigid rotation about an edge axis, p0 = 0.01, omega = 0.65 (M up to
about 1), nx1 = 8, tlim = 1. L1 of v / p:

| nx | average | upwind | positive |
| --- | --- | --- | --- |
| 16 | 2.4715e-3 / 3.7149e-3 | 2.5185e-3 / 3.6217e-3 | same as upwind |
| 32 | 1.2543e-3 / 1.4293e-3 | 1.2578e-3 / 1.4120e-3 | same as upwind |
| 64 | 1.1003e-3 / 1.1391e-3 | 1.1009e-3 / 1.1384e-3 | same as upwind |

The mass drift is identical in all three (-7.503e-05 at n64, from the radial boundary). The
seam adds no drift in any mode. Seam cells do not drift or generate spurious velocity
beyond the interior level. This matches earlier results: the seam halo is second order, and
a blast crossing a seam behaves as it does on a Cartesian grid.

The healthy base (rst 30, 0.4 rotations, default vs positive) ran 124153 vs 124152 cycles
with the same dt (14.05 vs 14.04 s). At the end, mass differs by 3e-8, total E by 5e-7, and
the KE components by 3e-4..2e-3. The KE components themselves drift 1-2 % over the interval.

## Q3: code location and fix

- Defect: `src/bvals/flux_seam_cc.cpp`, `RecvAndUnpackFluxSeamCC`, the average at
  rt-integration lines 473/475 (branch lines 546/548). The comment at the top of the file
  (lines 26-27) describes it.
- Contributing: `src/bvals/bvals_cc.cpp`, the along-seam 3-point resample of the
  cell-centred conserved state (rt-integration lines about 436-450). It makes the two
  panels' seam Riemann problems differ.
- Fix (branch `cs-seam-0926`): `<mesh>/cs_seam_flux` in `src/mesh/mesh.cpp`, member
  `Mesh::cs_seam_flux`, and a separate kernel `SeamFluxRecvUpw` in `flux_seam_cc.cpp`.
  - a = this side's outward mass flux, b = the neighbour's.
  - a > 0, b < 0: this side is the donor, so it keeps its own flux.
  - a < 0, b > 0: take the neighbour's flux.
  - Otherwise the mean.
  - `positive` additionally sets wa = b/(a+b) when a < 0 and b < 0 (mass flux exactly
    zero), and handles the zero-flux edges.
  - The neighbour's weight is 1 - wa, so both panels write the same number and the update
    stays conservative.
  - MHD goes through the same function. Its CT/EMF path is untouched.
- Gates (CPU, 2 MPI ranks, `cs_test`, 200 cycles): `rigidrot`, `blast` and `mhd` with the
  key off are **bitwise identical** to rt-integration cc529150, with and without an explicit
  `cs_seam_flux = average`. With the key on, all three run and differ, as they should.

Recommendation: use `positive` for dhj production and, once longer runs confirm it, make it
the default. `upwind` alone already cured all 8 crash arms. `positive` also covers the one
case (both panels claim inflow) where the mean still drains a cell.

## Status and next steps

See `/viper/ptmp2/jinma/seam_0926/RESUME.md` for the job list. Still open:
1. A longer check (several rotations) of `positive` from a base restart before making it
   the production default. Compare seamstat.py numbers and the T-floor/dfloor counters.
2. Optional: the vertex-adjacent Mach excess (Q1) is a separate, milder effect. With
   cs_vertex_fill_cc = 0 the cell-centred corner block is extrapolated.
