# Scope: equatorial-band static refinement for the cubed-sphere dhj hydro run (hyd4)

2026-09-23. Scoping only: no code was changed, nothing was built or run. Target: the hyd4 run
(`inputs/production/deep_hot_jupiter_cs_hyd4.athinput`, 128 x 6 x 32^2, blocks 128 x 16 x 16)
with a band around the equator at ~1024 cells around the equator (256 per panel edge, 3 levels
above 32) and the rest of the sphere unchanged. Cost numbers come from the model in
`docs/dev/dhj_resolution_plan.md` (sections 2-4) and the block-counting script described in
section 5. Line numbers are at HEAD dccdf506.

## 0. Verdict

**Feasible in principle, and not worth building for this goal.** The mesh side is roughly
half there. Hydro SMR on the cubed sphere runs and is conservative, but it has only ever been
tested with RADIAL level boundaries. Three facts decide against it:

1. **The octree also doubles the radial resolution at each level** (section 2). A 3-level band
   has 1024 radial cells, and the radial CFL, which sets dt today, then drops dt 8x for the
   whole mesh. AthenaK has one global dt (`mesh.cpp:1005`, an MPI_MIN over all blocks).
2. **The band contains the narrowest horizontal cells.** The equatorial seam midpoints
   (xi = +-45 deg, eta = 0) are 0.724 r dxi wide. So even without point 1, the horizontal dt
   in the band equals that of a uniform 256^2 grid.
3. **Refined regions are quantised by the block grid, need 2:1 buffers, and must stop short of
   the polar seams** (section 3). The refined region comes out wider than asked.

Result (section 5): every band layout that reaches 256 per panel costs **137x to 1417x** the
current work per rotation. A uniform 256^2 grid costs **164x** and needs no new code except a
state remap. The only band layout that is cheaper is root 128^2 plus one level to
|eta| < 22.5 deg, at 137x: 17 % cheaper, for several weeks of new code (section 6).
**Recommendation: do not build equatorial SMR. Run a uniform grid started from an interpolated
developed state** (section 7).

## 1. Geometry

- Panels 0, 1, 2, 4 are the equatorial panels (normals +x, +y, -x, -y) and panels 3 and 5 are
  the poles. On every equatorial panel b = +z, so **eta = x3 is the latitude-like coordinate**
  (`src/coordinates/cubed_sphere.hpp:44-77`). The angles are equiangular,
  xi = (pi/4) x2 and eta = (pi/4) x3 (`src/coordinates/coordinates.cpp:614-619`).
- Latitude on an equatorial panel is lat = atan(tan(eta)/sqrt(1 + tan^2 xi)), which follows
  from `PanelToCart` (`cubed_sphere.hpp:85`). A band |eta| < E is therefore NOT a latitude
  circle. Its edge is at latitude E at the panel centre and atan(tan E/sqrt 2) at the seams:
  - E = 11.25 deg: 11.25 at the centre, 8.0 at the seams;
  - E = 22.5 deg: 22.5 and 16.3;
  - E = 33.75 deg: 33.75 and 25.3.
- The polar seams are at eta = +-45 deg, i.e. latitude 45 deg at the panel centre and 35.26 deg
  at the cube vertices. **A band touches only the 4 equatorial panels and no cube vertex as
  long as every refined level stays inside |eta| < 45 deg.** Given the 2:1 buffers (section 3),
  the finest level must stay inside about |eta| < 31 deg, which is about +-23 deg of latitude at
  the seams.
- The equatorial seams (xi = +-45 deg) cross the band. Along them the eta axes of neighbouring
  panels agree, so a band with the same eta extent on all four panels meets itself at the same
  level across those seams.
- The mesh requires nx2 = nx3 (`src/mesh/mesh.cpp:154`) and the same number of root blocks
  along x2 and x3 (`src/mesh/build_tree.cpp:97-106`).

## 2. What refinement means in AthenaK

- **The octree is isotropic and cannot be changed.** In 3-D every refined node gets 8 children
  (`src/mesh/meshblock_tree.cpp:73-75`, `Refine` at 157-166). No anisotropic or tangential-only
  option exists. A level-L band therefore has nx1 = 128 x 2^L: 256, 512 and 1024 for L = 1, 2, 3.
  - The only escape is a radially coarser root, e.g. nx1 = 64 with the band at L1 = 128. That
    would halve the radial resolution everywhere outside the band, where the requirement is
    already unmet (median 5.0 cells per H_p, `dhj_resolution_plan.md` section 1.1).
- **Levels from 32 per panel:** 3 levels reach 256 (32 -> 64 -> 128 -> 256). A root of 64^2
  needs 2 levels, and a root of 128^2 needs 1.
- **dt:** the binding limit today is radial (dt1 = 20.4 s at r/a_p 1.257; horizontal 62.6 s at
  the 0.724 vertex width; `dhj_resolution_plan.md` section 2). Assume, as that plan does, that
  the binding downdraft is found in the band too. Then the model is
  dt = min(20.1/2^L, 62.6 x 32/n_panel), which gives 10.1, 5.0 and 2.5 s for L = 1, 2, 3.
- **How to define the band:**
  - `<refined_region>` takes x1/x2/x3 ranges and a level (`build_tree.cpp:160-176`). For
    x1 = full, x2 = [-1, 1], x3 = [-E/45, E/45], it adds the node **on every panel**
    (`build_tree.cpp:305`, loop over `npanels`). So it also refines a strip through both poles.
    That strip meets the equatorial panels' polar seams at level 0 and is refused
    (section 3.1). **A panel selector is missing** (work item W1).
  - `user_ref_func` / AMR is not an option: adaptive refinement is refused on the cubed sphere
    (`build_tree.cpp:122-133`).

## 3. Cubed-sphere SMR status

### 3.1 What works (verified in earlier sessions, memory `cubed-sphere-smr`)

- **Hydro SMR with radial level boundaries.** Test input: `inputs/tests/cubed_sphere_smr.athinput`.
  - Mass drift -9.6e-14 and energy drift 1.0e-13 on rigid rotation.
  - Runs at 1/2/3/6 MPI ranks and on 1 and 2 MI300A.
  - Commits f14e075d and 608d43e1.
- **Area-weighted flux correction on all three face directions**
  (`src/bvals/flux_correct_cc.cpp:47-59` and 290-344), and **volume-weighted restriction**
  (`src/mesh/mesh_refinement.cpp:1191-1197`).
- **Exact solid angle.** It makes the geometry additive across levels: one static cycle gives a
  spurious radial momentum of -4.0e-18 (151da303).
- **Cross-panel neighbour search at depth** (`meshblock_tree.cpp:605-618`, via `FindMeshBlock`),
  and **mixed-level cross-panel edge/corner exchange with the seam transform**
  (`bvals_cc.cpp:152` `cs_coar`, 316-319 and 422-426).
- **Cube-vertex corner fills on the coarse array** (180a9b3e). A band inside |eta| < 45 deg
  never reaches a vertex, so this is not exercised here.
- **Restart tree rebuild per panel** (`build_tree.cpp:544-576`, panel id in LogicalLocation),
  with the same seam check after a restart (`build_tree.cpp:610`). No gate for a refined
  cubed-sphere restart is recorded in the notes.
- **Stretched radius under refinement.** The poly stretch maps face positions globally
  (`src/coordinates/grid_stretch.hpp:55-66`, `coordinates.cpp:534`), so fine faces nest in
  coarse ones.

### 3.2 What is missing or untested for an equatorial band

| gap | evidence |
|---|---|
| **A level boundary lying on a seam is refused** ("Prolongation and restriction have no seam path") | `build_tree.cpp:39-80`; `flux_seam_cc.cpp:47-53` shares buffers on the assumption that it never happens |
| **Level boundaries crossing the panel TANGENTIALLY have never been run.** Every cubed-sphere SMR input refines x2 = x3 = [-1, 1] | `inputs/tests/cubed_sphere{,_mhd,_resist,_resist_narrow}_smr.athinput` |
| 2:1 balance is enforced **inside a panel only**. Neighbours across a panel seam are skipped because `panel` is not periodic | `meshblock_tree.cpp:174-245` |
| The doubly-ghost x2x3 edge/corner seam buffers (slots 40-55) are only first order at a level boundary. At a band edge on an equatorial seam these are exactly the mixed-level diagonals | memory `cs-crosslevel-seam-halo-first-order` (measured on the FC path; the CC path is unmeasured) |
| Prolongation is limited-linear in index space. It is not well-balanced for a stratified state (below) | `src/bvals/prolongation.cpp` ProlongateCC; `hydro_tasks.cpp:723-735` |

**Quantisation and the seam constraint.** A level-L region is made of sibling pairs, and its
level-(L-1) parents need level-(L-1) neighbours (`Refine`, `meshblock_tree.cpp:157-245`). Each
level therefore grows by one parent block, rounded up to the grandparent size. Let N0 be the
number of root blocks per panel edge. With N0 = 2 (today: 16-cell blocks, 45 deg), refining
anything makes the whole root block L1 up to eta = 45 deg. That meets polar level 0 on a seam
and is refused. Legal layouts need smaller root blocks:

| root / block tangential | levels | finest band |eta| < | L1 reaches | legal |
|---|---|---|---|---|
| 32^2 / 16 (N0 = 2) | any | - | 45 | **no** |
| 32^2 / 8 (N0 = 4) | 3 | 5.6 | 22.5 | yes (narrowest possible) |
| 32^2 / 8 (N0 = 4) | 3 | 11.25 | 45 | **no** |
| 32^2 / 4 (N0 = 8) | 3 | 11.25 or 22.5 (max 25.3) | 22.5 / 33.75 | yes, with 128 x 4 x 4 blocks |
| 32^2 / 4 (N0 = 8) | 2 | 33.75 | 45 | **no** |

This comes from a 1-D model of the eta direction that reproduces `Refine`'s 2:1 rule (section 5
script). A +-30 deg band at 3 levels is not legal with blocks of 4 or more cells. Allowing it
means lifting the seam refusal (W5).

## 4. Physics modules under SMR

| module | status under a refined band | evidence |
|---|---|---|
| **ck two-stream RT (spherical tm sweep)** | **Silently wrong.** Every block sweeps its own column [is, ie] and treats ie+1 as the top of the atmosphere (B at the top, stellar flux in). A band block holds 1/2^L of the column. There is no MPI and no gather in the sweeps. `two_stream_column_partition.hpp` partitions a column over a thread TEAM, not over blocks. No guard exists: `multilevel` never appears in `two_stream_rt.hpp` or the pgen | `two_stream_rt.hpp:2565, 2594-2622` (top slot); `two_stream_column_partition.hpp:8-30` |
| **pseudo-spherical beam** | Same column assumption: it deposits along the block's own column from ie+1 down to icut | `two_stream_rt.hpp:2407, 2777-2793` |
| **radiative conduction, rad_implicit_x1** | **Guard is insufficient.** It checks `mb nx1 == mesh nx1` (true at the root) but not that a refined block covers the whole radius. The tridiagonal is column-local, so a band block would solve half a column with a false boundary | `src/diffusion/conduction.cpp:259-270` |
| rad_implicit_ang | refused under multilevel (off in hyd4) | `conduction.cpp:339-344` |
| **WB polytropic (wb_x1)** | Reconstruction is cell-local, and the potential comes from each block's own x1v, so it runs. **Risk:** ghost cells at a TANGENTIAL level boundary are prolongated linearly. For an exponential profile the error is about (dr_c/H)^2/32, i.e. 0.4 % at the H_p minimum (2.9 cells/H_p). That sits against a physical day-night gradient of a few % per cell, so band edges could drive spurious flows | `hydro_fluxes.cpp:173-230`; `dhj_resolution_plan.md` 1.1 |
| floors, EOS table, C2P | cell-local: fine | - |
| user radial BCs | applied before `Prolongate`; `FillCoarseInBndryCC` restricts the ghosts into the coarse array | `hydro_tasks.cpp:705-735` |
| initial condition | pointwise from an r-table (`zarr_init`/`logparr_init`): fine | `deep_hot_jupiter_rt.cpp:1025-1052` |
| dt diagnostic | per cell: fine; the global dt is MPI_MIN | `mesh.cpp:1005, 1043` |
| load balance | cost 1.0 per block; fine for equal-size blocks | `build_tree.cpp:360` |
| restart | tree rebuilt per panel. Restarting a UNIFORM run into a refined tree is impossible, because the tree comes from the file. Needs a remap (W7) | `build_tree.cpp:544-576` |
| **bin output / analysis** | The level is written (`binary.cpp:190`) but the **panel is not**. `dhjcs.py:53` assigns the panel as `m // (nmb/6)`, which is wrong once panels hold different block counts. `level_interp` assumes one block holds the whole column | `src/outputs/binary.cpp:135-190`; `docs/handover/scripts/dhjcs.py:49-53, 81-94` |

## 5. Cost

Model, all from `dhj_resolution_plan.md`:
- dt0 = 20.1 s (radial);
- 8.91e5 cycles to rot 50 at dt0;
- per-GPU rate 6.56e6 zone-cycles/s on one node (arm C x 0.968), divided by 1.5 when going
  off-node (measured 1.56x at 6 blocks per GPU; **unmeasured at these loads**);
- memory about 12M zones per APU;
- 2 APUs per node.

For comparison, hyd4 in production runs at 15.5 cycles/s on 2 GPUs (100 cycles per 6.46 s,
`bench/cs_hyd4_prod/log.out.11943719`, cycles 21500-21700). That is 0.27 h per rotation at
dt0, against 0.25 h from the model.

Other assumptions:
- The band spans the full radius.
- The polar panels stay at the root level.
- W = zones/zones0 x dt0/dt, i.e. the work per rotation relative to today.
- The zone rate is taken as independent of column length. Arm E says longer columns cost
  1.24x per zone at 1.5x radial, so the band rows are **optimistic**.

Zone counts come from exact block counting, not area fractions.

| grid | zones (M) | dt (s) | W | h/rot, 1 node | to rot 50, 1 node (h) | h/rot, 4 nodes | to rot 50, 4 nodes (h) | h/rot, 16 nodes | to rot 50, 16 nodes (h) |
|---|---|---|---|---|---|---|---|---|---|
| uniform 32^2 (today) | 0.8 | 20.1 | 1 | 0.25 | 15 | 0.09 | 6 | 0.02 | 1 |
| uniform 64^2 | 3.1 | 20.1 | 4 | 1.01 | 59 | 0.38 | 22 | 0.09 | 6 |
| uniform 128^2 (512 around the equator) | 12.6 | 15.7 | 21 | 5.2 | 305 | 1.9 | 114 | 0.48 | 28 |
| **uniform 256^2 (1024 around the equator)** | 50.3 | 7.8 | **164** | mem | - | 15.5 | 909 | **3.9** | **227** |
| 32^2, L1 band 11.25 (blocks 4^2) | 1.7 | 10.1 | 4 | 1.1 | 64 | 0.41 | 24 | 0.10 | 6 |
| 32^2, L1 band 33.75 (blocks 4^2) | 3.5 | 10.1 | 9 | 2.3 | 134 | 0.85 | 50 | 0.21 | 12 |
| 32^2, L2 band 11.25 (L1 to 22.5) | 10.0 | 5.0 | 51 | 12.8 | 752 | 4.8 | 280 | 1.2 | 70 |
| 32^2, L2 band 22.5 (blocks 4^2) | 18.2 | 5.0 | 93 | 23 | 1375 | 8.7 | 512 | 2.2 | 128 |
| **32^2, L3 band 5.6** (blocks 8^2) | 39.3 | 2.5 | **400** | mem | - | 38 | 2212 | 9.4 | 553 |
| **32^2, L3 band 11.25** (blocks 4^2) | 72.4 | 2.5 | **736** | mem | - | 69 | 4070 | 17 | 1017 |
| **32^2, L3 band 22.5** (blocks 4^2) | 139 | 2.5 | **1417** | mem | - | mem | - | 33 | 1959 |
| 64^2, L2 band 11.25 | 39.8 | 5.0 | 203 | mem | - | 19 | 1121 | 4.8 | 280 |
| 64^2, L2 band 22.5 (blocks 8^2) | 72.9 | 5.0 | 371 | mem | - | 35 | 2050 | 8.7 | 512 |
| **128^2, L1 band 22.5** | 41.9 | 7.8 | **137** | mem | - | 12.9 | 758 | 3.2 | 189 |
| 128^2, L1 band 33.75 | 56.6 | 7.8 | 185 | mem | - | 17 | 1023 | 4.4 | 256 |

"band X" means level L covers |eta| < X deg (latitude X at the panel centre). Readings:
- **Your ~170x for uniform 256^2 checks out:** 164x (64x zones, dt 2.57x).
- **A 3-level band costs more than uniform 256^2 even at +-5.6 deg,** because 8x radial means
  8x dt and 8x zones per unit area in the band.
- The +-10/20/30 deg bands at 3 levels are 736x, 1417x, and not legal.
- At 2 levels from 32^2 (512 per equator) the bands cost 51-93x, against 21x for uniform 128^2.
  **SMR loses at every target resolution.**

Script: the eta-direction tree model and the table are in this session's scratchpad
(`band.py`, `tab.py`). Rerun them with a different dt model if a measured band dt arrives.

## 6. Work items, if the band is built anyway (cheapest legal target: 128^2 root + L1 band)

| # | item | effort | risk |
|---|---|---|---|
| W1 | `<refined_region>` panel selector (a `panels` list in the 3-D loop, `build_tree.cpp:305`) | 0.5 d | low |
| W2 | guards: refuse ck RT, the beam and rad_implicit_x1 on radially split blocks, i.e. a refined block whose x1 range is not the mesh range | 0.5 d | low. **Worth doing anyway:** both are silently wrong today under any radial split |
| W3 | radial column gather for the ck sweep, the beam and the rad_implicit_x1 tridiagonal across 2^L blocks and ranks. Template: the rad_m1 gathered line solve, uniform decompositions only (dccdf506; `rad_m1_implicit.cpp:534` excludes multilevel) | 1.5-3 wk | **high**: 72 % of the cycle is the column sweep (`dhj_resolution_plan.md` 3); the extra communication is on the critical path |
| W4 | tangential in-panel level boundaries: a static stratified gate (zero acceleration, as in 151da303), conservation, and a well-balanced prolongation (log-space or HSE-anchored) | 1-2 wk | medium-high |
| W5 | level boundary ON a seam (prolongation and restriction with the seam transform). Needed only for bands wider than the table above or for coarse polar panels next to L1 | 2-4 wk | high (the cross-level seam halo took several sessions, memory `cs-crosslevel-seam-halo-first-order`) |
| W6 | panel id in bin output; `dhjcs.py` per-panel mapping, multi-level columns and regridding | 2-3 d | low |
| W7 | remap a developed state onto the new grid (the pgen reads a coarse restart or bin and interpolates within each panel; the panel frames are resolution-independent) | 3-5 d | medium |
| W8 | GPU efficiency of 4^2 or 8^2 tangential blocks (16^2 vs 32^2 already costs 8.5 %, arms B/C); 4^2 is unmeasured | timing only | medium |

W1 + W2 + W3 + W4 + W6 + W7 is about 4-7 weeks, and it buys at most 17 % (137x vs 164x).

## 7. Alternatives

- **Latitude-stretched horizontal grid: not available, and structurally limited.**
  - The equiangular map is hard-wired: `M_PI/4 * CellCenterX` appears at 14 sites in
    coordinates.cpp, 10 in bvals_cc.cpp, 8 in bvals_fc.cpp and 4 in flux_seam_cc.cpp.
  - The seam resample assumes a uniform angle step, `dang` (`bvals_cc.cpp:360-372`), and that
    both charts share the seam-normal coordinate exactly (`bvals_cc.cpp:331-337`).
  - Both tangential axes of a polar panel are along-seam axes of equatorial seams. So longitude
    resolution cannot be local: 1024 around the equator needs nx2 = nx3 = 256 on every panel
    (`mesh.cpp:154`). Only eta on the equatorial panels could be stretched. That breaks the
    seam-normal match at the polar seams, and it gives latitude resolution only.
- **Uniform high resolution from an interpolated developed state: RECOMMENDED.**
  - Every module works unchanged on a uniform grid: whole columns per block, 128 x 16 x 16
    blocks.
  - The only new code is W7.
  - 256^2: about 3.9 h per rotation on 16 nodes (32 APUs, 1.6M zones each; the arm C load).
    5 rotations take about 20 h, i.e. about 620 GPU-h.
  - 128^2 fits on 1 node (6.3M zones per APU): 5.2 h per rotation there, or about 1.9 h per
    rotation on 4 nodes.
  - Caveat: the dt model assumes that finer grids do not develop faster downdrafts. The plan
    flags this as unmeasured.

## 8. Recommended first step

Before writing any code, measure what the model assumes. Take a **multi-node GPU timing of
uniform hyd4 at 128^2 and 256^2** on the apu partition (apudev is capped at 1 node):
- start from scratch, 150 cycles, outputs off, the same binary as bench/cs_hyd4_prod;
- 128^2 on 1, 2 and 4 nodes; 256^2 on 4, 8 and 16 nodes;
- two interleaved repeats.

This replaces the assumed 4.4e6 zone-cycles/s per off-node GPU, which is the largest
uncertainty in every row above, and settles whether 256^2 for a few rotations is 20 h or much
more. Then build W7 (the remap) and W2 (the guards).
