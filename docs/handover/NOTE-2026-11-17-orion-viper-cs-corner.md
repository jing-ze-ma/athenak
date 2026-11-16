# NOTE 2026-11-17 Orion -> viper: part 2 of TASK-2026-11-17 (cubed-sphere corner unit test for the RSG E corner jet)

## Answer
- **No key points at the corner jet in hydro-only.**
  - There were 60 CPU runs: 5 blob types x 3 placements x 4 arms.
  - In 13 of the 15 blob x placement cases, the four arms are **bitwise identical**.
  - The other 2 are the 1000x-contrast blob at the corner and at the seam. There the guard acts, changing rho by at most 1e-4 relative.
- **No corner focusing or jet.**
  - A blob centred on the cube corner sends out the same front as one centred on a panel or on a seam.
  - Front radius agrees within one 0.02 bin. Peak amplitude agrees within about 2 % at N=64 and about 6 % at N=32.
  - The corner value is not systematically larger.
- **Recommendation: no GPU twins of RSG for these hydro keys** (guard 4 vs 0, 3cell on vs off).
  - The jet must come from something this test leaves out. See "What this test cannot see" below.

## Code / setup
- Binary: fork `race3-prod-1118-viper` 0e7f5ca26 (it has all five keys), CPU build, `-DPROBLEM=cs_test`, MPI, gcc 13, `-ffp-contract=off`.
- Pgen addition: a 3-D blob option for cs_test iprob 12, 17 lines.
  - It lives on a local throwaway branch on Orion and is **not pushed**. Patch: /orion/ptmp/jinma/cscorner_1117/blob.patch.
  - Keys `<problem>/blob_w, blob_r, blob_d, blob_p, blob_v, blob_tophat`.
  - The profile is a Gaussian (or a top hat of radius w) in 3-D Cartesian distance from the point `blob_r * (blast_x, blast_y, blast_z)`:
    - rho = 1 + blob_d g;
    - p = 1 + blob_p g;
    - v_r = blob_v g.
- Background: uniform rho = p = 1, gamma 5/3, at rest, no gravity.
  - r in [1, 3]; inner boundary reflect, outer boundary outflow.
  - Blob at r0 = 1.6, w = 0.2 (5 tangential cells at N=64).
- Hydro as in the RSG deck: plm, hllc, `fofc = true`, nghost 3, rk2, cfl 0.3, `cs_seam_flux = shared`, `cs_seam_resample = linear`.
  - am_transport and layout_invariant are on by default.
- Grid:
  - Batch 1: N = 64 per panel, nx1 = 64, 2x2 MeshBlocks per panel (so `cs_seam_blockend = extend` is active), t = 0.30.
  - Batch 2: N = 32, t = 0.10, with `<hydro>/fofc_report = 1`.
- Placements are on the **0/1/5 corner** of the RSG report: C = (1,1,-1), E = seam midpoint (1,0,-1) between panels 0 and 5, P = centre of panel 5 (0,0,-1).
- Arms, all naming every key explicitly:

  | arm | `cs_seam_rho_guard` | `cs_vertex_3cell` | `cs_seam_blockend` |
  |---|---|---|---|
  | prod | 4 | true | extend |
  | g0 | 0 | true | extend |
  | v3off | 4 | false | extend |
  | both | 0 | false | extend |

  - `cs_vertex_mono` and `cs_cross_mono` are `<rad_m1>` keys (rad_m1_implicit.cpp). They do nothing in a hydro-only run, so this test cannot cover them.
- Blob types:

  | id | blob | blob_d | blob_p | blob_v | N | t |
  |---|---|---|---|---|---|---|
  | S | smooth Gaussian | 3 | 3 | 0 | 64 | 0.30 |
  | H | top hat | 30 | 30 | 0 | 64 | 0.30 |
  | V | outward v_r pulse | 0 | 0 | 1.0 (Mach 0.8) | 64 | 0.30 |
  | X | hot top-hat cavity | 0 | 100 | 0 | 32 | 0.10 |
  | D | dense hot top hat | 1000 | 1000 | 0 | 32 | 0.10 |

- Cost: all 36 + 24 runs ran concurrently, 4 ranks each, on the Orion holds 230761 and 230775 (heavily shared). Batch 1 took 24 min wall, batch 2 took 3 min.
- Run directories: /orion/ptmp/jinma/cscorner_1117/runs/{S,H,V,X,D}_{C,P,E}_{prod,g0,v3off,both}. Analysis: ana.py, cmp.py, plot.py there. No restart dumps were written.

## Do the keys act at all? (cmp.py: max |difference| vs prod of rho and v_r over the whole final dump)

| blob | C | E | P |
|---|---|---|---|
| S, H, V, X | bitwise in all 3 arms | bitwise | bitwise |
| D (rho x1000) | g0, both: d rho 0.082 (rho ~ 1e3), d v_r 7e-5; v3off bitwise | g0, both: d rho 0.090, d v_r 2.5e-3; v3off bitwise | bitwise (blob has not reached a seam) |

- **`cs_vertex_3cell` is inert in hydro.**
  - In hydro it changes only the FOFC trial state at a vertex (hydro_fofc.cpp CsVertexFofcDrop).
  - `fofc_report` counted **0 FOFC cell-stages in every run**, including the X and D blasts.
  - In RSG it also switches the M1 vertex VET terms (M1Vtx), a radiation path.
- **`cs_seam_rho_guard` fires only at a jump of more than 4x across the seam resample.**
  - With the linear resample, a 31x top hat (H) never triggers it.
  - A 1000x one (D) does, with a 1e-4 relative effect and no change in the metrics below at 3 digits.

## Front anisotropy (prod arm; the other arms are identical, see above)
Column definitions:
- sf: distance from the blob centre of the peak of max(dp) in 0.02 bins. It is measured in the outward cone (within 20 deg of the blob axis) and in the lateral ring (within 10 deg of perpendicular).
- pk: that peak dp/p0.
- axis v_r and axis drho: max over cells within 0.06 rad of the blob axis, outside r0.
- lat max/med: dp at the lateral front, maximum over median (azimuthal imprint).

| blob | pos | sf out | sf lat | pk out | pk lat | axis v_r | axis drho | lat max/med |
|---|---|---|---|---|---|---|---|---|
| S | C | 0.59 | 0.59 | 0.170 | 0.171 | 0.160 | 0.147 | 1.046 |
| S | E | 0.59 | 0.59 | 0.170 | 0.170 | 0.160 | 0.146 | 1.058 |
| S | P | 0.59 | 0.59 | 0.170 | 0.168 | 0.160 | 0.147 | 1.048 |
| H | C | 0.71 | 0.71 | 0.670 | 0.685 | 1.196 | 2.338 | 1.084 |
| H | E | 0.69 | 0.71 | 0.673 | 0.683 | 1.203 | 2.317 | 1.078 |
| H | P | 0.71 | 0.69 | 0.659 | 0.701 | 1.182 | 2.214 | 1.089 |
| V | C | 0.57 | - | 0.175 | 0.083 | 0.589 | 0.101 | 1.064 |
| V | E | 0.57 | - | 0.175 | 0.082 | 0.584 | 0.101 | 1.108 |
| V | P | 0.57 | - | 0.175 | 0.082 | 0.578 | 0.102 | 1.110 |
| X (N=32) | C | 0.63 | 0.59 | 2.72 | 3.12 | 1.425 | 0.514 | 1.140 |
| X (N=32) | E | 0.63 | 0.63 | 2.93 | 3.14 | 1.456 | 0.579 | 1.169 |
| X (N=32) | P | 0.63 | 0.61 | 2.91 | 2.97 | 1.455 | 0.549 | 1.162 |

- **Corner/centre ratio of axis v_r:**
  - S 1.00, H 1.012, V 1.019 at N=64;
  - X 0.98 at N=32, where the corner is *weaker*.
- **Peak dp, corner/centre:**
  - outward: S 1.00, H 1.017, X 0.94;
  - lateral: H 0.98, X 1.05.
- The front speed is the same to one 0.02 bin (<3 %) at all three placements.
  - In the inward cone at small s the profiles differ more, because of the reflecting wall at r = 1, 0.6 from the blob (same for all placements).
- **Spontaneous corner signal** when the blob is far away (P runs): over all radii, the excess of v_r within 0.06 rad of the corner over the same-radius median is <= 6e-10 (S), and <= 1e-12 for the other blobs. The corner stays quiet to round-off.
- D (rho x1000) has moved only about 2 cells by t = 0.1, so its front metrics measure how the top hat is sampled on the grid, not propagation. It is used only for the guard check above.

Plots:
- docs/handover/cscorner-1117/cs_corner_profiles.png: S/H/V, cone profile and axis v_r, all 12 position x arm curves per panel. They overlap.
- docs/handover/cscorner-1117/cs_corner_hotblob_n32.png: X.

## What this test cannot see (where the jet must come from)
The production RSG has the following, all absent here:
- Radiation:
  - `cs_vertex_mono` and `cs_cross_mono`;
  - the M1 vertex VET terms that `cs_vertex_3cell` also switches.
- Gravity plus well-balancing at the seam (`wb_x1 polytropic`, wb_seam).
- `merge_lat` lateral merging of 8:4:2 at r < 5.7e12.
- General EOS.
- Floor-level density contrasts (1e-20 floor against a 1e-9 envelope). These are what the guard was built for, far above the 1000x tried here.
  - Caltech guard 0 = clean for 1822 d, against viper guard 4 = jet. This test does not rule out a guard role at those contrasts **with** gravity. It only shows that, in plain hydro, the guard neither creates nor amplifies a corner front.

Cheapest next CPU discriminators, if wanted:
1. A stratified, hydrostatic run with noise (cs_test iprob 13 + vnoise), with guard 4 vs 0. Optionally add merge_lat.
2. An M1 corner test with `cs_vertex_mono` / `cs_cross_mono` / 3cell on vs off.

Orion can do either in under 1 h on request.
