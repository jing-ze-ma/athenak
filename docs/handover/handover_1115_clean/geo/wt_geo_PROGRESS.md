# clean-cc-geo progress (cell layout on spherical polar), 11-15 PST

Worktree /viper/ptmp2/jinma/clean_1114/wt_geo, branch clean-cc-geo off clean-cc-bstress 2a68179c7 (local).
Tests: /viper/ptmp2/jinma/clean_1114/geo/gt.py (inputs, exact solutions, scorer), runs in geo/runs.

## 1. Survey: what in the cell layout assumed Cartesian
| part | before | spherical polar (this branch) |
|---|---|---|
| divergence (CleanCellDivergence, CleanCellSourceDivergence) | h idx_d (f+ - f-), uniform dx | h (A+ f+ - A- f-)/V from Coordinates area/volume |
| momentum curvature terms | none | angular-momentum conserving stress fluxes (r faces of y_theta, y_phi weighted r_f/r; theta faces of y_phi sin_f/sin) + well-balanced source -c (D_tt + D_pp)(A_r+ - A_r-)/(2V) E in y_r, -c D_pp (A_t+ - A_t-)/V E in y_theta (cancels a uniform isotropic stress exactly) |
| RC (Jacobian + flux) | 0.75/0.25 idx (uniform 2dx central, dx compact) | CleanCellRhieChow: central g over x_{i+1} - x_{i-1}, compact over the centre distance (dxface: dr, r dtheta between centroids) |
| dissipation weight tau = chi dx | mb_size dx (dtheta on sp!) | physical centre distance |
| cell preconditioner blocks (CleanCellPreconditionerSetup) | s = h idx | per-row area weights + curvature diagonal block; line solves geometry-free |
| boundary faces (CleanCellFaceJacobian/Source) | geometry-free per unit area | unchanged (areas enter through the divergence); beam needs vet_sc -> not on sp |
| positivity check |F| <= cE | basis invariant | unchanged (orthonormal local basis) |
| closure | vet_sc (VetShortCharSweep, Cartesian sweep; refused on sp/cs by SphericalS1Check) | NOT AVAILABLE: sp rays are vet_gd (rad_m1_vetgd.cpp), see below. Test closure clean_closure = sphere |
| gas coupling | a_b = v_b + D_bc v_c, dm from div | local-basis velocity, dm includes the curvature terms (gas feedback not tested on sp) |
FATALs (CleanInit): cubed sphere; sp without clean_layout = cell + clean_closure = sphere; use_polar_boundary;
SMR/AMR; clean_closure = sphere on a Cartesian mesh.

Closure on sp (open): vet_gd writes D = K/J relative to vet_col's uniaxial tensor into tau_ten LAT slots (folded into
slot 0 by VetLatBuild), the half range r+-_a in vgd_hr (M1_VGD_HRC1 = C1, HRJ = J with vgd_thx), and has no formal
solution below vlat_icut.  A CleanClosureFromGd (D full tensor in the local basis, h+- = r+-, C1, J) is needed; vet_col
(column optical depth) has to be checked against the no-tau-in-blend rule.

## 2. Commits (clean-cc-geo, local)
- 5840de770 sp geometry + clean_closure = sphere; 3ea368946 geometry arrays as plain views (build fix);
  51b860d62 comment to 90 columns.  Binary /viper/ptmp2/jinma/builds/bin/athena_none_cpu_3ea36894_geo2.

## 3. Results (CPU, login node; geo/runs/<case>/<tag>_<N>, scored with gt.py)
All runs hesdirk2 with 1 BE step (the first) and 0 stage fallbacks unless noted; clean_tol 1e-14.

(a) free-streaming split pulse, D_rr = 1 (sphere closure R = 0), r in [1, 24], E0 = 0.1 + exp(-((r-10)/0.3)^2), F0 = 0,
c dt/dr = 0.5, T = 4; exact r^2 E = (G(r-ct) + G(r+ct))/2, r^2 F/c = (G(r-ct) - G(r+ct))/2, scored on r in [5, 20]
(domain of dependence), L1 relative to the pulse part.  clean_pos_limit = false (D = 1 is the realizability edge;
with the limiter or a 1e-3 background the stages fail on 1 negative cell and the run falls back to BE, L1 0.34-0.56).
| N | L1(E) | L1(F) | order |
|---|---|---|---|
| 512 | 1.970e-1 | 5.747e-2 | |
| 1024 | 5.726e-2 | 1.670e-2 | 1.78 |
| 2048 | 1.536e-2 | 4.479e-3 | 1.90 |
| 4096 | 4.413e-3 | 1.287e-3 | 1.80 |
2-D wedge (512 x 64, theta reflect): L1 1.976e-1 / 5.763e-2 (= 1-D), F_theta = 0 exactly, theta spread of E 3.7e-15.

(b) T4 spherical diffusion, rho kappa_s = 4e4, c = 10, D = 1/3 (R = 100), E0 = e^-10 + exp(-40 (r-3)^2) on [1, 5],
T = 580.8, dt = T/512 at N = 64 scaled with dr; exact u = rE heat-kernel solution.
| N | L1(E) | L1(F) | orders |
|---|---|---|---|
| 64 | 1.485e-2 | 4.346e-2 | |
| 128 | 3.995e-3 | 1.336e-2 | 1.89, 1.70 |
| 256 | 1.017e-3 | 3.706e-3 | 1.97, 1.85 |
| 512 | 2.549e-4 | 1.000e-3 | 2.00, 1.89 |
2-D lateral mode (x (1 + 0.5 cos(pi (theta-th0)/dth)), N x N/4), self-convergence of E against 256 x 64:
32: 4.097e-2, 64: 1.264e-2 (1.70), 128: 2.740e-3 (2.21) [exact order 2 gives 2.07, 2.32 against a 256 reference].

(c) T6a black-body sphere R = 1 (bath 4 at r = 1, flux 1), thin, steady (400 steps of 4 dr/c), sphere closure R = 1:
| N | L1(E) | L1(F) | orders |
|---|---|---|---|
| 64 | 3.008e-3 | 6.686e-3 | |
| 128 | 1.772e-3 | 2.548e-3 | 0.76, 1.39 |
| 256 | 1.110e-3 | 1.063e-3 | 0.67, 1.26 |
Boundary limited: N128 first cell E/E_ex 1.090, F r^2 0.770, second 1.002 / 0.975, third 0.969 / 1.020, interior
0.9997-0.9999.  Cause: the upwind stress against the isotropic bath ghost (D_g = 1/3, y_g = 0) at x1min is not the
half-range stress of the emitting surface (same issue as the open-face stress on Cartesian; boundary code not touched).

Energy (closed: reflect both ends, sum E V over the run): difc 1-D N128, 2047 steps: 8.4e-15; latc 2-D 32x8: 8.9e-16,
64x16: 2.2e-16.
Cartesian regression (rw P1/3 tau10 N256, 60 cycles, cc13x 2a68179c vs geo2): max |diff| 7.5e-15 of 1.01 (the RC
coefficients and the tau distance are now computed from the centre distances: round-off, not bitwise).

## 4. Cubed sphere plan (no code)
1. Multi-block first: production sp/cs meshes have many MeshBlocks; the clean solve is single block.  Every Apply,
   slope-weight and closure step needs a ghost exchange of the cell vectors (cln_x slots K/X/P/S..., slope weights,
   D tensor, face low-order flags) through MeshBoundaryValuesCC, global dot products, block-local line preconditioner.
2. Basis: cell vectors in an orthonormal local frame (r_hat, e_xi, r_hat x e_xi) per cell; D from vet_gd's cs kernel
   (VetGdCsMoments, about r_hat) rotated into it (CsTanFrameInit gives the frame).
3. Momentum rows in Cartesian-flux form: per face, the stress P n A converted to Cartesian components and summed
   (exact linear-momentum conservation, no Christoffel terms, well-balanced for a uniform stress by construction since
   sum_f n_f A_f = 0); project the Cartesian divergence back onto the local frame.  E row: plain sum A_f fE / V.
4. Seams: one face flux per seam face, built from the two cells' states in a common frame (the neighbour's y and D
   rotated by the panel-to-panel transform that the cs vector halo already applies), canonical seam areas
   (cs_area_eff) and centre distances (cs_dxf_eff) for RC and tau; the PLM slopes across seams from the rotated ghosts.
   Low-order marking must be mirrored on both sides of a seam face.
5. Tests: radial pulse and T4 on the full cs shell (all 6 panels; y_tangential stays round-off at seams and corners),
   sph_shell lat_l mode with lat_axis = z / edge / vertex (self-convergence), closed-domain energy, then T6a once the
   boundary stress is fixed.
