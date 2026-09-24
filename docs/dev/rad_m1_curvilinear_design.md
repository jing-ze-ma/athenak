# Implicit M1 on the curvilinear grids (spherical polar and cubed sphere): design

Status (2026-09-24): S0 and S1 DONE (sp wedge, Eddington), S2 in progress (branch m1-sp2).
POLES DEFERRED by the user: sp work is on a wedge (theta away from 0 and pi) until further
notice; everything that needs the pole (sect. 1.5 pole stiffness, 2.2, the polar-ring
preconditioner of S4) waits. See sect. 8 for what was built and learned. File:line
references in sects. 0-7 are to be02c647 (the design date) and have drifted since.

Scope: `<rad_m1>/transport = implicit` with `closure = eddington | m1 | vet_sc`,
`implicit_enthalpy = plm`, `implicit_vimp`, `time_scheme = hesdirk2`, BiCGStab with the
radial line preconditioner, `implicit_halo_mpi` / `implicit_halo_direct`, on

* **sp**: spherical polar, `x1 = r`, `x2 = theta`, `x3 = phi` (`Mesh::use_spherical_polar`,
  polar boundary `use_polar_boundary`), orthogonal, radially stretched in production;
* **cs**: gnomonic equiangular cubed sphere, `x1 = r`, `x2 = xi`, `x3 = eta`, six panels
  (`Mesh::use_cubed_sphere`), angular metric `g = [[1, c],[c, 1]]`, `c = cos alpha =
  Coordinates::cos_cell`, `|c| <= 1/2`.

Both grids have `x1` radial and orthogonal to the two angular directions. That single fact
is what makes the design cheap: the exact radial line solve, the radial face equation and
the radial halo are the SAME on sp and cs. Only the transverse faces differ.

---

## 0. Two things that are wrong TODAY (fix first, whatever else is done)

1. **No refusal.** `rad_m1` never checks `use_spherical_polar` / `use_cubed_sphere`
   (grep: the only geometry guards are `multilevel`, plus the halo-direct and Krylov
   `halo_mpi` guards at `rad_m1_implicit.cpp:504-506, 2665`, `rad_m1_krylov.cpp:198`).
   Explicit and implicit M1 on sp/cs therefore run silently with Cartesian arithmetic:
   divergence `/dx` from `mb_size` (`rad_m1_update.cpp:44-49`, `rad_m1_fluxes.cpp:240,295,
   351`), 77 uses of the uniform `mbsize.dx1..3` in `rad_m1_implicit.cpp`. On sp the
   `mb_size.dx2` is `Delta theta` in radians, so the transverse transport is wrong by
   `r` (and `r sin theta`) outright.
2. **The halo slot trap.** `MeshBoundaryValuesCC` treats components `v == IVY (2)` and
   `v == IVZ (3)` of ANY array as the tangential vector pair: seam `TransformMomentum`
   (`bvals_cc.cpp:170-200, 405-435`) and the polar sign flip (`bvals_cc.cpp:295-299`).
   For `pbval_u` that is right by luck (`M1_F2 = 2`, `M1_F3 = 3`, `rad_m1.hpp:42-46`),
   provided F is stored COVARIANT on the cs (as the hydro momentum is). For the implicit
   scratch halos it is wrong: `M1HaloCompT` (`rad_m1_implicit.hpp:440-457`) puts
   `M1_IW_F1` (radial, a scalar under both maps) and `M1_IW_KT` (a scalar) in slots 2, 3,
   so they would be sign-flipped at the pole and sheared at a seam, while the true
   tangential pairs `N2,N3` (slots 6,7), `A2,A3` (9,10), `V2,V3` (12,13) would pass
   untransformed. Same for `thq` (`M1_NHALO_Q = 4`, slots 2,3 = F1, KT).

**Stage 0** (small, 0.5 agent-day): a startup fatal for `rad_m1` on sp/cs until a stage
lifts it; and a per-exchange **component-role table** in `MeshBoundaryValuesCC`
(`role[v] in {scalar, tan_a, tan_b}` with the partner slot), defaulting to the current
`v == IVY / IVZ` rule so every existing exchange is bitwise unchanged; the rad_m1 scratch
exchanges pass their own table.

---

## 1. The equations in the code's form

### 1.1 What the code solves (Cartesian, `rad_m1_implicit_design.md` sect. 1, 11-13)

Backward Euler (or the hesdirk2 stages) for `E` in cells and the face-normal reduced flux
`F0_f` on faces, `F = F0 + A`, `A = v E + v.P` (enthalpy flux, upwinded, `plm` by
deferred correction), closure `P = D E` lagged per Picard pass,
`D = (1-chi)/2 I + (3chi-1)/2 n n`:

```
F0_f' = th_f [ F0_f^n - c^2 dt ( D_dd,R E'_R - D_dd,L E'_L )/dx_d
               - c dt v_f g0_f - c^2 dt (sum_{e!=d} d_e P_de)_f ]      th_f = 1/(1+c dt k_f)
E'_i + (dt/dx) sum_d [ (A'+F0')_{+} - (A'+F0')_{-} ] = E^n_i + dt c G_E,i'
```

The x1 couplings form the tridiagonal `(TA,TB,TC)` solved exactly per column (the
preconditioner, `ImplicitTridiagSolve` `:2149`); the transverse face couplings are the
`CJM..CKP` entries (7-point M-matrix); `implicit_offdiag = operator` adds the
off-diagonal Eddington terms as the 19-point stencil (`M1StIdx` `:2948`,
`ImplicitStencilBuild` `:2988`, slots 15..18 are the (j,k) edges).

### 1.2 The coordinate-free split that makes both grids easy

Write the uniaxial tensor as an isotropic part plus a rank-one part,

```
P = p I + q n n,     p = (1-chi) E/2 >= 0,     q = (3chi-1) E/2 >= 0   (chi in [1/3,1])
div P = grad p + div(q n n)
```

* `grad p` has NO curvature term on any grid. Its normal component at a face is a
  two-point difference over the metric distance between the two cell centres, plus (cs
  only) the non-orthogonal cross term. Positive coefficients on both sides: M-matrix.
* `div(q n n)` carries ALL the curvature. Its normal-normal part along a face's own axis
  can be written as a two-point difference with an INTEGRATING FACTOR (below); the rest
  (its tangential-derivative parts and the curvature terms of the other components) is
  lagged into the right-hand side exactly as `sum_{e!=d} d_e P_de` is today.

For `closure = eddington` (`chi = 1/3`, `q = 0`) the whole curvature disappears from the
flux equation: the face flux is `-c^2 dt th_f d_n(E/3)`. That is why stage 1 is cheap.

### 1.3 Energy equation (both grids)

Replace `/dx_d` by areas and volume, as the hydro does (`hydro_update.cpp:61-83`):

```
E'_i + (dt/V_i) sum_f s_f A_f (A' + F0')_f = E^n_i + dt c G_E,i'
```

`V = Coordinates::volume`, `A_f = area.x1f/x2f/x3f` (`coordinates.hpp:69-70`), built by
`CoordSphericalPolar` (`coordinates.cpp:1580`) / `CoordGnomonicEquiangle` (`:476`),
including the radial stretch. Every face contributes `+A_f nu_f/V_i` to one diagonal and
`-A_f nu_f/V_i` to one off-diagonal, so the row structure and the M-matrix property of
sect. 11 of the implicit design are unchanged; the column sums become `V`-weighted
(conservation: `sum_i V_i E_i` changes only by boundary fluxes and sources, telescoping to
round-off because both cells of a face use the same `A_f F_f`).

### 1.4 Radial face (identical on sp and cs)

`r` is orthogonal to both angles, and `tr D = 1`, so with `n_r` the radial component of `n`:

```
(div P)_r = d_r p + (1/r^2) d_r( r^2 q n_r^2 ) - q (1 - n_r^2)/r + [d_t (q n_r n_t) terms]
```

Discretised at face `i+1/2` (`r_f` face radius, `dr_f = x1v(i+1) - x1v(i)`):

```
F0_f' = th_f [ F0_f^n - c^2 dt { (p'_R - p'_L)/dr_f
                                + (r_R^2 Q_R E'_R - r_L^2 Q_L E'_L)/(r_f^2 dr_f) }
               - c dt v_f g0_f - c^2 dt RHS_f ]
Q = (3chi-1)/2 n_r^2  (lagged),   p' = (1-chi)/2 E'  (chi lagged)
RHS_f = lagged: -q(1-n_r^2)/r + the tangential-derivative terms of q n_r n_t
```

Properties:

* **M-matrix at any dr and dt**: both brackets are two-point differences with positive
  weights. (The naive form `d_r(D_rr E) + (3D_rr - 1)E/r` with a face-mean `E` puts
  `-(3D-1) dt/(2r)` on the upwind off-diagonal and loses the M-matrix for `dr > r`, and
  it is only first-order exact for `E ~ r^-2`.)
* **Exact steady free streaming**: `chi = 1, n = r_hat` gives `F0 ~ -d_r(r^2 E)/r^2`, whose
  discrete steady state is `r^2 E = const` to round-off (gate T-S2).
* Diffusion limit: `q -> 0`, the Cartesian 3-point operator in `r` with the `r^2` areas.
* `D_rr` of the face, `th_f` and `dr_f` enter `TA/TB/TC` only; the line solve, the
  gathered-stack path and the PCR path (`:1998-2580`) are unchanged in structure.

### 1.5 sp transverse faces

`theta` face at `j-1/2` (arc length `dl = r (theta_j - theta_{j-1})`, from `dx2`/`dxface`):

```
(div P)_theta = (1/r) d_theta p + [ (1/(r sin)) d_theta( sin(theta) q n_th^2 ) ]
                + lagged: (1/r^3) d_r(r^3 q n_r n_th) + cot(theta)(-q n_ph^2)/r + d_ph terms
```

The pole face (`theta = 0, pi`) has `area.x2f = 0`: it drops out of the divergence and
its `F0` is set to 0 (reflecting in the face-normal sense, exactly as the hydro's pole
face carries no mass flux). `cot theta` is evaluated only at interior faces and cell
centres, never at the axis. The bracket is a two-point difference with the positive weights `sin theta_{j-1}`,
`sin theta_j` (M-matrix, as the `r^2` of sect. 1.4). For Eddington the bracket and the
lagged terms are all zero.

`phi` face: `dl = r sin(theta_j) dphi`, `(div P)_phi = (1/(r sin)) d_phi p + lagged`.
Periodic in phi; no new machinery.

Pole halos: the polar ghost rows (`bvals_cc.cpp:295-299`: `j` mirrored, `phi + pi`,
tangential pair sign-flipped) are the right transform for every lagged tangential pair
once the role table of sect. 0 exists; `E, chi, KT, G0, F1(radial), a_1, v_1` are
scalars across the pole. This is consistent with the hydro, which uses the same ghost.

**Pole stiffness (the main sp risk).** In an optically thick polar row the phi coupling
`nu_phi ~ c dt/(3 k (r sin theta dphi)^2)` exceeds the radial one by
`(dr/(r sin theta dphi))^2` (x 100-1000 in the polar rows at production resolution). The
row stays an M-matrix (positivity holds) but the radial line preconditioner does not see
it, and BiCGStab iteration counts will be set by the polar rows. Mitigations, in order of
cost: (a) nothing, measure (the Krylov operator already has the exact phi terms);
(b) a second preconditioner sweep: cyclic tridiagonal in phi on the `n_pol` polar rows
only (needs the full phi ring in one MeshBlock; sp production blocks usually are full in
phi near the poles, check the layout); (c) the hydro's polar averaging idea applied to
the radiation (phi-average `E` over the first ring after the solve: conservative, kills
the stiff mode, costs accuracy only at `l ~ nphi/2`).

### 1.6 cs transverse faces

On a `xi` face (unit normal `n_xi`, area `area.x2f`, arc length `dl`, `s = sin alpha`,
`c = cos alpha` at the face, `Coordinates::sin_face_xi/cos_face_xi`) the normal
derivative of a scalar is (same as `cs_implicit_transverse.md` sect. 2)

```
d_n p |_xi = ( d_l,xi p - c d_l,eta p ) / s
```

* `(p_R - p_L)/(dl s)` is the implicit two-point part (M-matrix, the `CJM/CJP` slots);
* `-(c/s) d_l,eta p` at the face is the CROSS term: the four (j+-1, k+-1) neighbours, i.e.
  exactly the 19-point stencil edge slots 15..18 that `implicit_offdiag = operator`
  already allocates. Two options: (i) lagged into `TRHS` (M-matrix kept, positivity
  guaranteed, Picard contraction factor `|c| <= 1/2`); (ii) inside the Krylov operator
  (exact linear solve, positivity NOT guaranteed where `c` is large and `E` jumps). Build
  (i) as default, (ii) behind `implicit_offdiag = operator`, same as today's split.
* `div(q n n)` on the tangent plane has Christoffel terms of the gnomonic chart. Do NOT
  hand-derive them: evaluate the lagged part as a face sum of Cartesian vectors, the
  way `SrcTermsCurvilinearWB` (`coordinates.cpp:1350`) and its cache `wb_geom`
  (`BuildWBGeometry` `:1221`, 60 entries per `(m,k,j)`: cell triad, four tangential face
  normals and triads) do for the hydro momentum: build `P = pI + q nn` as a Cartesian
  tensor from the cell triad, form `(1/V_dual) sum_faces (P . N) A` about the face, and
  project on `n_xi`. The same kernel serves the sp polar rows and removes every `cot
  theta` (it is what `sp_cart_all_momentum` does for the hydro).
* Only the implicit two-point part and the lagged part depend on the chart; the radial
  face (sect. 1.4) is unchanged.

**Components and the gas coupling.** Store the cell F covariant on the unit panel basis,
like the hydro momentum (then `pbval_u` is already correct, sect. 0). The face-normal
`F0` on a `xi` face is `F . n_xi = s F^xi` (contravariant), so the cell vector is rebuilt
from the two face means by `F^xi = F0_xi/s`, then lowered with `g`; the same algebra as
`TransformFieldToDstNormals` (`cubed_sphere.hpp:237-316`). The radiation force handed to
the hydro (`rad_m1_coupling.cpp`) is then in the hydro's own basis. On sp F is in
physical orthonormal components, as the hydro momentum.

---

## 2. Seams (cs) and poles (sp) in the implicit operator, the halos and the Krylov vectors

### 2.1 cs seams

What exists: the cell-centred exchange across a seam does the index permutation, the
along-seam RESAMPLE (interpolation, because a ghost centre of one chart is not a cell
centre of the other) and `TransformMomentum` on the tangential pair
(`bvals_cc.cpp:140-200, 405-435`, `cubed_sphere.hpp:158-236`); the cube-vertex corner is
not exchanged but filled by `cs_vertex_fill_cc`. `flux_seam_cc.cpp` makes a seam face
flux single-valued by averaging the two outward fluxes (faces pair one-to-one, no
resample needed; momentum flux via `TransformMomentum`).

Design:

1. **Two-point seam coupling without resample.** The implicit normal flux across a seam
   face needs only the two cells that share the face. Those pair one-to-one (the
   `CS SEAM FLUX MISMATCH` gate in `pgen/cs_test.cpp` shows the faces coincide to
   round-off), and by the reflection symmetry of the cube about the edge plane the two
   half-distances centre-to-face are equal. So the seam face uses the neighbour's OWN
   cell value, not the resampled ghost: `nu_f (E_B - E_A)` with bit-identical `nu_f` on
   both sides, and the seam face is single-valued and conservative by construction, and
   the operator's seam coupling is symmetric. This needs a "no resample" mode per
   exchange object (today only the global diagnostic env `CS_NORESAMP_CC`,
   `bvals_cc.cpp:66`) for the Krylov vector `krw` and the `thq` scalars, and the
   `ImplicitHaloCopy` must write the first ghost layer from it. Verify the half-distance
   symmetry numerically in stage 3 (one gate), do not assume it.
2. **Resampled ghosts for everything wider**: the cross term (sect. 1.6) and the lagged
   `div(q nn)` read diagonal and second-layer ghosts; they come from the ordinary
   resampled exchange of `thw` with the role table (the tangential pairs `N2/N3`,
   `A2/A3`, `V2/V3` transformed as covariant vectors; `n` is a unit vector: store it
   covariant, renormalise after the transform). `rad_tr_halo_faces_only`-style dropping
   of the x2x3 edge buffers must be refused on the cs (the cross term reads them).
3. **Persistent face state.** `F0` on transverse faces is state (restart, hesdirk2 stage
   memory). On a seam face each panel stores its own copy; after the accepted solve run
   the `flux_seam_cc` average (face-normal F0 is a scalar up to the sign `sgn` it already
   handles) so both copies agree bitwise; with (1) they already agree and the average is
   a no-op check.
4. **Krylov vectors.** Only one Krylov vector is exchanged per operator application
   (`ImplicitKrylovHalo` `:2743`, 1 component, scalar): no transform, no resample (mode 1).
   Global reductions are unchanged. BiCGStab does not need symmetry, but mode (1) keeps
   the seam part of the operator symmetric, which helps.
5. **halo_mpi / direct halos** (`ImplicitHaloDirect` `:2707`, init `:2661`, the
   `impl_halo_mpi` default `:498-506`) are excluded on seams because the direct copy
   kernel knows only same-orientation neighbours. Extension: per (MeshBlock, direction)
   store the seam index map `(aj,bj,ak,bk,swap)` from `GetPanelBoundary` and, for the
   tangential slots, the `SeamXform` geometry (only needed in `thw`, not in `krw`/`thq`
   under mode 1). The `krw`/`thq` direct copy across a seam is then a permuted copy,
   which is the case that matters for cost (it runs 2x per Krylov iteration). `halo_shell`
   (`:925-929`) stays off on the cs. Medium effort, pure performance; do it after the
   physics gates, and measure on apudev.
6. **Cube vertices**: the 3-panel corner is read only by the cross term and the lagged
   terms, same exposure as the explicit `radcond2/3` and the cs ADI (accepted there).

### 2.2 sp poles

* Pole face flux `F0 = 0`, area 0 (sect. 1.5). No line crosses the pole (the lines are
  radial), so the singular-line problem of the ADI (`cs_implicit_transverse.md`: sp keeps
  its fatal) does not arise.
* Ghosts across the pole: the existing polar exchange (mirror in `j`, `phi + pi`,
  tangential sign flip) with the role table; the Krylov vector is a scalar (plain mirror).
* `impl_halo_mpi` / direct halo across the pole: add the polar map (mirror + `phi+pi`
  shift: a permuted copy) to the direct kernel, same work as the seam maps.
* Lagged terms in the two polar rows: use the Cartesian face-sum (sect. 1.6 last bullet),
  which has no `cot theta`, mirroring `sp_cart_polar_momentum` for the hydro.
* Radiation force: consistent with the hydro's sp momentum formulation; if the hydro
  runs `sp_cart_all_momentum`, the force must be handed in the same Cartesian-projected
  form in the polar rows (`SrcTermsSphericalPolarCartRows` `coordinates.cpp:1799`).

---

## 3. The VET formal solution on curved grids

Today (`rad_m1_vet.cpp:1-80`): grey SC, rays with FIXED directions in the Cartesian
frame (double-Gauss in `mu_1 = Omega.x1` times uniform azimuths about x1), swept layer by
layer in x1, the upwind intensity bilinearly interpolated on the PLANE of the upwind
layer's cell centres with periodic x2/x3 wrap; exact multi-block/multi-rank version by
band exchange (`VetSweepMB` `:2222`), full tensor option (`VetFullTensor` `:3017`).
LIMITS: periodic x2/x3, uniform mesh. What carries over: the layer-by-layer radial sweep,
the upward/downward pair of launches, the band exchange pattern, the uniaxial projection
(`vet_axis = flux`), the moments. What does not: fixed ray directions (a straight ray
changes its local `mu` from shell to shell) and the plane geometry.

| option | what | cost per VET update | accuracy | parallel ordering |
| --- | --- | --- | --- | --- |
| (A) local-frame SC | rays at fixed LOCAL `(mu, psi)` about `r_hat`; trace back to the upwind SHELL, interpolate laterally on the shell (bilinear, as now), in `mu` (the ray's `mu` at the upwind shell from the conserved impact parameter) and in azimuth after parallel-transporting the local frame (rotation angle from the two triads; across a seam from the panel bases) | ~4x today (16 reads per ray-cell vs 4) | 2nd order in space, 1st-2nd in angle; angular interpolation smooths beams | identical to today: radial layers, band exchange laterally; seams: the band goes through the seam map + azimuth rotation; poles: band through the polar map |
| (B) long characteristics | trace each (cell, direction) ray to the boundary | `N_cells N_ang N_r` ~ 2e11 per update at the He-presn grid (1.2e8 cells, 16 ang, ~100 steps): too slow per step; rays cross many ranks | best (no angular diffusion) | ray passing across ranks; not a sweep |
| (C) curvilinear discrete ordinates | FV `S_N` with angular redistribution terms (`(1-mu^2)/r d_mu` and the frame-rotation terms) | like (A) | positive only with step differencing; ray effects | full 3-D wavefront (KBA); phi and the closed cs panels make the dependency graph CYCLIC -> lagged transverse inflow, i.e. no longer exact |
| (D) column-wise | per radial column, 1-D SPHERICAL SC in `(r, mu)` with the column's own profile (laterally homogeneous), `mu` interpolation at each shell; `D = diag(f_K, (1-f_K)/2, (1-f_K)/2)` about `r_hat`, or about the M1 flux axis | `N_col N_r N_mu`, ~2e9 at the He-presn grid (8-16 `mu`): a few % of a step | exact for spherical symmetry (captures sphericity: `f_K -> 1` at large r, the thing Eddington/M1 get wrong in an extended atmosphere); misses lateral anisotropy (plumes, day-night) | column-local, no lateral communication at all; radial stack as `VetSweepMB` does for x1 splits |

Recommendation: build **(D) first** as `vet_col` (the ray/`mu` bookkeeping of (A) without
the lateral part: (A) = (D) + lateral interpolation + frame transport), validate it
against the 1-D spherical Python reference (Chapman/Hummer-Rybicki sphericity), then **(A)**
only if a lateral-anisotropy case needs it (irradiated dhj day-night, strong plumes). Keep
(B) as an offline Python reference (`tests_ck_sph/real_rays.py`/`rr_lib.py` already trace
real rays through the cs prod3 atmosphere; reuse for VET validation). Do not build (C).

Two-stream precedent: the curved-column r^2 handling already exists for the two-stream
RT (area-weighted intensities `J = A I`, areas relative to the column top; see
`two_stream_rt.hpp:333-357` and the `rt_plane_parallel` branch `:1689-1712`); `vet_col`
should use the same relative-area scaling (the absolute-area pivot failure of the He
mode-3 fix is the lesson).

Bottom boundary for upward rays on a sphere: the diffusion intensity of the bottom cell as
now; core-intersecting rays (`mu` > `mu_core` at the inner shell) end at the inner
boundary, the rest (tangent rays) come from the opposite side of the shell; in (D) they
are mirrored in `mu` at the turning point `r = p` inside the column (standard 1-D
spherical SC).

---

## 4. Lifting LIMIT 4 (radially split columns)

State: `implicit_partition = gather` is parsed (`rad_m1_implicit.cpp:694-705`), with its
checks (`:765-780`), `ImplicitPartitionInit` (`:1069`) and `ImplicitGatherSolve`
(`:1998`) in the tree from 3b phase 1 (implicit design sect. 8). The gather moves
ASSEMBLED rows `(TA,TB,TC,TR)`, so once the curvilinear coefficients are in the rows it is
geometry-agnostic; its checks are about block counts (verify the stretched radial grid
passes `:769`). Options on curved grids:

1. **Block-local preconditioning** (each block solves its own radial segment; the
   stack-neighbour coupling is left to BiCGStab through the ordinary halo). Trivial
   (skip the gather), no new communication. Risk: in optically thick layers the
   cross-block-face error mode is smooth in r and poorly damped; iteration counts grow
   with the number of blocks in a stack and with the radial stiffness. Small effort
   (1 agent-day incl. the measurement).
2. **Gather** (existing): one gather + scatter per preconditioner application (2 per
   BiCGStab iteration) along the stack; latency-bound when the stack spans ranks.
   Small effort to validate on sp/cs.
3. **Distributed partitioned line solve** (SPIKE/Wang, the algorithm of
   `two_stream_column_partition.hpp` lifted from threads to blocks): each block
   eliminates its segment with symbolic end unknowns, exchanges 2 rows per column with the
   stack, every block of the stack solves the `2 nblk` reduced system redundantly,
   back-substitutes. One small message round per application, exact to round-off.
   Medium effort (4-6 agent-days).

Recommendation: LIMIT 4 is NOT on the critical path for the global runs planned so far
(the He presn grid keeps the full radial extent, 192 cells, in one block: 192x80x80,
memory `he4-presn-global-plan`). Measure (1) vs (2) on the spherical diffusion test with
2 and 4 blocks per stack on apudev; build (3) only if (2) costs more than ~15 % of a step
and (1) more than doubles the iteration count. `vet_sc`/`vet_col` with split columns: the
stack is swept in turn (`VetSweepMB` already does this for x1 splits).

---

## 5. What exists vs what is new

| item | exists (file:line) | new | effort |
| --- | --- | --- | --- |
| geometry arrays (V, A, dx, face dl, stretch) | `coordinates.hpp:69-83`; `coordinates.cpp:476, 1580` | read them in the implicit kernels instead of `mbsize.dx*` (77 sites in `rad_m1_implicit.cpp`, the f2face/f3face/tcell kernels `:1692-1797`, `ImplicitTransTheta :1486`, asm `:5898`, aphll `:5641`, src `:5805`, stencil `:2988`, od cache `:2603`, vimp `:4736`), iw slots for per-cell geometry or on-the-fly reads | M (3 d) |
| area/volume divergence | hydro `hydro_update.cpp:61-83` | same form in the E row and the post-solve update; explicit M1 (`rad_m1_update.cpp`, `rad_m1_fluxes.cpp`) keeps the fatal | S |
| radial face with integrating factor | Cartesian x1 face (sect. 11 operator) | `r^2` weights and lagged `-q(1-n_r^2)/r` | S (1 d) |
| lagged curvature terms via Cartesian face sum | `SrcTermsCurvilinearWB :1350`, `wb_geom :1221`, `SrcTermsSphericalPolarCartRows :1799` | an M1 kernel on the same cache for `div(q nn)` projected on face normals | M (3 d) |
| cs non-orthogonal face gradient | `cs_implicit_transverse.md` sect. 2; `conduction_transverse.cpp:155-157, 342-420`; 19-pt stencil `:2948-3130` | cross term in TRHS (default) or in stencil slots 15..18 | M (2 d) |
| seam/pole ghosts | `bvals_cc.cpp:140-200, 295-299, 405-435`; `cubed_sphere.hpp` | component-role table; per-exchange no-resample mode | S-M (1.5 d) |
| seam face single-valuedness | `flux_seam_cc.cpp` | apply to transverse `F0` state after the solve | S |
| F covariant (cs) / physical (sp), force to hydro | hydro conventions, `TransformFieldToDstNormals` | cell F from face normals, coupling basis in `rad_m1_coupling.cpp` | M (2 d) |
| BiCGStab, line preconditioner, gather, PCR | `:1998-2580, 3889-4580` | nothing structural; polar-ring preconditioner optional | 0 / M (3 d) |
| direct / halo_mpi halos on seams, poles | `:2661-2745`, guards `:504-506, 2665`, `rad_m1_krylov.cpp:198` | permuted direct copy per direction | M (3 d) |
| time2 (hesdirk2) | `rad_m1_time2.cpp` | none beyond the operator; face state on seams | S |
| restart of face state | 3b sect. 8/11 | transverse `F0` on curved faces; seam copies | S |
| `vet_col` (D) | `VetSweepAng :1175`, moments `:914-1111`, `VetSweepMB` stack order | 1-D spherical SC per column, `mu` interpolation, radial axis projection | M (4 d) |
| local-frame SC (A) | band exchange `:324-722`, `VetSweepMB` | shell geometry, `mu`/azimuth interpolation, frame transport, seam/pole band maps | L (10-15 d) |
| LIMIT 4 | gather `:1069, 1998`; `two_stream_column_partition.hpp` | measure (1)/(2); SPIKE (3) if needed | S / M (5 d) |
| boundary conditions on a shell | `bvals/physics/rad_m1_bcs.cpp`, `SetImplicitX1BC :4580` | inner: imposed luminosity `L/(4 pi r_in^2)`; outer: vacuum/Marshak with the `r^2` face | S |

S = under 1 agent-day, M = 2-5, L = over 5.

---

## 6. Staged plan and tests

Every stage: Cartesian gates bitwise (all curvilinear branches appended as overwrites
behind the geometry flag; the cs ADI lesson: dividing by an exact 1.0 or hoisting into a
branch changes FMA contraction under hipcc); CPU for correctness, cost only on apudev with
the validated env (`HSA_XNACK=1`, `HSA_NO_SCRATCH_RECLAIM=1`), same binary, interleaved.

| stage | content | effort | gates |
| --- | --- | --- | --- |
| **S0** DONE | fatal on sp/cs for rad_m1; component-role table in `MeshBoundaryValuesCC` | 0.5 d | every existing exchange bitwise (hydro, MHD, two-stream, M1 Cartesian) |
| **S1** DONE (wedge) | sp, `closure = eddington`, one radial block per column, no SMR, `halo_mpi` off: geometry in rows, area/volume, ~~pole face, polar ghosts via roles~~ (deferred), shell BCs, physical-component F and force, stretched r | 4 d | T-S1 spherical diffusion; T-S3 Marshak in a thin shell; T-sym; T-bit |
| **S2** IN PROGRESS | sp, M1/Kershaw closures: radial integrating factor, lagged curvature by the Cartesian face sum, ~~polar rows~~ (deferred) | 3 d | T-S2 free streaming 1/r^2; T-S4 star-like atmosphere (M1 closure); T-S5 radwave in a shell |
| **S3** | cs: radial reused; transverse face gradient with cross term; no-resample seam coupling; seam `F0` state; covariant F; vertex corners | 6 d | all S1/S2 tests on the cs, plus seam conservation and half-distance symmetry, l = 2 mode across seams |
| **S4** | direct / `halo_mpi` halos across seams and poles; polar-ring preconditioner if S1/S2 iteration counts demand it | 3-6 d | bitwise equal to the ordinary exchange; apudev timing |
| **S5** | `vet_col` (option D) on sp and cs | 4 d | T-S6 VET sphericity vs 1-D reference; T-S4 with vet_col |
| **S6** | LIMIT 4 on curved grids: measure block-local vs gather; SPIKE if needed | 1-6 d | T-bit across stack splits (gather: bitwise vs one block) |
| **S7** (optional) | local-frame SC (option A) | 10-15 d | T-S6 plus a lateral test (beam/shadow on a shell), vs the offline long-characteristics rays |

Total to a production-usable cs + vet_col: ~20 agent-days (S0-S5); +S6/S7 as needed.

### Tests

* **T-S1 spherical diffusion (Eddington, constant `rho kappa`, thick)**: (a) steady
  shell between `r_in`, `r_out` with imposed `L`: `L(r) = 4 pi r^2 F` constant to 1e-10
  (telescoping), `E = a + b/r` to `O(dr^2)`; (b) time-dependent point-source Green's
  function `E ~ (4 pi D t)^(-3/2) exp(-r^2/4Dt)` on a shell with Dirichlet from the
  analytic solution; L1 convergence at 2nd order in `r`, time 1st (BE) / 2nd (hesdirk2).
* **T-S2 free-streaming point source** (`kappa = 0`, M1 closure, inner `F = cE`): steady
  `r^2 E` constant to round-off with the integrating-factor face; transient front at `c`.
* **T-S3 Marshak wave in a shell**: `dr_shell << r` reproduces the Cartesian Marshak gate
  (existing tests_m1 reference) to the Cartesian error; thick-shell variant vs a 1-D
  spherical run (sp with `nx2 = nx3 = 1` if the mesh allows, else the S1 code on a narrow
  wedge).
* **T-S4 radiative-equilibrium star-like atmosphere** (grey, gas coupled, extended):
  `L_out/L_in = 1` to 1e-4 at equilibrium, `T(tau)` vs a 1-D spherical Python reference
  with the exact sphericity factors (Chapman / Hummer-Rybicki); Eddington vs M1 vs
  `vet_col` distance to the reference measured, not asserted.
* **T-S5 radwave (T10) in a thin spherical shell**: tangential propagation on sp (away from
  and across the pole) and on the cs (across a seam and a cube vertex); dispersion
  relation vs the Cartesian T10 at small curvature; decay rate to the Cartesian error.
* **T-S6 VET sphericity**: `f_K(r)` of `vet_col` vs the 1-D reference; `vet_col` on a
  spherically symmetric state gives the same `D` on every column to round-off.
* **T-sym**: a spherically symmetric state stays symmetric: spread of `E` over each shell
  < 1e-12 relative on sp (pole rows included) and on the cs (all panels, seams, vertices).
* **Seam conservation (cs)**: global `sum V E` + boundary fluxes + sources, drift per step
  at round-off with the no-resample seam coupling.
* **T-bit restart/rank**: restart mid-run bitwise; 1 vs N ranks bitwise (serial vs MPI);
  `halo_mpi`/direct halo on vs off bitwise; gather vs single block bitwise; GPU vs GPU
  rerun bitwise.
* **Positivity**: `E > 0` at `c dt/dx ~ 1e6` in a thick-thin shell with the M-matrix
  defaults (cross term lagged).

---

## 7. Main risks

1. **Silent wrong answers today** (sect. 0): no refusal, and the scratch halos would
   transform scalar slots 2,3. S0 closes both.
2. **Polar stiffness** (sect. 1.5): radial-only preconditioning vs the phi coupling in
   thick polar rows; Krylov counts may be dominated by `n_pol` rows. Measure in S1.
3. **Seam half-distance symmetry** (sect. 2.1): the no-resample seam coupling is only
   single-valued if the two half-distances and areas agree bit for bit; if they agree
   only to round-off, fall back to averaging the two coefficients (still conservative).
4. **cs cross term and positivity**: exact in the operator loses the M-matrix; lagged
   keeps it at a Picard-count cost bounded by `|cos alpha| <= 1/2`.
5. **Stretched radial grid**: the implicit code assumes uniform `mb_size.dx` in 77 places;
   every one must move to per-cell spacing, and the Cartesian path must stay bitwise.
6. **Force balance with the hydro's well-balanced source**: the radiation force in a
   hydrostatic star must be discretised consistently with the hydro's WB/Cartesian
   face-sum momentum source (`force_reference wb_arad`, `sp_wellbalanced_src`,
   `cs_wellbalanced_src`), otherwise a static radiative envelope drives spurious flows.
   Gate: T-S4 with gas momentum on, rms v at round-off level over 100 steps.
7. **VET on curved grids**: `vet_col` ignores lateral anisotropy (dhj day-night); the full
   (A) is large and adds angular diffusion.

---

## 8. Progress log and lessons (2026-09-24)

**S0** (on rt-integration before S1): rad_m1 is fatal on sp/cs unless a stage lifts it; the
per-exchange component-role table in `MeshBoundaryValuesCC` exists.

**S1** (merged 4e428a30; tests_m1/runs_5a_sp_s1/README.md):
* Target grid: the He presupernova model `inputs/hydro/he4_presn_sp.athinput` (bench/wt_he4) is
  ALREADY a wedge: theta = pi/2 +- 45 deg, phi 90 deg, both periodic; r stretched with
  `use_grid_stretch_r_poly`, nx1 = 96, one block along r. So per-cell spacing (risk 5) was part
  of S1: on sp the rows use `Coordinates::area`, `volume` and `dxface`; the uniform `mbsize.dx`
  paths are untouched.
* Files: `rad_m1_implicit.cpp` (sp overwrite blocks), `rad_m1.cpp` (guard), new
  `rad_m1_sph.cpp` (`SphericalS1Check`), `rad_m1_newdt.cpp` (smallest physical width),
  `pgen/tests/rad_m1_tests2.cpp` (`m1_test = sph_shell`).
* Allowed on sp after S1: transport implicit, closure eddington, time_scheme be, one block
  along x1, no SMR, halo_mpi off (default off on sp), implicit_flux central, recon dc,
  trans_limit none, no vimp, offdiag forced none. Fatal: poles (polar boundary or theta
  reaching 0/pi), cs, theta stretch, radial stretch on a Cartesian mesh, and every other
  option. halo_mpi / hesdirk2 / aphll / plm_dc / trans_limit / vimp / tau closure have no
  geometric blocker but divide by the uniform dx or are untested: each needs its own gate.
* Results: T-S1 order 2.00 uniform and 1.99 on the stretched He grid, r^2 F const to 8e-11;
  T-S3 shell - Cartesian = 1/R exactly; T-sym 2e-14 at tight tolerance; restart bitwise;
  BiCGStab 1.8-2.4 inner iterations per solve at tau 2000.
* Open: the plm enthalpy correction is not exact on a stretched grid; risk 6 (radiation force
  vs the hydro's well-balanced source) is tested in S2 (T-S4 with gas momentum).

**Lesson (hipcc FMA, stronger than the sect. 6 note).** The first S1 version wrote the
Cartesian and sp coefficients through per-face select variables (`nu2p = sph ? A/V : nu2`,
Cartesian value identical). CPU was bitwise, but on the GPU the 3-D Cartesian box differed
at round-off from the first step: hipcc contracts FMAs differently once the operands come
through a select. The fix (ced90b9a): every Cartesian expression stays TEXTUALLY as before,
and the sp form is a separate `if (sph)` block that overwrites the result. A ref-vs-ref GPU
rerun confirmed the difference was the code, not noise. Rule for S2-S7: never rewrite a
Cartesian expression, not even into an equivalent select; add overwrite blocks, and gate
the Cartesian path bitwise on the GPU as well as on the CPU.

**Related tools now on rt-integration.** `<rad_m1>/implicit_op_check = K` (72dbb7fd)
compares every operator variant on one vector per row at 1e-12 and checks the halo ghosts;
use it on sp/cs operators as they are added (every variant must agree with the plain
`ImplicitStencilOp`). tst/test_suite/rad_m1 has CPU slab, operator-check and restart tests.
