# Clean radiation scheme (`rad_m1/scheme = clean`): design note

Status: DESIGN, for user review before any code (2026-11-14). No source changed.
User decision (11-14 ~01:45 PST): write a new clean scheme next to prod-1118 instead of adding more fixes.
**Gate 1: the full IDORT suite at 2nd order in space and time.** Beamed tests may be 1st order in space.

---

## 0. Why a new path, and what "clean" means

prod-1118's `rad_m1` is 78k lines with about 300 `rad_m1/*` keys. Most keys added since October do not model physics. Each one
covers a symptom of one of four defects:

| Defect | Symptom keys it produced |
|---|---|
| D1. The implicit transport operator is not positivity-preserving: central face flux in anisotropic regions, plus the full Eddington-tensor (VET) off-diagonal terms | `face_fct`, `face_mmat`, `zalesak_cap_rounds`, `k1_redist`, `time2_posfix`, `vet_col_lat_posfix`, `implicit_pos_*`, E floors, `adm_poslim`, `dtr_poslim` |
| D2. Closure data (D, face weight) are not smooth in time: the weight switches, and rotation jumps | `blend_dtc`, `x2_fields 15`, `xfade`, half-range face memory (`cmem`) |
| D3. Lateral data from `vet_gd` is wrong at beam edges | x1-only face masks (`implicit_thin_js_gd`, `_hcface_gd`, `_bq_gd`, `implicit_hr_face_gd`), `vet_col` handover + LAT slots |
| D4. The local gas–radiation exchange solve is not robust (silent root-find failure, T⁴ tangent leading to E < 0, floor flips) | `opac_lag_passes`, `predictor_tcap`, `fb_local*`, `stage_local_be`, `tsolve_robust` |

"Clean" means: each part below is derived from the equations; each mechanism has exactly one implementation; and the scheme
reads no compensation key. The new path does not read any of the keys in the table. If the clean scheme cannot pass a test
without one of them, that is a finding to report, not a key to add.

## 1. Equations (unchanged)

Grey moments: E in cells and F on faces (staggered, as now), gas energy e, O(v/c) lab-frame coupling as in the current code, and
the VET closure P = D E with D from a formal solution. The exchange is κ_P(T) and κ_E(T). EOS and opacity tables are as now.

## 2. The four parts

### 2.1 Time: H-ESDIRK2 (keep) + ONE extrapolation rule for all closure data
- Integrator: the existing H-ESDIRK2 (docs/dev/rad_m1_time2_design.md). It is 2nd order uniformly in τ in the model: 1.95–2.01
  on 12 radwave cases, L-stable, stiffly accurate. It is not the source of the order losses.
- Every quantity derived from the formal solution (D, the half-range factors h±, the face weight w of 2.2) is ONE field set
  𝒞. 𝒞 is evaluated at the step start from the formal solution and extrapolated linearly to the stage time:
  𝒞* = 𝒞ⁿ + (Δt/Δt_prev)(𝒞ⁿ − 𝒞ⁿ⁻¹). The model shows this is what gives 2nd order. Today it applies to 15 fields
  (`x2_fields`) and not to the weight (hence `blend_dtc`). In the clean path it applies to all of 𝒞, with no options.
- 𝒞ⁿ and 𝒞ⁿ⁻¹ must come from the same direction set. At a rotation step the formal solution is evaluated on both sets
  (this is what `xfade` does, as one extra sweep per rotation). Optimisation for later: run the second sweep only there.
- Rotation interval: at least 1 domain light-crossing time (measured order threshold 0.6–1.2 w/c on the stiff pulse). Every
  production already meets this.

### 2.2 Space: one face flux with one ray-based weight
- Face flux F_f solves the staggered moment equation with a convex face blend:
  F_f = (1 − w) F_f^central + w F_f^hr.
  - F^central: the 2nd-order staggered flux from ∇·(D E), the present diffusion/central face.
  - F^hr: the half-range upwind flux c (h⁺ E_L − h⁻ E_R), with h± from the rays. It is 1st order and positive.
- **The weight comes from angular information only** (user rule; no cell τ): w = W(a), where a is the ray anisotropy at the face.
  - The candidate is a = C1 = ⟨|I − S|⟩/J, the existing `fse` measure.
  - W is a fixed smooth ramp: 0 for near-isotropic fields (diffusion, linear waves, the rw grid), giving 2nd order;
    1 for beamed fields, giving 1st order and positivity, which the user allows for beams.
  - W is part of 𝒞, so it is extrapolated in time (2.1). It is never switched.
- Off-diagonal D terms: positivity fails here (D1) and at cs vertices. Clean treatment, to be decided by the 1-D/2-D
  prototype (3.1):
  - (a) Split D = D_diag + D_off. D_diag goes into the implicit operator, which is then an M-matrix after the F
    elimination. D_off goes in as a flux limited against the diagonal part (two-point flux plus a limited cross term).
  - (b) Lag D_off by one Picard pass. This is not positive in general, so it is the fallback only.
- **One positivity mechanism, as a safety net only.** With w → 1 in beams, the operator is positive wherever the central flux is
  not used. If a smooth region still produces E < 0 (it should not; the tests count events), apply ONE Zalesak FCT pass on the
  central part of the flux, relative to the half-range solution. Its firing count is printed every cycle and must be 0 on the
  smooth IDORT tests.

### 2.3 Exchange: exact local backward-Euler per stage
- Per cell, with neighbours frozen inside a Picard pass: eliminate E analytically,
  E′ = (Ẽ + c γΔt κ_P a T′⁴)/(1 + c γΔt κ_P).
  E′ > 0 by construction for any T′ > 0, so no floor is needed.
- The gas energy equation becomes one monotone scalar equation g(T′) = 0. Solve it by Newton with a bracket (Illinois in
  log T), with κ evaluated at the trial T. Failure is a hard error with a counter, never "converged".
- The `local` variant of m1replay (fork/replay-1114) is this solver. It is tested offline on the hard-cell corpus (RSG cell
  (43,31,351), R136 silent failures, Milne cold start) before it goes into the code.
- Energy is booked from the assembled row (e + E conserved to round-off per cell); this is checked in the 0-D tests.

### 2.4 Closure: one SC family, one direction family
- Short characteristics along a fixed direction set: `vet_gd` on sp/cs, and the same builder (vet_sc with
  `vet_dirset = gd`) on Cartesian. Geodesic level 2/3, as the user prefers; no HEALPix.
- One face-mask rule on all faces. The x1-only masks (D3) are not carried over. If the lateral data is wrong at beam edges,
  that error belongs to the SC interpolation and is fixed there.
- Ray effects in thin regions (spherical pulse: lateral spread grows with resolution, order 1.1) are an angular
  discretisation error. Rotation (2.1) is the cure. It is not counted against the spatial order of the moment solver.

## 3. Build order (each step has a pass/fail test before the next starts)

1. **1-D prototype (Python, CPU, hours):** E/F staggered + H-ESDIRK2 + 2.2 blend + 2.3 exchange, closure from a 1-D
   two-direction-set SC. Tests:
   - rw linear wave at 6 (P, τ) corners: space+time order ≥ 1.9, and time-only order ≥ 1.9;
   - Marshak wave and Su–Olson (positivity, front position);
   - a 1-D beam (positivity, 1st order allowed).
   Decides W's ramp and the D_off option (2.2).
2. **0-D exchange in C++ (shared header, the replay harness):**
   - corpus pass rate 100 %;
   - e + E to round-off;
   - Milne cold start monotone.
3. **AthenaK Cartesian path `rad_m1/scheme = clean`:**
   - new files `src/rad_m1/clean/*.{hpp,cpp}`, listed in src/CMakeLists.txt;
   - reuses the SC formal-solution kernels, the Krylov solver + preconditioner, and unify step 4's x1 row builder;
   - target < 5k new lines;
   - reads no compensation key (fatal if one is named together with `scheme = clean`).
4. **IDORT gate 1 (Cartesian), section 4.**
5. sp/cs ladders with production-like decks, then 200-cycle windows on real RSG/HG/R136 states vs prod-1115/prod-1118, then
   the corner gate with radiation, then the accretor and MHD.

## 4. IDORT gate 1: what "pass" means

The reference for each metric is the PROD19/FIN3 numbers in /viper/ptmp2/jinma/blend_1022/PHASE2_RESULTS.md and
idort_gate_1111/rwcmp.md (ref columns).

| Test | Metric | Pass |
|---|---|---|
| rw linear-wave grid, 18 (P_rad/P_gas, τ) cases, N 16..256 with Δt ∝ Δx | order p(ρ), p(E) on the last two pairs; L1 at 256 | p ≥ 1.8 in every case; L1(256) ≤ 1.2 × ref |
| rw time-only ladder (N fixed, Δt/2ᵏ), 6 corner cases | p(E) | ≥ 1.8 |
| homogeneous sphere core_u | L1(r > 1.1), E/J − 1 per shell | ≤ ref |
| non-LTE scattering nltek1–8 | L1 vs ε | ≤ ref |
| shadows shd, shdg, hcyl | E_L1, umbra levels, hcyl E/F L1 | ≤ 1.2 × ref; umbra no worse |
| crossing beams cb, cbg, cbh, ws66 | peak, FWHM, lum_top, E_L1; positivity | beams may be 1st order in space; no E < 0, no NC, no FATAL; cbh must RUN (it is FATAL on PROD) |
| every test | NON-CONVERGED, FCT firings, BE redos, dt retries | 0 NC; 0 FCT firings on rw/core/nltek |

## 5. Risks and open choices (for the user)

1. **W ramp.** Where isotropic turns into beamed, w is between 0 and 1, and order there is somewhere between 1 and 2. The
   rw grid at τ 0.01–0.1 is the test. If it shows order loss, the choice is the ramp shape, never τ.
2. **D_off positivity on cs vertices:** option (a) needs a 2-D prototype step.
3. **Cost:** extrapolating all of 𝒞 costs memory (2 copies); the second sweep at rotations costs about 18 % on R136
   (measured for xfade). Both are acceptable for gate 1. Speed comes after correctness.
4. **Production RSG/HG/R136 keep their current binaries** until step 5 passes.

---

## Appendix A. CleanSolve discretisation (Cartesian, one stage)

One stage solves U = Ũ + h R_imp(U) with h = γΔt (or Δt for BE). Ũ is the stage's old vector from the time2 driver.
Unknowns: E_c (cells), F_f (faces, one per direction), gas T_c. The velocity v is frozen at the stage start inside the solve.
The O(v/c) source terms are the code's mixed-frame terms, linear in (E, F) for fixed v.

**Face equation** (x-face f = i+½; y and z alike):

  F_f (1 + h c ρκ_F,f) = F̃_f − h c² [ ((D_xx E)_{i+1} − (D_xx E)_i)/Δx + ∂_y(D_xy E)|_f + ∂_z(D_xz E)|_f ] + v-terms

- ρκ_F,f is the face average of the neighbours.
- The cross derivatives use the 4-cell average at the face.

The flux that enters the E equation is the blend

  𝔽_f = (1 − w_f) F_f + w_f c (h⁺_i E_i − h⁻_{i+1} E_{i+1}).

- w_f = W(max(C1_i, C1_{i+1})), extrapolated with 𝒞.
- W is a C¹ smoothstep from 0 at C1 = c_0 to 1 at C1 = c_1.
- c_0 and c_1 are fixed constants, set once from the IDORT rw/beam tests and never tuned per run.

**E equation:**

  E_i (1 + h c ρκ_P,i) + h Σ_d (𝔽_{d,+} − 𝔽_{d,−})/Δ_d = Ẽ_i + h c ρκ_P,i a T_i⁴ + v-terms

Substituting F_f from the face equation gives a linear system in E alone.
- The diagonal-D part (diffusion, weight 1 − w) and the half-range part (upwind, weight w) both give non-positive
  off-diagonals and a positive diagonal, so they form an M-matrix.
- Only the (1 − w)-weighted cross terms D_xy, D_xz can break this. In beams w → 1 and in near-isotropic fields D_xy → 0,
  which is why W is a function of anisotropy.

**Exchange, nonlinear iteration (k):**
1. The local exact T: for the current E_k, solve ρ c_v (T − T̃) = −h c ρκ_P(T) (a T⁴ − E_k) + work terms for T > 0.
   Its right-hand side is monotone decreasing in T, so there is exactly one root. Use Newton plus a bracket in log T.
   This T_k is consistent with E_k.
2. Linearise a T⁴ ≈ a T_k⁴ + 4 a T_k³ (T − T_k) around that consistent T_k, eliminate T per cell, and solve the E system
   (Krylov) for the step δ.
3. Update E_{k+1} = E_k + λ_i δ_i, with λ_i = min(1, 0.9 E_k,i/|δ_i|) where δ_i < −0.9 E_k,i, so that E stays positive.
   - λ only shapes the iteration path. The fixed point is the exact stage solution, which is unique and positive.
4. Converged when max |δ|/E < tol and the gas residual with κ(T) < tol.
   - A cell whose T root find fails is a hard failure. It is counted and never passed as converged.
5. F_f from the face equation using E_{k+1}. Gas energy and momentum get minus the radiation source integrated over the
   stage, so total energy and momentum are conserved to round-off.

**Boundaries:** the boundary face flux is half-range: 𝔽_b = c h⁺ E_in − F_inc.
- F_inc is the incoming half-range flux of the formal solution's boundary data: zero for vacuum, a beam inflow, or an
  isotropic wall.
- Reflecting walls: 𝔽_b = 0.

---

## Appendix B. As built (11-14 overnight), and what changed from Appendix A, with the evidence

Branch `clean-1114` (local worktree /viper/ptmp2/jinma/clean_1114/wt), file `src/rad_m1/rad_m1_clean.cpp`.
Progress and numbers: /viper/ptmp2/jinma/clean_1114/PROGRESS.md.

1. **Lab-frame face flux.** The face flux is the lab-frame F with the source −cχ(F − A), A = vE + v·P, so
   F = θ(F̃ − h c²∇·(DE)) + (1−θ)A.
   - A first version evolved the comoving F0 with a separate enthalpy flux. That solves the "vetqsc" model: at τ = 1000 it was
     0.7 % too fast and damped 23 % too little.
   - The rw reference and the exact transfer solution are the lab form (vetqs_d).
2. **Implicit gas velocity.** The velocity in the drag and in A is implicit: a fixed point on dm inside the Newton loop.
   - With v lagged at the stage start, every radiation-dominated rw case (P_rad/P_gas ≥ 1/3) was first order.
3. **Closure from the state.** The formal solution uses the state's own source. `vet_source = relaxed_esc` and
   `vet_scatter_j = relaxed` (in the production keys) depend on dt and shift the closure in time by O(dt).
   - Measured: time self-convergence 1st order with them, 2.00 without, equal to production (P_rad/P_gas 1/3, τ 10).
   - Frozen-closure control: 2.00.
4. **Ray–moment consistency** (high-order/low-order) in the half-range flux: + q w c (H_f − upwind of J).
   - H_f is a limited upwind (MUSCL, minmod) reconstruction of the rays' half-range fluxes (`clean_hface = muscl`), clamped to
     the donor bound.
   - q = (min(E,J)/max(E,J))² is the quasi-steadiness gate.
   - Without the gate, beams stay sharp but every stage goes negative.
5. **Cells the rays leave empty** (J < 10⁻³ E) take the closure of their ray-described neighbours (8 fill passes), isotropic
   beyond.
   - A closure from the cell's own flux sustained a spurious sideways stream along the injection row (cbg).
   - Isotropic alone random-walked the beam radiation (cb_r4).
6. **Linear-solve accuracy.** Values within 10 × tol × max|B| of zero are zero to solver accuracy (near-vacuum cells next to
   beams). Real negatives are counted, and a stage with negatives falls back to backward Euler.

Gate-1 status (CPU, single block):
- Pass:
  - rw grid: 13/18 cases at 2nd order with errors 0.9–1.04 × PROD19. Three marginal fails: thin τ = 0.1 E-order 1.3–1.4 (reference
    1.55); P 10/3 τ 0.1 at 1.16–1.40 × reference. Two τ = 0.01 cases still running.
  - nltek1–8.
  - shd (0.79 ×), shdg (0.87 ×).
  - cbg (1.14 ×), ws66 (1.22–1.26 ×).
- Fail:
  - hcyl: 1.97 × reference (E ≈ 4 % below J outside the cylinder).
  - cb_r4 (4-ray ring): 2.97 ×, beams 15 cells wide (the consistency gate is off at beam edges).
- Not run: core_u (needs a Cartesian deck), cbh.

**Final form of item 4 (commit f68856a17): consistent half-range face factors.**
- h+_f = H+_f/J_L and h−_f = H−_f/J_R.
- H±_f are the rays' half-range fluxes reconstructed to the face (MUSCL, minmod; `clean_hface = muscl`), clamped to |h| ≤ 1.
- The face flux c(h+_f E_L + h−_f E_R) equals the rays' face flux when E = J.
- h+_f ≥ 0 ≥ h−_f keeps the M-matrix, so it is positive by construction: no gate, no limiter.
- The additive correction is superseded: it needed the gate (negatives without it), and its FCT limit starves at large CFL. It
  stays behind `clean_hr = additive` / `clean_gate` for comparison only.
- Status table: /viper/ptmp2/jinma/clean_1114/PROGRESS.md.
