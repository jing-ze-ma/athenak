# Implicit transport for `<rad_m1>` (stage 3): design note

Status: DESIGN (2026-09-21).  Companion of `rad_m1_design.md` (stages 1-2, explicit
transport); read its sections 1, 3, 4, 10-12 first.  Survey behind the choices below:
memory note `rad-m1-implicit-gpu-survey` (summarised in `docs/handover/HANDOVER-2026-09-21.md`).

```
<rad_m1>/transport = explicit        # stages 1-2: PD-ARS, sub-cycled, light-speed CFL
                   = implicit_x1     # stage 3a (this note): implicit along x1 columns
                   = implicit        # stage 3b: 3-D, Krylov + line preconditioner
```

**Why.**  Explicit transport needs `dt = 0.4 dx/chat`.  The reduced speed of light is valid
only for `chat >> v tau_max`: fine for the He box cut at `tau = 100` (K = 10: 9 substeps),
useless for the B star (FeCZ at `tau` 579-3268) and the global He4 model, and even in the
He box it inflates the radiation heat capacity by `c/chat` (error ~3/K) and is a suspect
for the growing 72 s oscillation of the 1-D column (stage 2a).  The implicit scheme runs at
the TRUE `c` with `dt = dt_hydro`, i.e. radiation CFL `c dt/dx ~ 1e3-1e4`.

**Scope of 3a**: the column (x1) solve, complete for 1-D problems, and later the line
preconditioner of 3b.  Out of scope here: the Krylov solver, multigrid, AMR level
boundaries, curvilinear geometry.

---

## 1. Formulation: eliminate the flux at the faces

What the literature supports (Enzo FLD, Reynolds et al. 2009: Schur complement to a scalar
radiation equation, MG-preconditioned CG, flat iteration counts to 4096 cores; Olivier et
al. 2024: eliminate the first moment, the Schur complement is a diffusion-like operator;
DSA theory: tau-independent iteration counts when the preconditioner matches the
transport discretisation) and what ARK-RT (Bloch et al. 2021, the only implicit M1 on
Kokkos) warns against (an assembled collocated `(E,F)` system with AMG: slower on GPU
than on CPU).  So: a SCALAR operator for `E`, matrix-free, with `F` living on faces.

Backward Euler over `dt`, true `c`, lab-frame moments, comoving opacities, closure
`P = chi(f) E` lagged (Picard).  The normal flux is a FACE variable.  Split as in the
explicit scheme, `F = F0 + A`, `A = v E + v.P` (enthalpy flux, upwinded with the face
velocity).  The `F0` equation at face `i+1/2`, with the stiff comoving source implicit:

```
F0_f' = theta_f [ F0_f^n - c^2 dt (chi_{i+1} E'_{i+1} - chi_i E'_i)/dx ]
theta_f = 1 / (1 + c dt (rho kappa_F)_f)          (rho kappa_F)_f = arithmetic face mean
```

Substituted into the energy equation this gives one linear equation per cell,

```
E'_i + (dt/dx) [ (A' + F0')_{i+1/2} - (A' + F0')_{i-1/2} ] = E^n_i + dt c G_E,i'
D_f = theta_f c^2 dt / dx      ->  c /((rho kappa_F)_f dx)   for c dt rho kappa >> 1
```

i.e. tridiagonal in `E'` along x1, with off-diagonals `-(dt/dx) D_f chi_{i+-1}` (negative)
plus the upwind advection: an M-matrix, so `E' > 0` for any `dt` (the collocated HLL block
system has no such guarantee).  Properties that follow without any `alpha` device:

* thick limit: `D_f chi -> c/(3 rho kappa dx)`, the compact 3-point diffusion operator,
  exact face diffusivity (harmonic mean of `D`) at an opacity jump: the T3b weakness of
  the explicit scheme (a cell-centred `F` divided by the cell's own opacity, leaking into
  `alpha avg(F)`) cannot occur, since no cell-centred `F` enters the transport;
* thin limit: `theta -> 1`, a backward-Euler discretisation of the M1 wave system on a
  staggered grid: strongly damped for `c dt >> dx`, which is the intended behaviour: each
  step returns the quasi-static radiation field, and the STEADY state is the exact steady
  M1 solution (`div F = sources`, `F0 = -c div P/(rho kappa)`);
* no odd-even mode (staggered), no Riemann dissipation to correct.

It is first order in time (as Jiang 2021 and AREPO-IDORT).  Second order (BDF2 /
theta-scheme) is a later option; the radiation field is quasi-static on `dt_hydro`.

**State.**  In implicit mode the persistent variables are `E` (cell) and the face-normal
fluxes `F0_f` (x1 faces in 3a).  The cell-centred `F_i` of `u0` is DERIVED after the solve
(mean of the two faces, plus `A`) for the closure (`f`, `n`), outputs, the explicit-mode
BC code and diagnostics.  Restart must carry the face fluxes (T8 gate).

## 2. Matter coupling inside the solve

Per cell, the energy exchange `G_E = rho (kappa_P a T'^4 - kappa_E E0')` is nonlinear in
`T'`.  Newton on the local gas equation, eliminated into the `E` system (the Enzo /
AREPO-IDORT route, and what makes a tabulated EOS work): with `e(T)`, `c_v` from
`EOS_Data::ThermoAt`,

```
rho e(T') - rho e^n = -dt c rho (kappa_P a T'^4 - kappa_E E')            (local)
linearise about the iterate T_k:  a T'^4 ~ a T_k^4 + 4 a T_k^3 (T' - T_k)
=> T' - T_k = [ r_k + dt c rho kappa_E (E' - E_k) ] / ( rho c_v + 4 dt c rho kappa_P a T_k^3 )
```

which adds a positive term to the diagonal and a source to the right-hand side (the
M-matrix property is kept).  Outer loop: solve the tridiagonal system for `E'`, update
`T` by the SAFEGUARDED scalar root find of the explicit scheme (rtsafe; on the table `c_v`
is not the slope of `e(T)`), re-evaluate `chi`, `theta`, opacities (optional), upwind
directions; stop at `max(|dE|/E, |dT|/T) < tol` (default 1e-8, max 30).  The gas energy is
then SET from the converged exchange so that `E_gas + E_rad` changes only by the boundary
fluxes and the work term (algebraic conservation, as in stage 1).  Velocity-dependent
pieces of `E0` and the work term `vbar.delta(rho v)` stay explicit (non-stiff).  Momentum:
`delta(rho v) = dt rho kappa_F F0'/c` from the face fluxes (half of each face to the two
cells: conservative), minus `dt rho a_rad_ref` under `force_reference = wb_arad`.

Placement: ONE implicit solve per hydro step in the operator-split lists (no
sub-cycling, no PD-ARS); `<rad_m1>/strang` halves optional later.

## 3. Boundaries

Face-flux form makes them simple: imposed flux (bottom): `F0_f = F_in`; free surface
(top): Marshak `F_f = c q E_top` with `q` from the closure at the free-streaming-outgoing
state (`q = 1/2` for `f = 1/2`; expose `<rad_m1>/marshak_q`); reflecting: `F_f = 0`;
periodic: cyclic tridiagonal (Sherman-Morrison).  Across MeshBlocks/ranks the line is
partitioned with the existing column-partition machinery
(`src/utils/two_stream_column_partition.hpp`) or, in 3a, restricted to one block per
column with a fatal otherwise (say which in the report).

## 4. 3-D (stage 3b, outline only)

Operator `A(E) = E - dt div[ theta (c^2 dt grad(chi E)) ] + dt div A_upw + coupling`,
matrix-free (one ghost exchange per application), off-diagonal Eddington-tensor terms
lagged in the Picard loop; BiCGStab (non-symmetric through advection and `chi`),
preconditioned by the x1 line solve of 3a, optionally ADI transverse lines
(`conduction_transverse.cpp`).  Known weak point: no coarse space, iteration counts grow
with the MeshBlock count (ARK-RT: doubled at >= 4 ranks); remedy = the ported multigrid
(`src/multigrid`, 7861549c) extended with a variable-coefficient stencil and an x1 line
smoother, used as the preconditioner.  Global reductions are the GPU cost (2-3 per
BiCGStab iteration); Chebyshev smoothing needs none.

## 5. Test plan for 3a

Reuse the stage 1 gates (`tests_m1/`), at radiation CFL `c dt/dx` = 1, 10, 1e2, 1e4:

**(I1) thick pulse** `tau_cell` = 10, 1e3, 1e6: diffusion rate within 2 % at CFL <= 10;
error vs `dt` first order; **Nyquist** mode: damped, no growth at any CFL.
**(I2) T3b opacity jump**, cases A and B: absolute flux error (explicit: 1.40 % / 0.33 %
low); target < 0.1 %.  **T3c top-hat**: <= 5 %.
**(I3) equilibration** at `dt` up to 1e5 thermal times, ideal gas and tabulated He EOS:
equilibrium to 1e-10, invariant to round-off, Picard iterations reported.
**(I4) Marshak wave** (constant `c_v`): L1 vs the S_N reference at CFL 1, 10, 100;
iterations per step.
**(I5) moving medium**: T4 dynamic (`beta tau` = 14): centre within 1 cell, width 2 %;
T4b: drift < 1e-12.
**(I6) free streaming**: beam/pulse front at CFL 1 and 100: report the damping honestly
(this scheme is not meant to propagate fronts at large CFL); steady grey atmosphere
against the semi-analytic M1 solution of `bench/m1_stage2/ic/build_ic.py` (discrete
steady state reached in O(1) steps at CFL 1e4).
**(I7) consistency**: implicit at CFL 0.4 vs the explicit scheme on T3 and T6 (agree to
truncation error).  Conservation over 200 steps; restart bitwise; 1 vs 4 MeshBlocks.
**(I8) the 1-D He column at true c**, no RSLA, no sub-cycling: does the 72 s oscillation
grow (rate vs the explicit K-scan of the stage 2a diagnosis)?  Cost per hydro step vs
explicit K = 10 (13 substeps) and vs the two-stream.
**Bitwise**: `transport = explicit` (default) reproduces `tests_m1/runs_1cB/RESULTS.txt`
and box G1 to the last digit.

Open risks: (R1) Picard convergence when opacities are updated inside the loop at large
`dt` (kappa-T feedback at the Fe bump); (R2) accuracy of backward Euler for the radiation
response on the pulsation time scale (phase error ~ `dt/P`, harmless for `P` = 72 s and
`dt` = 0.16 s); (R3) the advected enthalpy flux inside the implicit operator at
`beta tau >> 1` (upwind keeps the M-matrix; accuracy first order); (R4) the derived cell
`F` near `f -> 1` (closure direction from face means).

---

## 6. Findings of 3a (implemented 2026-09-21; supersede the text above where they differ)

Branch `m1-implicit-x1`.  Implementation: `src/rad_m1/rad_m1_implicit.{hpp,cpp}` plus
delimited hooks in `rad_m1.{hpp,cpp}`, `rad_m1_tasks.cpp`, `rad_m1_newdt.cpp`,
`driver.cpp`, `restart.cpp`/`pgen.cpp`/`outputs.hpp` (the face array in the restart
file) and `src/CMakeLists.txt`.  Gate tables and the exact commands:
`tests_m1/runs_3a/RESULTS.txt` and `runs_3a/run_gates.sh`.

**The formulation holds.**  The tridiagonal M-matrix, the face elimination and the
local Newton elimination of the emission term all behave as sect. 1-2 says.  `E` stayed
positive in every run; no floor was ever active on the solution; the Picard loop
converged in 1-4 iterations on everything but the Marshak front (6-14) and the He column
(13.5).  The predicted thick-limit property is measured: at an opacity jump of 1000 the
absolute flux error is **1.9e-3 %**, against 1.3 % for the explicit scheme (700x), and
the peak deviation at the jump falls from 33 to 5.5e-6.  The diffusion rate is within
0.1 % of `2D` at `tau_cell` = 10, 1e3, 1e6 at every CFL from 0.4 to 1e4, and the Nyquist
mode is damped at **exactly** the backward-Euler rate `ln(1 + lam dt)/(lam dt)` (five
digits, at CFL 0.4 to 100): first order in time, confirmed, with no growth anywhere.

**Three things the note did not say, that the implementation had to decide.**

1. *The gas energy must be set from the source the ASSEMBLED ROW applied*, `q = SRCR -
   SRCB E'`, not from `rho kappa_P a T'^4`.  The two differ by the Picard remainder of
   the linearisation, and taking the second leaves `e_gas + (c/chat) E` drifting
   2.8e-11 / 1.4e-9 / 1.3e-7 over 2000 steps of T5 at `dt` = 1e3 / 1e5 / 1e7 coupling
   times.  With the first the drift is 2e-13 and the total energy of the dynamic T4 run
   is conserved to 2.0e-16, better than the explicit scheme's 1.8e-13.
2. *A Dirichlet end cell* (`implicit_bc_x1min|max = efix`, which freezes `E` in the
   boundary cell) had to be added.  A pure-scattering column with an imposed flux on
   both ends leaves the operator SINGULAR -- there is no source to anchor the level of
   `E` -- so T3b cannot be run in the four boundary types of sect. 3 alone.  `efix` is
   the face-flux-form equivalent of the fixed-`E` ghost the explicit tests use.
3. *The Marshak boundary needs the incident bath*, `F_f = +-c q (E - E_bath)`
   (`implicit_ebath_x1min|max`, default 0 = the plain free surface of sect. 3).

**Momentum.**  Each x1 face hands `dt (rho k_t)_f F0'_f/c` to the gas, half to each of
its two cells and all of it to the single interior cell at a physical boundary.  That is
exactly what the implicit face source removed from the radiation, so the *exchanged*
momentum is conservative to round-off.  The cell-centred `F` that `u0` carries also holds
the advected enthalpy flux `A`, which is not an exchange, so the diagnostic invariant
`sum(rho v + F/(chat c))` drifts 2.4e-8 on T4 (explicit: 1.1e-13).  Using the CELL
opacity with the mean of the two faces would be more natural for a cell-centred momentum
but is not conservative; the face rule was chosen for that reason.

**The two weaknesses, both foreseen and both real.**

* *Free streaming* (risk R4 was about the derived `F`; this is worse than that).  A
  free-streaming pulse at `f = 1` loses 63 % of its amplitude in one box crossing at
  CFL 0.4 and 92 % at CFL 100, against 5 % for the explicit scheme.  The operator is a
  backward-Euler, non-upwinded discretisation of the wave system, and it damps the thin
  limit hard *even below CFL 1*.  `implicit_x1` is therefore not a drop-in replacement:
  it is for optically thick columns, and `transport = explicit` stays the default.
* *The advected enthalpy flux* is donor-cell upwinded inside the operator (risk R3).  On
  T4 dynamic the pulse width comes out 1.9 % too wide against the explicit 1.0 %,
  independently of `dt`.  Reconstructing `A` at the face while keeping the M-matrix
  (a flux-limited correction in the RHS, lagged in the Picard loop) is the obvious
  3a-follow-up.
* *Risk R1 is visible*: the Picard count grows with `dt` at a strongly nonlinear
  emission front (Marshak with `a_rad` = 1e30): 6.6 -> 9.9 -> 14.1 iterations at CFL
  1 -> 10 -> 100, and 8 of 259 solves hit `implicit_maxit` = 30 at CFL 100.  It never
  hit the cap on the stellar column (13.5 mean, 23 max).

**I8, the result the milestone exists for.**  The 1-D He FeCZ column at the true `c`,
one solve per hydro step, bottom = imposed flux, top = Marshak: **the oscillation does
not grow**.  The envelope decays (`gamma` = -2.9e-5/s) where the explicit `K = 10` arm
e-folds in ~195 s and reaches 239 `v_MLT`; the drifts at 975 s are smaller than at
400 s; `F1/F_in` stays inside 0.9968-1.0012 over all 84 cells where the explicit arm
spans 0.77-1.40; the residual force is 5.1e-3 `g0` against 4.2e-1.  The residual
`max|v1|` = 10.3 `v_MLT` is steady and sits in the BOTTOM cell, at the imposed-flux
boundary.  Cost: 4.74e-3 s/step against 9.51e-3 s/step for explicit `K = 10` (and
4.43e-4 s/step hydro-only), i.e. the implicit scheme at the true `c` is **twice as fast
as the reduced-speed-of-light explicit one**.  The stage-2a suspicion that the RSLA
drove the growing oscillation is supported.

**Gaps.**  (i) The semi-analytic grey atmosphere of `bench/m1_stage2/ic/build_ic.py` was
not run: there is no test problem generator for a static plane-parallel column with
fixed gas, and I8 covers the same physics with the real column.  (ii) The partitioned
line solve was not implemented -- more than one MeshBlock along x1 is a startup fatal.
(iii) No MPI run; the solve is rank-local by construction and the Picard count is
`MPI_MAX`-reduced so that ranks stay in step, but that path is untested.
(iv) `F_2 = F_3 = 0` always in this mode, and `implicit_allow_multid` gives independent
columns (verified bitwise for 1 vs 4 MeshBlocks along x2).  (v) No GPU run.

---

## 7. Findings of 3a2 (implemented 2026-09-21; supersede sect. 6 where they differ)

Branch `m1-implicit-2`.  Milestone 3a2 attacks the five limits listed at the end of
sect. 6.  Gate tables: `tests_m1/runs_3a2/RESULTS.txt`.  New options, ALL defaulting to
the 3a behaviour so that `runs_3a/RESULTS.txt` is reproduced by an input file that does
not name them:

```
<rad_m1>/implicit_flux     = central | ap_hll | berthon      (LIMIT 1)
         implicit_recon    = dc | plm_dc                      (LIMIT 1)
         implicit_recon_w  = <= 0 (auto 1/(1+chat dt/dx)) | w  (LIMIT 1)
         implicit_recon_lag= picard | step                     (LIMIT 1)
         implicit_res_floor= 0 | x  (Picard residual scale)     (LIMIT 1)
         implicit_bmom_half= false | true                       (LIMIT 3)
         implicit_partition= none | gather                      (LIMIT 4, see gaps)
<problem>/m1_test = atmosphere                                  (LIMIT 2)
```

### LIMIT 1 -- the thin limit.  Three corrections the note did not anticipate.

The face flux can now be written `G_f = A_up + alpha F_HLL + (1-alpha) F_diff`, linear
in `E'` because the reduced flux `f` and the wave speeds `b_L <= 0 <= b_R` are LAGGED.
With `F_L = c f_L E'_L`, `F_R = c f_R E'_R`,

```
G_HLL = [b_R chat f_L E'_L - b_L chat f_R E'_R + b_R b_L (E'_R - E'_L)]/(b_R - b_L)
      = c_L E'_L + c_R E'_R,
c_L = b_R (chat f_L - b_L)/(b_R - b_L) >= 0,   c_R = b_L (b_R - chat f_R)/(b_R - b_L) <= 0
```

both signs following from the HLL consistency condition `b_L <= chat f <= b_R`, which
holds on the M1 admissible set.  **The M-matrix argument is by COLUMNS, not rows.**  Each
face contributes `+nu c_L` to the diagonal of its left cell and `-nu c_L` to the lower
off-diagonal of its right cell, and `+nu c_R` to the upper off-diagonal of the left cell
and `-nu c_R` to the diagonal of the right cell; every column of the matrix therefore sums
to exactly `1 + SRCB >= 1` while all off-diagonals are `<= 0`.  That is strict column
diagonal dominance, hence a nonsingular M-matrix and `E' > 0` at any `dt`, for the HLL
part, for the diffusion part (same argument) and for any convex blend of them.  Row
diagonal dominance does NOT hold and is not needed.  The code carries `fmax`/`fmin` guards
on `c_L`, `c_R` anyway, because `alpha < 1` rescales the advective and dissipative parts
differently and can turn `c_R` positive; the clamp only makes the flux more upwind and
leaves conservation exact (one number per face, used with opposite signs).

*(a) The literal blend counts the diffusion TWICE.*  At the piecewise-constant states the
MATRIX is built from, the HLL dissipation `b_R b_L (E_R - E_L)/(b_R - b_L) -> -c dE/3` IS
the physical diffusion -- that is what `alpha` is constructed to make true -- so adding
`(1 - alpha) F_diff` on top doubles it: gate I1 came out at **2.0021 x** the analytic
`2D` at `tau_cell` = 1e3, at every CFL.  Weighting the dissipation `alpha^2` (the explicit
scheme's `alpha2` form) repairs `tau_cell` >= 1e3 (1.0010) but leaves **+10.7 %** at
`tau_cell` = 10 and **destroys T3b** (97 % face-flux error at the opacity jump).  What
works is Berthon's form, `implicit_flux = berthon`: `alpha F_HLL` with NO `F_diff` at all,
both parts of the HLL flux weighted by `alpha`.  The advective part then supplies the
`1/(0.866 tau)` the dissipation alone is short of, and I1 comes out 0.984 / 0.9998 /
1.0000 at `tau_cell` = 10 / 1e3 / 1e6 at every CFL.  `ap_hll` is kept as the literal
blend of the design brief, with its numbers, and is NOT recommended.

*(b) The lagged reduced flux must not be the cell mean of the two faces* (design risk R4,
worse than the note supposed).  In free streaming the upwind face flux is `c E_{i-1}`, so
the cell mean is `c (E_{i-1}+E_i)/2` and the derived `f` is `(1 + E_{i-1}/E_i)/2`, i.e.
0.5 rather than 1 on the steep side of a pulse.  The wave speeds reopen to `+-c/sqrt(3)`,
the HLL flux turns CENTRED and the I6 pulse is flattened to its box mean in one crossing
(amplitude ratio **0.0014**).  Each FACE flux is normalised instead by the `E` of the cell
it comes FROM, and the cell value is the mean of the two ratios: exactly 1 for an upwind
free-streaming face.  That single change takes I6 from 0.0014 to 0.561.

*(c) The plm deferred correction needs an admissibility clamp AND a weight.*  The
correction is the explicit difference between the plm and the dc flux at the previous
iterate, so the matrix stays the low-order M-matrix.  Two failures had to be fixed:
the reconstructed face energy can exceed the DONOR cell's own (`plm` puts
`E_i (r-1)/(r+1)` on top of `E_i` for a geometric ratio `r`), and `c` times that is faster
than the donor can emit, so on an exponentially falling background the drain cascades --
the I6 background collapsed onto the floor and the peak grew 7x.  Clamping the corrected
flux to `[-c E_R, +c E_L]` fixes it.  Second, a deferred correction is a fixed-point
iteration with contraction factor `~ 2 nu/(1 + nu)`, `nu = chat dt/dx`: it DIVERGES above
CFL ~ 1 (measured, I6 amplitude 3.85 at `implicit_cfl` = 10).  Weighting the correction by
`w = 1/(1 + nu)` (`implicit_recon_w <= 0`, the default) makes the factor `2 nu/(1+nu)^2
<= 1/2` at every CFL and the fixed point a convex blend of the dc and plm fluxes -- still
monotone, second order where `w -> 1`.

**I6, free-streaming pulse, 128 cells, one box crossing** (peak-to-trough amplitude,
final/initial; explicit plm = 0.953):

```
 implicit_cfl                 0.4      1        10       100
 central (3a)                 0.369    0.319    0.285    0.076
 berthon / ap_hll, dc         0.561    0.492    0.145    0.079
 berthon / ap_hll, plm_dc     0.818    0.695    0.190    0.079
```

so the 3a figure of 0.37 at CFL 0.4 becomes **0.82**, against 0.95 for the explicit
scheme.  `berthon` and `ap_hll` are identical here because `alpha = 1` with no opacity.
`plm_dc` costs `implicit_maxit` Picard passes against 2 (see "Picard" below).

### LIMIT 2 -- steady-state accuracy in a thin layer.  The result reverses the brief.

`<problem>/m1_test = atmosphere` (`src/pgen/tests/rad_m1_tests2.cpp`,
`inputs/tests/rad_m1_atmosphere.athinput`, gate `tests_m1/t9_atmosphere.py`): a static
plane-parallel pure-scattering column, prescribed gas, `rho = 0.128 exp[(1-z)/0.113]`,
`kappa = 1`, 128 cells, imposed flux at the bottom face and a free surface at the top.
`tau` at the top cell centre is 5.0e-4, `tau_bot` = 100.8.  The steady state is the
closed-form M1 solution `chi(f)/f = q0 + tau`, `E = F/(c f)`.

**The consistent pair.**  The Marshak condition `F = c q E_top` fixes the surface reduced
flux to `f(0) = q = <rad_m1>/marshak_q` EXACTLY, so `q0 = chi(q)/q`; nothing else in the
problem determines it.  For `q = 1/2`, `chi(1/2) = 0.46481586` and **`q0 = 0.92963172`**.
`chi(f)/f` has a minimum 0.89806 at `f` = 0.82504 and the solution stays on the decreasing
branch for any `q0` above it.

```
 max|E/E_exact - 1| over the column / over the top 5 cells      L1(E)
 implicit central, CFL 100 and 1e4      1.35e-3 / 1.35e-3       5.0e-4
 implicit berthon, CFL 100 and 1e4      3.68e-1 / 1.21e-1       2.2e-1
 implicit ap_hll,  CFL 100              3.53e-1 / 5.50e-2       1.9e-1
 explicit plm (tau_bot ~ 10 column)     4.86e-1 / 4.86e-1       8.4e-2
 implicit central (same tau_bot ~ 10)   1.35e-3 / 1.35e-3       5.1e-4
 steps to |dE/E| < 1e-10 (central): CFL 1e4: 16-32   CFL 100: 1024-2048
                                    CFL 1: > 8192 (still 1.9e-3 per doubling)
```

**The 3a `central` flux passes the 1 % target with a factor 7 to spare, and every upwind
form fails it by a factor 40.**  The reason is exactly the reason `central` damps a
propagating front: the face-eliminated form has NO numerical dissipation, while the
upwind HLL flux carries `0.289 c dx`, which exceeds the physical `c/(3 rho kappa)`
wherever `tau_cell < 1.15` and then sets the gradient of `E` itself.  Below `tau ~ 20`
the berthon column is systematically too bright, reaching +37 % at `tau` = 0.3.
`plm_dc` does not rescue it (the automatic weight is 1/101 at CFL 100; forcing `w = 1`
diverges).

**LIMIT 1 and LIMIT 2 therefore pull in opposite directions**, and for the stellar
columns this milestone exists for, LIMIT 2 is the one that matters: `implicit_flux =
central` stays the recommended default, and `berthon` is for problems where a thin
front must PROPAGATE.  The explicit scheme's 49 % surface error on the same problem
(its dark ghost lets `f -> 1` at the free surface instead of imposing `f = 1/2`) is a
separate point in favour of the implicit face-flux boundary.

### Other gates with the new fluxes

```
 I2  T3b opacity jump x1000, max|E/E_ex-1| / face-flux error at the jump
     explicit plm   1.05e-1 / 33         central (3a)  1.13e-4 / 5.6e-6
     berthon        5.88e-3 / 2.0e-3     ap_hll        8.2e-3  / 9.8e-1  FAIL
 I4  Marshak wave, 128 cells, L1 vs the S_N reference (explicit 0.0019)
     central dc    0.0144 (CFL 1) 0.0135 (CFL 100)      <- 3a, reproduced
     berthon dc    0.0104          0.0098
     berthon plm_dc 0.2647 FAIL    0.0089
     The 0.5 % target of the brief is NOT met by any implicit variant; berthon dc is
     the best at 1.0 %, a 28 % improvement on 3a.  plm_dc is unusable at a sharp
     emission front at CFL 1.
 STEP 0  transport = explicit: T3 tau_cell 1e3 plm ratio 1.000015, identical to
     runs_1cB and runs_3a.  The explicit path is untouched.
```

### Picard statistics with the new flux

`berthon`/`ap_hll` with `dc` converge in **2 passes** on every test that 3a needed 2-4
for; the lagged `f` and wave speeds cost nothing.  `plm_dc` never meets `implicit_tol`
= 1e-8: the plm limiter switches on a handful of cells and the iteration is a small
limit cycle.  The SOLUTION is converged -- `implicit_maxit` 30 and 100 give the same I6
amplitude to five digits -- but the strict test is never passed, so `plm_dc` costs
`implicit_maxit` passes per step.  Under-relaxing the correction by 1/2 across passes
(implemented, always on) does not break the cycle; nor does scaling the residual by the
column peak (`implicit_res_floor`, added for this and kept because it is the right test
for a column with a large dynamic range).  Evaluating the correction once per step
(`implicit_recon_lag = step`) makes the loop converge in 2 passes but is an EXPLICIT
anti-diffusion and blows up (I6 amplitude 7.3 at CFL 0.4); `picard` is the default.

### LIMIT 3 -- the residual flow at the imposed-flux bottom boundary.  Found and fixed.

It is candidate (a) of the brief.  A PHYSICAL BOUNDARY FACE hands its whole
`dt (rho k_t)_f F0_f/c` to its one interior cell, while every interior cell receives half
of each of its two faces; the bottom cell therefore feels 1.5 face-shares of radiative
force where every other cell feels 1.0, and the extra half share is a STEADY force that
the well-balanced reference `arad_ref` -- built from the cell-centred flux -- does not
carry.  `<rad_m1>/implicit_bmom_half = true` gives the boundary face HALF, like any other
face; the missing half leaves the domain with the radiation, which is where it goes, and
the momentum EXCHANGED inside the domain stays conservative to round-off.

On the 1-D He FeCZ column at the true `c`, 1000 s, against the 3a arm rerun on the same
binary (`v_MLT` = 1.86e4 cm/s):

```
                              3a            bmom_half = true
 |v1| in the BOTTOM cell      10.30 v_MLT   0.47 v_MLT          (22x)
 whole-column max |v1|        10.30         1.05, and it MOVES to the TOP cell
 |v1| in the TOP cell          4.05         0.96                (4.2x)
 |F1bot/F_in - 1|             9.3e-4        1.65e-4
 F1top, F1mid deviation       2.0e-5, 1.9e-5   2e-6, 1e-6       (no regression: better)
 residual force               5.13e-3 g0    2.00e-4 g0          (26x)
 column KE at 1000 s          4.11e28       3.36e26             (122x)
 cost / Picard passes         30.5 s, 13.46 26.2 s, 11.66
```

The remaining 1 `v_MLT` sits in the TOP cell, i.e. at the free surface, which the gate
allows; and it is still falling at 1000 s.  Candidates (b) and (c) are ruled out: the
imposed flux is carried to 1e-3 in BOTH arms, and the flow never propagated away from
`i = 0`.

### LIMITS 4 and 5 -- not done

*(4) The partitioned line solve is NOT implemented.*  More than one MeshBlock along x1 is
still a startup fatal.  `<rad_m1>/implicit_partition` is parsed (`none | gather`) so that
input files can already name it, and `gather` fatals with "NOT IMPLEMENTED".  The design
to implement: `two_stream_column_partition.hpp` does NOT fit -- it partitions one column
over a Kokkos THREAD TEAM inside one MeshBlock (Schur-style segment condensation into a
reduced block-tridiagonal system over the segment boundaries), not over MeshBlocks or
ranks.  For the scalar tridiagonal with Picard-lagged coefficients the simplest design
that is BITWISE independent of the partition is to gather each column's assembled
`(a,b,c,r)` rows onto the rank that owns its lowest block, run the IDENTICAL serial
Thomas sweep there, and scatter the solution back: the arithmetic order is then literally
unchanged, so 1, 2 and 4 blocks agree to the last bit by construction.  Cost: two
messages per Picard iteration per column group (a gather of `4 nx1` reals and a scatter
of `nx1`), and a serial bottleneck of `nx1` that is irrelevant at the 84-512 cells these
columns have.  A segment-condensation variant removes that bottleneck but cannot be
bitwise, and would have to be quantified (<= 1e-14 relative, non-accumulating).

*(5) No GPU run.*  Nothing new needs anything the 3a kernels did not: `ifw` is an ordinary
`DvceArray5D`, the new kernels are plain `par_for` over `(m,k,j,i)` capturing only Views
and scalars by value, `PLM()` is `KOKKOS_INLINE_FUNCTION`, and the one new host-side
reduction (the column peak of `E` for `implicit_res_floor`) has the same shape as the
Picard residual reduction 3a already ran on the device.  No known blocker; UNTESTED.

### Recommended defaults

```
 implicit_flux      central   KEEP.  It is the only form that passes the atmosphere gate,
                              which is what the stellar columns need.  Use `berthon` only
                              where a thin front must PROPAGATE (1.5x better than central
                              on I6 at CFL 0.4, 2.2x with plm_dc), never `ap_hll`.
 implicit_recon     dc        KEEP.  plm_dc helps only the free-streaming gate, costs
                              implicit_maxit Picard passes instead of 2, and fails the
                              Marshak gate at CFL 1 (L1 0.26 against 0.010).
 implicit_recon_w   -1 (auto) KEEP; a fixed 1.0 is divergent above CFL 1.
 implicit_recon_lag picard    KEEP; `step` is an explicit anti-diffusion and blows up.
 implicit_res_floor 0         CHANGE to ~1e-8 for any column with a large dynamic range:
                              the pure |dE|/E test is dominated by cells orders below the
                              peak.  Left at 0 so that 3a is reproduced.
 implicit_bmom_half false     **CHANGE to true.**  This is a bug fix, not a tuning knob;
                              it is left false only so that runs_3a/RESULTS.txt is
                              reproducible from the same source.
```

---

## 8. Findings of 3b (implemented 2026-09-21; supersede sect. 6-7 where they differ)

Milestone 3b was scoped as five phases (defaults + restart, the partitioned x1 line
solve, transverse transport, the He box in 2-D/3-D, GPU).  **Phases 0 and 1 are done and
gated; phases 2, 3 and 4 are NOT started.**  Gate tables and the exact commands:
`tests_m1/runs_3b/RESULTS.txt`.

### PHASE 0 -- defaults, and the restart gate 3a left open

`<rad_m1>/implicit_bmom_half` now DEFAULTS TO `true` (sect. 7, LIMIT 3: it is a bug fix,
the boundary face handing its whole momentum to one interior cell; He column bottom-cell
`|v1|` 10.3 -> 0.47 `v_MLT`).  The key is kept, so `false` still reproduces
`runs_3a/RESULTS.txt` from the same source, and
`tests_m1/runs_3a/he_box_m1_1d_impl.athinput` now names `true`.

**Restart is BITWISE.**  The 1-D He column (`implicit_x1`, true `c`, `central`/`dc`),
300 s, restarted from the restart file the uninterrupted run itself wrote at t = 150.09 s:
the hydro history, the six-column user history and the final `hydro_w` binary dump are
byte-identical over t = 150-300 s.  (The restart must be taken from a dump the
uninterrupted run wrote; stopping a separate run at `tlim` = 150 truncates its last step
and is a different trajectory -- that is a property of the driver, not of the solve.)
So the face-flux array `f0x1` carried by the restart file is the complete persistent
state of the implicit transport: nothing else has to be saved.

### PHASE 1 -- `<rad_m1>/implicit_partition = gather`, the line solve across MeshBlocks
and ranks (LIMIT 4 of sect. 6, closed)

What is implemented is exactly the design recorded in sect. 7: each column's assembled
rows `(a,b,c,r)` are gathered onto the block with the lowest x1 logical location, the
IDENTICAL serial Thomas recurrence is run there over `part_nblk*nx1` rows, and the
solution is scattered back.  Three things the design note did not say:

1. *The gather alone is not enough.*  The ASSEMBLY of a row next to a block face reads
   the lagged closure `w`, the enthalpy coefficient `a`, `g0`, `v1`, the comoving reduced
   flux and the transport opacity of the neighbouring cell, and the FACE update reads the
   new `E` there as well.  A per-Picard-pass x1 ghost exchange (`ImplicitX1Halo`) was
   added, of exactly those six lagged quantities after the lag kernel and of `E` after
   the line solve.  It moves the VERY NUMBERS the neighbour computed -- never a quantity
   recomputed from a hydro or opacity ghost, which would be bitwise-fragile -- and that is
   what makes the partition bit-exact.  The transport opacity was copied into a new
   component `M1_IW_KT` of the work array for the same reason: every quantity read at a
   neighbouring cell now lives in one array.
2. *The convergence test has to be global.*  3a reduced only the final iteration count.
   With a partitioned column the ranks must take the same number of passes or the gather
   deadlocks, and even with rank-local columns a per-rank test makes the answer depend on
   the decomposition.  One `MPI_Allreduce(MAX)` of one double per Picard pass was added
   (unconditional, so rank-local multi-column runs became decomposition-independent too).
3. *Periodic x1 across several MeshBlocks is a CLEAN FATAL.*  The cyclic (Sherman-Morrison)
   sweep of 3a wraps inside one block; a wrapped gathered column would need the halo to
   wrap as well.  One block along x1, or a non-periodic x1 pair.  SMR/AMR and a non-integer
   block count along x1 are fatal as before.

**Communication.**  Per non-root block and Picard pass: 2 gather/scatter messages
(`4*nk*nj*nx1` reals up, `nk*nj*nx1` down) plus up to 4 halo messages
(`nq*nlay*nk*nj` reals, `nq` = 6 or 1, `nlay` = min(nghost,2)), plus one `MPI_Allreduce`
of one double per pass for the whole communicator.  Same-rank stack members and same-rank
x1 neighbours are plain device copies and send nothing.

**Cost.**  The root sweeps `part_nblk*nx1` rows serially per column; the kernel is
parallel over (root, k, j), so this is free as long as there are many columns, and it is
the serial bottleneck when there are not.  Measured on this serial CPU: the 1-D He column
(84 cells, ONE column, the worst possible case for a serial sweep) costs 1.58 / 1.62 /
1.64 s of CPU for 1 / 2 / 4 blocks along x1, i.e. **+2.3 % and +3.5 %**; the Marshak wave
(128 cells, one column) 0.145 / 0.185 s, +27 % at 4 blocks.  The scalable successor,
which is NOT implemented and cannot be bitwise, is Schur condensation to the
block-interface unknowns (two local sweeps plus one reduced tridiagonal system of
`part_nblk` rows per column) or cyclic reduction.

**Gates -- all BITWISE, 1 vs 2 vs 4 MeshBlocks along x1 and 2 vs 4 MPI ranks against the
serial single-block run** (final `.tab`, all four moments, `maxabs` difference exactly 0,
and identical Picard statistics):

```
 thick pulse tau_cell 1e3, implicit_cfl 100, reflecting x1
    implicit_flux = central | berthon   x   implicit_recon = dc | plm_dc   (4 cases)
 Marshak wave,  implicit_cfl 100, marshak/marshak                 (Picard 14.14 mean)
 static grey atmosphere, implicit_cfl 100, flux/marshak           (Picard 1.18 mean)
 MPI: 2 and 4 ranks x 4 blocks along x1, berthon/plm_dc pulse, Marshak, atmosphere
```

`bitwise` here is literal and by construction, not "to the solver tolerance": the
arithmetic sequence of the Thomas sweep and of every kernel is unchanged by the
partition.

**The 1-D He column is NOT a valid decomposition gate.**  With 1 vs 2 MeshBlocks along x1
its history differs at 1.9e-5 (relative, total x1 momentum) -- but the SAME input with
`transport = explicit` differs at 3.2e-4, i.e. the `box_convection` problem generator
itself is not block-decomposition invariant (the well-balanced wall continuation
`bc_mode = 3` and the column diagnostics are per-block).  The implicit partition adds
nothing to that: the Picard statistics are identical to the last digit (373 solves, mean
11.76139, max 23) for 1, 2 and 4 blocks.  The six user-history columns `F1top/F1mid/F1bot`
are summed over blocks by the pgen's history hook and come out `nblk` times too large;
that is a pgen diagnostic bug to fix before any multi-block He run is judged by them.

### Regression

`tests_gate_merge/postmerge.sh`: box G1 modes 3 and 0, 10/10 files IDENTICAL; the 1-D He
column smoke test reproduces the new reference line exactly (t = 300, F1top/Fin 1.0000097,
F1bot/Fin 0.9998352, V1max 1.8344e4, Picard mean 11.709).  `transport = explicit`, T3
`tau_cell` = 1e3 plm: ratio **1.000015**, identical to runs_1cB, 3a and 3a2.  I1 at
`tau_cell` = 1e3: central 0.999999 (CFL 1) / 1.000001 (CFL 100), berthon 0.999825 /
0.999827 -- the 3a and 3a2 numbers.

### PHASES 2-5 -- not done

*(2) Transverse (x2/x3) implicit transport, `implicit_solver = line_jacobi | bicgstab`:
not started.*  `transport = implicit` (as opposed to `implicit_x1`) does not exist; a
multi-dimensional run is still a set of INDEPENDENT x1 columns under
`implicit_allow_multid`.  The x1 halo built for LIMIT 4 is the piece a line-Jacobi outer
iteration would reuse: it already shows that a hand-rolled per-iteration exchange driven
from inside the Picard loop works and stays bitwise, so the transverse couplings can be
lagged through the same mechanism (extended to x2/x3 neighbours) without restructuring
the solve into a task loop.  *(3) The He box in 2-D/3-D, (4) GPU: not started.*  The 3a2
assessment of the GPU risk is unchanged and now also covers the new kernels, except that
the gather/scatter stages `part_sbuf`/`part_rbuf` through HOST mirrors for MPI: on a GPU
that is a device-to-host round trip per Picard pass and would be the first thing to
replace with GPU-aware MPI (the arithmetic stays on the device either way, which is what
keeps the partition bitwise).

---

## 9. Findings of 3c (implemented 2026-09-21; supersede sect. 6-8 where they differ)

Milestone 3c answers the question sect. 7 left open: LIMIT 1 (`berthon`, good at fronts)
and LIMIT 2 (`central`, good at a diffuse steady field) pull in opposite directions --
can a SMOOTH per-face blend have both?  Gate tables and the exact commands:
`tests_m1/runs_3c/RESULTS.txt` and `runs_3c/run_gates_3c.py`.  New options, all inert
unless `implicit_flux = blend`:

```
<rad_m1>/implicit_flux        = central | ap_hll | berthon | blend
          implicit_blend      = tau | f | tau_f          (W1 | W2 | W3)
          implicit_blend_fmode= max | mean               (f_face from the two cells)
          implicit_blend_mode = flux | dissipation       (what is blended)
          implicit_blend_tau0 = 1.0    implicit_blend_flo = 0.6   implicit_blend_fhi = 0.9
          implicit_recon_npass= -1     (freeze the plm correction after N Picard passes)
```

**The form.**  `F_f = (1 - w_f) F_central + w_f F_berthon`, `w_f` in [0,1] built from the
PREVIOUS Picard iterate at the face, so the row stays linear.  It needed no new
machinery: the 3a2 assembly already scales the central (face-eliminated) contribution by
`om = 1 - ifw(AL)` and adds `ifw(HCL) E_L + ifw(HCR) E_R + ifw(DG)`, so the blend is
`AL = w_f` with the berthon coefficients multiplied by `w_f`.  `w_f = 1` is therefore
BITWISE `berthon` and `w_f = 0` BITWISE `central` (the code multiplies by exactly 1.0 or
0.0), and the face update carries the same weight, so the stored `F0_f` -- the restart
state, the momentum deposit and the next iterate's reduced flux -- is the flux the row
applied.  Three weights were built and measured: `tau`, `w = exp(-(tau_f/tau0)^2)` (the
brief's suggestion, after Jiang 2021); `f`, a smoothstep in the lagged comoving reduced
flux at the face between `f_lo` and `f_hi`; and `tau_f`, their product.

**The M-matrix survives, by the column argument of sect. 7.**  Each face contributes
`+nu c` to the diagonal of one cell and `-nu c` to an off-diagonal of the other, i.e.
each contribution has ZERO column sum and non-positive off-diagonals -- true of the HLL
part (sect. 7) and, with `c = th c^2 dt w_i/dx`, of the face-eliminated part.  A convex
combination with the SAME `w_f` in both cells has the same structure, so the column sums
of the whole matrix stay `1 + SRCB >= 1`: nonsingular M-matrix, `E' > 0` at any `dt`.
Verified numerically: the floor was never active on a blend solution in any gate (the
worst case, the free-streaming pulse with `plm_dc`, keeps `min E = 5.6e-6` against a
floor of `1e-15`).

**The result: an f-gated blend is free, and the tau-only weight is harmful.**  W2 and W3
are BITWISE `central` on every diffuse problem -- the grey atmosphere (1.353e-3 at CFL
1, 1e2, 1e4), the Marshak wave (L1 0.0144/0.0135), the T3b opacity jump (1.126e-4, face
flux peak 5.6e-6), the thick pulse at `tau_cell` 10/1e3/1e6, T4 dynamic and T4b, T5 with
both EOS, and the 1-D He column (`F1top/Fin` 1.0000097, `F1bot/Fin` 0.9998352, `V1max`
1.83439e4, Picard 11.709 -- the 3b reference line to the last digit) -- and BITWISE
`berthon` in free streaming, where every weight is 1.  **W1 fails exactly where the
brief predicted it would**: the top of a grey atmosphere is thin AND diffuse
(`f -> marshak_q = 1/2`), so a tau-only weight sends it upwind and the gate goes from
1.35e-3 to 3.7e-1; it also degrades the T3b face flux 157x and, on the He column,
multiplies the residual force by 245 and stops the Picard loop converging at all.

**Free streaming, the gain.**  With `plm_dc` the blend reaches **0.922** of the initial
pulse amplitude after one box crossing at CFL 0.4 (explicit 0.953, `berthon` + `plm_dc`
0.818, `central` 0.369), and holds **0.786 at CFL 10** where `berthon` + `plm_dc`
collapses to 0.190.  The reason is that the deferred correction carries `w_f` as well,
so it is switched off in the diffuse tails where the unweighted plm correction fights
the solve.  The same weighting repairs the two `plm_dc` defects of sect. 7: the Marshak
wave no longer fails at CFL 1 (L1 0.0144 against 0.2647) and costs 6.6 Picard passes
instead of 23.8, because `w_f = 0` at the emission front removes the correction there.

**What the blend does NOT do.**  It does not recover `berthon`'s 28 % Marshak advantage
(L1 0.0104): that advantage comes from faces where the field is DIFFUSE and only
moderately thick, which is precisely where the atmosphere gate forces the blend to stay
central.  Measured, not assumed: opening the window to `f_lo = 0.3, f_hi = 0.6` leaves
L1 at 0.0144, and the tau-only weight gains 0.0001.  **The two limits are reconcilable
in free streaming and not in a semi-thick emission front.**

**Two variants rejected.**  `implicit_blend_mode = dissipation` (the full central flux
plus `w_f` times the HLL dissipation alone) keeps the M-matrix but does not help the
front -- 0.324 at CFL 0.4, BELOW plain `central`'s 0.369 -- and destroys the Picard
convergence (17-28 passes against 2).  `implicit_recon_npass` (freeze the plm limiter
after N passes, to break the `plm_dc` limit cycle of sect. 7) makes the pulse worse
(0.444 against 0.922 at npass = 3) and does NOT restore convergence: the limit cycle is
driven by the re-lagged closure and wave speeds, not by the limiter.  Both keys are kept
(defaults `flux` and `-1`) with these numbers recorded.

**Partition.**  `implicit_partition = gather` with the blend is BITWISE for 1 vs 2 vs 4
MeshBlocks along x1, on the atmosphere and the Marshak wave, with both W1 (whose weight
varies strongly across block faces) and W3, and with identical Picard statistics.  No
new halo was needed: the weight reads only `M1_IW_RF0` and `M1_IW_KT` of the
neighbouring cell, and both are already in halo A of the 3b exchange.

### Recommended defaults after 3c

```
 implicit_flux   central   UNCHANGED.  The decision rule for promoting a blend was
                           "at least as good as central everywhere AND improves the
                           free-streaming gate AND the Marshak gate"; the Marshak gate
                           is NOT improved, so the default stays.
 blend           use `implicit_flux = blend` + `implicit_blend = tau_f` (or `f`) for any
                 problem with both diffuse regions and a propagating front: it is
                 bitwise `central` in the diffuse regions and bitwise `berthon` in free
                 streaming, i.e. it costs nothing and removes the reason `berthon`
                 could not be used.  Add `implicit_recon = plm_dc` where the front
                 matters (0.922 against 0.561, and without the 3a2 Marshak failure).
 implicit_blend  NEVER `tau`: it fails the atmosphere gate, the opacity jump and the
                 He column.  The weight has to know that the field is BEAMED, not only
                 that the cell is thin.
 implicit_blend_fmode max.  `mean` halves the free-streaming gain (0.498 against 0.922
                      at CFL 0.4) because at the EDGE of a pulse one of the two cells is
                      still diffuse; on the atmosphere the two are identical.
 implicit_blend_mode  flux (dissipation is rejected, see above)
 implicit_recon_npass -1 (rejected, see above)
```

---

## 10. Findings of 3b phases 2-4 (2026-09-21; supersede sect. 6-9 where they differ)

Phases 2-4 of milestone 3b were scoped as: (A) two `box_convection` fixes, (B) transverse
x2/x3 implicit transport with `implicit_solver = line_jacobi | bicgstab`, (C) the He box
in 2-D/3-D, (D) GPU.  **Only PHASE A is done and gated.  B, C and D are NOT STARTED** --
`<rad_m1>/transport = implicit` does not exist, and a multi-dimensional run is still the
set of independent x1 columns of sect. 8.  Gate tables and the exact commands:
`tests_m1/runs_3b2/RESULTS.txt`.

### PHASE A1 -- the six `<rad_m1>` user-history columns under decomposition

The defect sect. 8 recorded is fixed.  `F1top`/`F1mid`/`F1bot` and `V1mid` were read at
each MeshBlock's OWN `i = ie / is / is+nx1/2`.  Tiled in x2/x3 that was already a proper
plane mean (`inc` divides by the MESH's `nx2*nx3`); STACKED ALONG x1 every block
contributed its own end plane, so the three flux columns came out `part_nblk` times too
large and mixed heights.  A per-block table built from the block's logical x1 location now
holds the LOCAL `i` of the GLOBAL bottom / mid / top plane, or -1.  Maxima (`V1max`,
`Fres`) and volume sums (`Etot`, `KEcol`) were always right.  Measured on the 1-D He
column at t = 20 s: 1, 2 and 4 blocks along x1 now agree to 2e-11, where the HEAD code
gave **exactly 2x** on 2 blocks and a `V1mid` of the wrong sign.  On ONE block the
reference line is untouched (postmerge: `F1top/Fin` 1.0000097, `F1bot/Fin` 0.9998352,
`V1max` 1.8343926e4, Picard 11.70892, box G1 10/10 identical).

### PHASE A2 -- `box_convection` is not x1-decomposition invariant, and it is not a bug

**The cause is round-off in the MESH CELL-CENTRE COORDINATES, amplified ~1e7 by the
momentum cancellation.**  `x1v` is evaluated from each MeshBlock's own
`(x1min, x1max, nx1)`, so the same physical cell centre comes out ~2 ulp apart under a
different decomposition (measured: 8.2e-14 relative).  The pgen interpolates the initial
column, the M1 IC and `a_rad_ref` at those coordinates, so the INITIAL STATE already
differs at 5e-15 -- before any step -- and after one cycle `mom1` differs at 1.8e-6
relative in EVERY cell, not at the block interface.  The control settles it: a
single-block run with `mesh/x1min` moved by ONE ULP reproduces the same signature at the
proportionally smaller size (coordinate 15x smaller, momentum 35x smaller).  The total
x1 momentum of a hydrostatic column starts at exactly zero and remains ~1e-8 of the force
terms that cancel into it, which is where the 1e7 comes from.  Nothing in `wb_phi_eff`,
the `a_rad_ref` fill, the M1 IC reader, the `bc_mode = 3` wall walk or the BC branch
treats a MeshBlock's x1 extent as the whole column; the implicit `gather` partition adds
nothing (identical Picard statistics for 1, 2, 4 blocks, and `transport = explicit` shows
the same divergence).

**So the "1 vs 2 vs 4 blocks bitwise" gate is not achievable for this problem and never
was.**  What holds: mass and total energy to round-off (1.4e-15 / 4e-16 over 60 s); the
total x1 momentum to 1.9e-5 relative, seeded by 5e-15.  The only way to make it bitwise
is to compute cell centres from the GLOBAL mesh extents and the global index instead of
from each block's edges -- a mesh-wide change, not made here, and the recommendation to
carry into any milestone that needs decomposition-invariant stratified columns.

Reported and NOT touched (a two-stream-only path): the non-M1 branch of `BoxConvHistory`
reads the emergent flux at the block's own `ie+1` face and its comment states the
assumption that one MeshBlock spans the whole x1 extent; it has the same defect A1 had.

### What phase B still needs

Unchanged from sect. 8: the x1 halo of LIMIT 4 is the mechanism to extend to x2/x3
neighbours (a per-iteration exchange driven from inside the Picard loop, deterministic),
the persistent face-normal fluxes need the same restart treatment as `f0x1` with a
backward-compatible read, and the linear-system residual of the full 7-point operator has
to be tested separately from the Picard residual or lagging the transverse couplings will
look converged when it has only stalled.

---

## 11. Findings of 3b phase B (2026-09-21; supersede sect. 4 and 8-10 where they differ)

Phase B is the TRANSVERSE implicit transport sect. 8 and sect. 10 left open:
`<rad_m1>/transport = implicit`, the same face-eliminated backward-Euler solve as 3a with
the x2 and x3 couplings added.  **Implemented and gated with
`implicit_solver = line_jacobi`; `bicgstab` is a clean fatal.**  Gate tables and the
exact commands: `tests_m1/runs_3b3/RESULTS.txt` and `runs_3b3/run_gates_3b3.py`.

```
<rad_m1>/transport        = explicit | implicit_x1 | implicit
          implicit_solver = line_jacobi | bicgstab     (bicgstab: NOT IMPLEMENTED)
          implicit_lin_tol= 1e-10   the max-norm residual of the FULL 7-point system
          implicit_maxit  defaults to 200 under transport = implicit (30 otherwise)
<problem>/pulse_1d, pulse_y0, pulse_z0   the 2-D/3-D thick pulse of the gate
```

### The operator

Per direction `d` and face `f`, exactly the x1 form of sect. 1, with the diagonal
Eddington component of that direction and the off-diagonal ones lagged:

```
F0_f' = th_f [ F0_f^n - c^2 dt (D_dd,R E'_R - D_dd,L E'_L)/dx_d - c dt v_f g0_f
               - c^2 dt (sum_{e != d} d_e P_de)_f ]
th_f  = 1/(1 + c dt (rho kappa_t)_f)        (rho kappa_t)_f arithmetic face mean
D_ab  = (1-chi)/2 delta_ab + (3 chi - 1)/2 n_a n_b     (lagged chi and flux direction n)
A_d   = v_d E + (v.P)_d = a_d E,  a_d = v_d + (v.D)_d,  upwinded with the face velocity
```

The off-diagonal divergence is evaluated with centred differences of the PREVIOUS pass'
`E` and closure (one-sided at a physical boundary) and is a pure right-hand-side term; it
enters the x1 face equation as well, so a multi-D x1 flux is no longer the 1-D one.

**The x1 couplings stay in the tridiagonal matrix** and are solved exactly per `(k,j)`
column, with the 3b `gather` partition when blocks are stacked along x1.  The x2/x3
couplings are LINE JACOBI: the neighbours' `E` is the previous pass', but their DIAGONAL
contribution `dT/dE_c >= 0` is kept on the matrix diagonal.  That is what preserves the
full 7-point M-matrix -- every face still contributes `+nu` to one diagonal and `-nu` to
one off-diagonal, so the column sums stay `>= 1` and `E' > 0` at any dt.  It is also why
the line-Jacobi iteration converges at all: what is lagged is strictly off-diagonal.

`transport = implicit` accepts `implicit_flux = central` only (fatal otherwise): the
asymptotic-preserving forms of 3a2/3c carry per-face coefficients (`ifw`) built for the x1
faces, and silently running them as `central` in x2/x3 would be a trap.

### Halo: the module's ORDINARY cell-centred exchange on a scratch array

The hand-rolled `ImplicitX1Halo` of 3b was NOT extended.  Instead the 13 lagged
quantities the transverse faces and the lagged off-diagonal terms read at a neighbour
(`E`, `chi`, `n_1..n_3`, `rho kappa_t`, `v_1..v_3`, `a_1..a_3`, `g0`) are copied into a
scratch 5-D array `thw` and exchanged with `MeshBoundaryValuesCC` once per Picard pass,
then copied back into the ghost zones of `iw`.  That buys periodic wrap, edge and corner
neighbours and MPI with no new protocol, and it is still bit-exact in the sense that
matters: a block reads the VERY NUMBERS its neighbour computed, so the two blocks that
share a face build that face's flux from identical inputs and the face flux is
single-valued (gate G4: bitwise).  Two exchanges per pass (before the assembly, and after
the accepted iterate for the face update).  Under `transport = implicit` the x1 halo is
not run at all -- this exchange carries its six quantities and seven more.

PHYSICAL (non-periodic) x2/x3 boundaries are not filled by the exchange: every kernel
branches on the MeshBlock boundary flag and imposes `F = 0` there (reflecting).  Marshak /
imposed-flux transverse boundaries are NOT implemented.

### Convergence: the Picard test is not enough

`implicit_tol` on `max(|dE|/E, |dT|/T)` and, separately, the TRUE residual of the full
7-point system below `implicit_lin_tol`.  The residual is free: after the line solve
`E^{k+1}` satisfies `D1(E^{k+1}) + TDIA E^{k+1} + U(E^k) = b` with
`U(E) = T(E) - TDIA E_c`, so the residual of the FULL system at `E^{k+1}` is exactly
`U(E^k) - U(E^{k+1})` -- the change of the lagged off-diagonal term between two passes,
normalised by the max norm of the right-hand side.  It costs one extra pass (the first
pass has no previous `U`), which is why the `implicit` Picard counts are one above the
`implicit_x1` ones on an x2-uniform problem.

### State and restart

`f0x2` (and `f0x3` in 3-D) are allocated only when `transport = implicit` and the mesh is
multi-D, and are written to and read from the restart file immediately after `f0x1`, in
the same order in writer and reader (one shared lambda each, so the order cannot drift).
A file that does not carry them is accepted with one warning and the arrays start at zero;
that negotiation is unambiguous only when the general-EOS `wtemp`/`wder` tail is absent
(it is a single extra length, added to the same total).

### Gates (serial CPU, `build_cpu_m1`)

All numbers in `tests_m1/runs_3b3/RESULTS.txt`.  Summary:

* **G1 isotropy and rate.**  2-D 64^2 and 3-D 32^3 thick pulse, pure scattering,
  `gas_feedback = false`, `tau_cell` 10 and 1e3, `implicit_cfl` 1 / 1e2 / 1e4.
  `d(sigma^2)/dt / 2D` is 0.9955-1.0000 in every measurable case (tolerance 2 %) and the
  spread BETWEEN directions is 0 to 1.7e-9, i.e. three to six orders better than the
  1e-3 the brief asked for.  Three (tau, CFL) combinations are marked NOT MEASURABLE:
  the diffusion length of the twelve steps needed for a rate fit then exceeds the box.
* **G2 anisotropic cells.**  `dx2 = 4 dx1` (the He box ratio) and `dx2 = dx1/4`, same
  accuracy (0.999999 in both directions), isotropy 1.6e-7 / 2.8e-8, and the pass count
  moves from 3 to 4.
* **G3 a problem uniform in x2** (`nx2 = 4`, `problem/pulse_1d`): `transport = implicit`
  reproduces `transport = implicit_x1` to **9.9e-14** of the column scale (`E` to 4.3e-14
  absolute on a peak of 0.91, `F1` to 7.2e-18 on a peak of 7.3e-5), with `F2` identically
  zero.  NOT bitwise, and cannot be: the multi-D closure builds `D_11` from
  `(1-chi)/2 + (3chi-1)/2 n_1^2` where the 1-D one uses `chi` directly, and the linear
  residual test adds one Picard pass.
* **G4 decomposition.**  1 x 2 MeshBlocks with periodic x1, and 2 x 2 MeshBlocks with
  REFLECTING x1 and `implicit_partition = gather` (periodic x1 across several blocks is
  the clean fatal of sect. 8): **bitwise**, difference exactly 0.000e+00 at
  `implicit_lin_tol` 1e-10 and 1e-13.  Total E is conserved to 5.9e-9 relative, which is
  the precision of the single-precision `bin` dump the sum is taken from, not a drift of
  the solve; the two decompositions give the same number to the last digit.
* **G5 restart**, 2-D 64^2, 2 x 2 blocks, restarted from the mid-run `rst` the
  uninterrupted run itself wrote at t = 9.92: the final double-precision slice is
  **BYTE-IDENTICAL**.  `f0x1` + `f0x2` are therefore the complete persistent state.
* **G6 regression**: `tests_gate_merge/postmerge.sh` -- box G1 modes 3 and 0, **10/10
  files IDENTICAL**; the 1-D He column reproduces the reference line (Picard mean
  11.70892, max 23, `V1max` 1.8343925542664529e4, `F1top/F1bot` 0.99982556 = the
  1.0000097 / 0.9998352 of the reference to eight digits).  `implicit_x1` thick pulse,
  `tau_cell` 1e3, `implicit_cfl` 100, `central`: ratio **1.000001**, Picard mean 2.000 --
  the 3a/3a2/3b/3c number.  Nothing on the `explicit` or `implicit_x1` paths changed
  arithmetically: every new branch is behind `trans`, which is false for both.

### Pass counts

The cost parameter is the DIFFUSION CFL of the transverse direction,
`nu_t = c dt/(3 tau_cell dx)`, not `c dt/dx`: at `tau_cell = 1e3` the 2-D pulse takes 3
passes at `implicit_cfl` 1, 8 at 1e2 and 131 at 1e4, and at `tau_cell = 10` it already
takes 131 at `implicit_cfl` 1e2.  That is the expected behaviour of a Jacobi iteration on
the transverse direction and it is the reason `implicit_maxit` defaults to 200 here.  It
is also the case for a Krylov wrapper (the optional `bicgstab`, not implemented).

### What is NOT done

`implicit_solver = bicgstab` (fatal); the asymptotic-preserving transverse fluxes;
Marshak / imposed-flux transverse boundaries; SMR/AMR; GPU; the 2-D He slab (G7 of the
brief) and therefore any statement about `box_convection`'s M1 wiring in 2-D.

## 12. Findings of 3b phase C (2026-09-21; supersede sect. 4 and 8-11 where they differ)

Phase C ships `<rad_m1>/implicit_solver = bicgstab`, the Krylov wrapper sect. 11 named as
the remedy for the 131-pass line-Jacobi cases, and takes the first look at the 2-D He
slab (G7 of the phase-B brief).  Gate tables and the exact commands:
`tests_m1/runs_3b4/RESULTS.txt` and `runs_3b4/run_gates_3b4.py`.

```
<rad_m1>/implicit_solver  = line_jacobi | bicgstab
          implicit_lin_tol   = 1e-10   the max-norm TRUE residual of the 7-point system,
                                       relative to the max norm of its right-hand side
          implicit_lin_maxit = 200     the BiCGStab iteration cap (new)
```

### The system the two solvers share

The Picard pass already assembles the row
`TA E_{i-1} + TB E_i + TC E_{i+1} = TR`, with `TB` carrying the transverse DIAGONAL
`M1_IW_TDIA` and `TR` carrying the lagged off-diagonal term `M1_IW_TRHS`.  Phase C stores,
in the same kernel, the four (six in 3-D) transverse OFF-DIAGONAL coefficients
`M1_IW_CJM..CKP` that `TRHS` was built from, so the full row is

```
A E = tridiag(TA,TB,TC) E + CJM E_{j-1} + CJP E_{j+1} + CKM E_{k-1} + CKP E_{k+1}
b   = TR + sum_nb C_nb E^k_nb
```

and ONE line-Jacobi pass is exactly `E <- M^{-1}(b - sum_nb C_nb E)` with
`M = tridiag(TA,TB,TC)`.  The two solvers therefore have the SAME matrix, the SAME
right-hand side and the SAME fixed point; only the number of passes differs, which is
what makes the gate comparison meaningful.  `b` is rebuilt from `TR` rather than
accumulated separately, so the identity holds by construction and not by agreement of two
derivations.

The off-diagonal Eddington divergence `sum_{e != d} d_e P_de`, the `F0^n` memory term and
`g0` stay FROZEN inside a pass (they are part of the nonlinearity the outer loop handles),
exactly as line Jacobi freezes them.

### The solver

Matrix-free BiCGStab, RIGHT-preconditioned by the exact x1 line solve (the phase-B Thomas
/ cyclic-Thomas / gathered-stack sweep, now extracted as `ImplicitTridiagSolve` and called
from both solvers -- no arithmetic changed, see the regression below).  Per inner
iteration: two preconditioner applications, two operator applications (one halo exchange
of one Krylov vector each, through a second `MeshBoundaryValuesCC` on a 1-variable scratch
array), and **three global reductions** (`rhat.r`, `rhat.v`, `t.s` with `t.t`), five
scalars; the measured average is 4.0-4.6 reductions per iteration because the max-norm
residual check and the occasional true-residual check add one each.

Stopping is on the TRUE residual: the recursive residual is monitored every iteration and,
when it passes `implicit_lin_tol`, ONE extra operator application checks `max|b - A x|`;
if the two disagree the recurrence is restarted from the true residual.  Breakdown
(`|rho|`, `|rhat.v|` or `|omega|` underflowing) restarts once from the current iterate
with a fresh shadow residual; a third breakdown in one pass falls back to a single
line-Jacobi update and is counted.  The end-of-run line reports outer passes, linear
solves, inner iterations mean/max/total, breakdowns, fallbacks and reductions.

MPI: the Krylov halo is a second exchange object sharing the tag space with the phase-B
one.  It is safe because the two are used strictly sequentially and every rank issues the
identical SEQUENCE of exchanges -- the Picard count, the BiCGStab count and every
breakdown decision come from global reductions -- so MPI's non-overtaking guarantee keeps
the streams apart.  This was reasoned, not measured: all of phase C is serial CPU.

### Gates (serial CPU, `build_cpu_m1`), numbers in `tests_m1/runs_3b4/RESULTS.txt`

* **G1/G2 accuracy.**  `d(sigma^2)/dt / 2D` and the isotropy are the phase-B numbers to
  six digits in every measurable 2-D and 3-D case at `tau_cell` 10 and 1e3 and
  `implicit_cfl` 1 / 1e2 / 1e4, and with `dx2 = 4 dx1` and `dx2 = dx1/4`.  The isotropy
  spread is 1e-9 to 1e-7, i.e. the two solvers differ only at the linear tolerance.
* **Iterations.**  The 131-pass line-Jacobi cases become 5-7 OUTER passes with 5-6 inner
  iterations each; the cheap cases (3-8 passes) become 3-5 outer passes with 0.3-1 inner
  iterations, because the transverse system is then already solved by the first
  preconditioner application.  Wall clock per step improves 7.7x on the expensive 2-D case
  and is a wash (+-10 %) on the cheap ones.
* **G4 decomposition** (2-D, 1 vs 2x2 MeshBlocks, `implicit_partition = gather`) and
  **G5 restart** on a single block: see RESULTS.txt.
* **G6 regression.**  `tests_gate_merge/postmerge.sh`: box G1 modes 3 and 0, **10/10 files
  IDENTICAL**; the 1-D He column reproduces the reference line exactly (Picard mean
  11.70892, max 23, `V1max` 1.8343925542664529e4, `F1top/Fin` 1.0000097,
  `F1bot/Fin` 0.9998352).  Extracting the line solve into its own method changed nothing.

### The 2-D He slab (G7 of the phase-B brief): what it showed

`tests_m1/runs_3b4/he_slab_m1_2d.athinput` is the 1-D column of `runs_3a` widened to
nx2 = 32 cells of the production width `dx2 = 4 dz`, x2 periodic, `transport = implicit`,
true `c`.  **`box_convection` needed no fix**: its M1 wiring (the `a_rad_ref` fill, the
IC reader, the x1 ghost fill and the six user-history columns) is already written over the
full (k,j,i) range and the global-plane history of phase A1 is decomposition-safe.

* **STATIC (no seed).**  `line_jacobi` does NOT converge here: the transverse diffusion
  CFL in the optically thin top of the slab is enormous, and the final 7-point residual
  sits at ~3e1 after the 200-pass cap on EVERY step.  `bicgstab` converges on every step
  (final residual ~1e-10) and is ~1.45x faster in wall clock at the same time.  The slab
  is as quiet as the 1-D column (numbers in RESULTS.txt).
* **SEEDED.**  With any horizontal structure at all the OUTER loop diverges: the lagged
  off-diagonal Eddington terms of the optically thin top are amplified by `(c dt/dx)^2`
  (here 7e3), and the Picard iteration on them has no fixed point.  BiCGStab solves each
  frozen linear system faithfully -- that is not the failure -- but the outer residual
  grows to ~2e3 and the pass count pins at `implicit_maxit`.  **This, not the transverse
  diffusion, is what blocks a production 3-D run**, and phase C does not fix it.

### What is NOT done

The lagged off-diagonal Eddington coupling in optically thin cells (the seeded-slab
blocker above); the asymptotic-preserving transverse fluxes; Marshak / imposed-flux
transverse boundaries; SMR/AMR; GPU; MPI (reasoned, not measured); a coarse space for the
preconditioner, which sect. 4 expects to matter once the MeshBlock count grows.

## 13. Findings of 3b phase D (2026-09-21; supersede sect. 4 and 8-12 where they differ)

Phase D set out to fix the blocker phase C left: with the production entropy seed the
2-D He slab's OUTER (Picard) loop diverges within 2 s of simulated time.  Phase C's
diagnosis was the LAGGED OFF-DIAGONAL Eddington terms `sum_{e!=d} d_e P_de`, amplified by
`(c dt/dx)^2 ~ 7e3` in the optically thin top.  **That diagnosis is wrong, and the
measurement that shows it is one line: dropping those terms entirely
(`implicit_offdiag = none`) does NOT make the seeded slab converge -- it still pins at
`implicit_maxit` with a 7-point residual of 1.6e3 -- while freezing the CLOSURE over the
step (`implicit_closure_lag = step`), with the off-diagonal terms left exactly where
phase C had them, converges in 4 passes with a residual of 2.4e-15.**  What has no
contraction is the Picard iteration on `chi` and `n`, i.e. on the closure itself: each
pass rebuilds the Eddington tensor from the flux of the iterate, and in a thin cell an
O(1) change of `n` changes the face flux by `c E`, which changes `E` again.

```
<rad_m1>/implicit_offdiag  = auto | lagged | operator | none
              auto     = operator where implicit_solver = bicgstab, lagged otherwise
                         (so every line_jacobi and implicit_x1 run is unchanged)
          implicit_closure_relax  = w in (0,1]   chi,n <- (1-w) old + w new per pass
          implicit_closure_relax_thin = false    relax only where theta_f > 1/2
          implicit_closure_lag    = pass | step  freeze chi,n at the START of the step
```

### (1) `operator`: the off-diagonal terms inside the Krylov operator

With the closure frozen inside a Picard pass, `d_e(D_de E)` for `e != d` is LINEAR in
`E'`, so it belongs in the matrix.  `ImplicitOffDiagOp(xc, yc, sgn)` accumulates exactly
the cell-row contribution the face fluxes carry, with the lagged `E` replaced by any
component of the work array; it is called twice per pass: once with `E^k` to ADD
`L_off(E^k)` to the right-hand side (undoing the assembly's lagged term) and once inside
every `ImplicitApplyOp` to put `L_off(x)` on the left.  Both use the same routine, so the
two cannot drift.  The stencil is 9-point in 2-D and 19-point in 3-D and needs NO extra
communication: the Krylov halo is the module's ordinary cell-centred exchange, which
fills edge and corner neighbours anyway.  Boundaries follow the assembly face by face
(imposed-flux x1 faces and reflecting x2/x3 faces carry no off-diagonal term; an EFIX
Dirichlet row is replaced whole).  It is a fatal with `implicit_solver = line_jacobi`.

The system is then NOT an M-matrix, and `E' > 0` is no longer guaranteed.  The smallest
`E` of every linear solve is reduced globally; a non-positive cell drops the REST of that
step to `none` and is counted (`<rad_m1> offdiag=... positivity fallbacks=N min E=...`).

### (2) closure under-relaxation, and (3) the start-of-step closure

`implicit_closure_relax = 0.3` on the seeded slab: still 200/200 passes, residual 6.6e2.
Under-relaxing the closure between passes does NOT recover a contraction (it only slows
the walk), and it is not a remedy.  `implicit_closure_lag = step` is: the Eddington
tensor is then EXPLICIT in time, as in a VET code that reuses the previous step's tensor.
The cost is formally first order in `dt` in the closure alone; the radiation field is
quasi-static on the hydro `dt` and the closure moves by `O(v dt/L)` per step, which is
what the gates below measure directly (7e-6 relative on the static He slab, 2e-6 on the
thick-pulse diffusion rate, 6.7e-2 of the peak on a free-streaming 2-D pulse at
`c dt/dx = 1`, where the field is NOT quasi-static).

### Gates (serial CPU; `build_cpu_m1` for the pulse, `build_cpu_box` for the slab)

Numbers, commands and raw logs: `tests_m1/runs_3b5/RESULTS.txt`.

* **G-oblique** (new; the test sect. 12 did not have).  2-D 64^2 thick-pulse pgen at
  `tau_cell = 0.1`, where the pulse free-streams and its flux direction sweeps every
  angle.  The exact solution stays circularly symmetric, so the `cos 4 phi` amplitude of
  `E` is a direct measure of the off-diagonal error: **explicit 1.7e-2, lagged 1.3e-2,
  operator 1.3e-2, none 2.2e-1** -- `none` is wrong by 22 % of the peak (and by 44 %
  against the explicit scheme in max norm).  `operator` and `lagged` agree **BITWISE**
  (max difference 0.000e+00), which is the proof that the operator form and the lagged
  form have the same fixed point and that the sign bookkeeping above is right.
* **G-static**: the static 2-D He slab at 100 s, against the phase-C reference
  (`F1top/Fin` 0.99998231, `F1bot/Fin` 0.99985481, `max|v1|` 1.056615 v_MLT).
  `operator` with the closure per pass reproduces all three to **eight digits** and takes
  the phase-C pass count (79.21 mean / 139 max), with 0 positivity fallbacks.  With
  `closure_lag = step` the three become 0.99998297 / 0.99985481 / 1.0566227, i.e. 7e-6
  relative, and the cost drops from 79.2 to **4.05 passes** and from 1.487 to 0.108
  s/step: **12x faster** than phase C at the same answer to six digits.
* **G-pulse**: `d(sigma^2)/dt / 2D` at `tau_cell` 10 / `implicit_cfl` 1 and `tau_cell`
  1e3 / `implicit_cfl` 1e4 is 0.996514 and 0.995500 under both `lagged` and `operator` --
  the phase-B/C numbers to six digits -- with isotropy 2e-9 to 7e-10.  With
  `closure_lag = step` they are 0.996512 and 0.995500 at 3 passes instead of 5.
* **G-seed**, the seeded 2-D He slab, `operator` + `closure_lag = step`: it now **runs
  to t = 200 s** (1393 steps, 0.54 s/step) where phase C reached 1.91 s at 11-13 s/step,
  with no dt collapse, 12.75 outer passes mean and a 7-point residual of 5.7e-11.  It is
  a SMOKE TEST and nothing more: 9 steps of 1393 still hit `implicit_maxit` (the same
  nine under `none`, so they are not the off-diagonal terms), and the DYNAMICS ARE NOT
  PHYSICAL -- `KE_2` grows at 0.157 /s over the first 20 s, 75x the convective rate
  `v_MLT/H_p` = 2.1e-3 /s, and by 200 s (0.42 turnover) the slab is supersonic
  (`max|v|/c_s` 2.6, `max|v2|` 740 v_MLT) and has lost 61 % of its emergent flux.  The
  positivity fallback fires on 1388 of 1393 steps, so after its first pass nearly every
  step runs in `none`.  A third arm, `lagged` + `closure_lag = step`, BLOWS UP at cycle
  14 (dt 8e-33, v 4e38): with the closure frozen the three modes order
  `operator > none >> lagged`, which is the phase-C mechanism -- a term amplified by
  `(c dt/dx)^2` in the thin top -- appearing as an instability instead of as a
  non-converging iteration.  That is why the positivity fallback goes to `none`.
  **The next suspect is no longer the transport solve**: the linear system, the lagged
  off-diagonal terms and the closure iteration are all converged to 1e-11.  The 2-D slab
  has two things the 1-D column has not: a TRANSVERSE radiative force `dm2` with no
  well-balanced reference (`force_reference = wb_arad` subtracts `rho a_rad_ref` in x1
  only) and no sponge; the 0.16 /s growth is in `KE_2`, i.e. exactly there.
* **Regression** (`tests_gate_merge/postmerge.sh`): box G1 modes 3 and 0 **10/10 files
  IDENTICAL**; the 1-D He column reproduces the reference line exactly (Picard mean
  11.70892, max 23, NON-CONVERGED 0, `V1max` 1.8343925542664529e4,
  `F1top/Fin` 1.0000097, `F1bot/Fin` 0.9998352).

### Recommended defaults

`implicit_offdiag = auto` (shipped): `operator` under `bicgstab`, `lagged` otherwise.
`operator` is never worse than `lagged` (same fixed point, bitwise on G-oblique) and is
the only form that can be solved implicitly; `none` is a fallback, not a choice -- 22 %
of the peak on an oblique flux.

`implicit_closure_lag` is shipped at `pass`, because `step` moves the static He slab in
the seventh digit and G-static is an eight-digit gate.  **Every 2-D/3-D production run
should set `implicit_closure_lag = step`**: it is what makes a seeded slab converge at
all, and it is 12x cheaper.  `implicit_closure_relax` stays at 1 (it is not a remedy).

### What is NOT done

The 9 non-converged steps of the seeded slab; the transverse radiative force without a
well-balanced reference, which is now the leading suspect for the slab's 0.16 /s `KE_2`
growth; a positivity-preserving form of the 9-/19-point operator (the fallback to `none`
fires on nearly every seeded step, and `none` is 22 % wrong on G-oblique); the
asymptotic-preserving transverse fluxes; Marshak / imposed-flux transverse boundaries;
SMR/AMR; GPU; MPI (still reasoned, not measured); a coarse space for the preconditioner.

## 14. Findings of 3b phase E (2026-09-21; the TRANSVERSE gas-radiation coupling)

Phase D left the seeded 2-D He slab running but unphysical: `KE_2` grew at 0.157 /s,
75x the convective rate `v_MLT/H_p` = 2.1e-3 /s, and the leading suspect was named as
"the transverse radiative force without a well-balanced reference".  Every 2-D/3-D gate
up to phase D had either `gas_feedback = false` or no transverse structure, so the
x2/x3 momentum exchange, the work term, the advective enthalpy flux `A_2 = v_2 E +
(v.P)_2`, the velocity-dependent `E0` corrections and the `wb_arad` reference had never
been measured at all.  Phase E builds the test that measures them, and it **exonerates
the transverse coupling**: it is bitwise as accurate along x2 as along x1.

### T10, the radiation-modified acoustic wave (the new built-in test)

`<problem>/m1_test = radwave` (`src/pgen/tests/rad_m1_tests2.cpp`,
`inputs/tests/rad_m1_radwave.athinput`, analysis `tests_m1/t10_radwave.py`).  One
wavelength of a linear sound wave in a uniform, optically thick,
radiation-pressure-significant medium, on a periodic box, with the direction selectable
(`<problem>/radwave_dir = x1 | x2 | x3 | xy`).  The SAME physical problem is then
propagated along the column direction of the implicit solve and along a transverse one,
so the transverse path is measured against an x1 path that milestones 3a-3b already gate.

With `tau` per wavelength >> 1 and `c >> c_s` the gas and the radiation are in
EQUILIBRIUM DIFFUSION and the mixture has Chandrasekhar's generalised exponents; with
`beta = P_gas/P_tot`,

```
Gamma_1   = beta + (4 - 3 beta)^2 (gamma-1)/[beta + 12 (gamma-1)(1 - beta)]
Gamma_3-1 = (Gamma_1 - beta)/(4 - 3 beta)
c_s^2     = Gamma_1 (P_gas + P_rad)/rho
```

and the right-travelling eigenmode of amplitude `A` is `drho/rho = A cos(phi)`,
`dv = n c_s A cos(phi)`, `dT/T = (Gamma_3-1) A cos(phi)`,
`dE = 4 E (Gamma_3-1) A cos(phi)`, `F = n [(4/3) dv E + (c/(3 kappa rho)) |k| dE
sin(phi)]`, the last term being the diffusive flux that damps the wave.  Code units are
chosen round: `rho0 = T0 = 1` (so `P_gas = 1`), `arad = 3` (so `E = 3`, `P_rad = 1`,
`beta = 1/2`), `gamma = 5/3`, whence `Gamma_1 = 1.4259259259`, `Gamma_3-1 = 0.3703703704`
and `c_s = 1.6887426837300739`; `c_light = 1000 c_s` and `kappa = 1e5` on a unit box give
`c/c_s = 1e3` and `tau_lambda = 1e5`, i.e. the He-box stiffness.  64 cells per
wavelength, amplitude 1e-4, 2 periods, hydro CFL 0.3, so `c dt/dx ~ 300`.

The measurement is the complex Fourier coefficient of the density perturbation at the
box fundamental, `Z(t) = <(rho - <rho>) exp(-i k.x)> -> (A/2) exp(-i omega t - Gamma t)`:
the unwrapped phase gives the phase speed and `ln|Z|` the decay rate.

### Results (serial CPU, `build_cpu_m1`; `tests_m1/runs_3b6/RESULTS.txt`)

```
arm            dir  transport            c_phase      c/c_s-1   decay/period
a_x1_implx1    x1   implicit_x1 (1-D)  1.6884023    -2.016e-04    -6.6561e-02
b_x1_impl2d    x1   implicit   (2-D)   1.68840254   -2.014e-04    -6.6561e-02
c_x2_impl2d    x2   implicit   (2-D)   1.68840251   -2.014e-04    -6.6560e-02
e_xy_impl2d    xy   implicit   (2-D)   1.68813290   -3.610e-04    -9.2813e-02
d1_x1_expl2d   x1   explicit   (2-D)   1.67027218   -1.094e-02    -5.4967e-01
d2_x2_expl2d   x2   explicit   (2-D)   1.67027218   -1.094e-02    -5.4967e-01
f1_x1_thin     x1   implicit, tau=1e2  0.999975039  -4.079e-01    -1.2388e-01
f2_x2_thin     x2   implicit, tau=1e2  0.999975039  -4.079e-01    -1.2388e-01

cross-arm:  b vs a  dc/c = 1.37e-07,  dGamma/Gamma = 4.5e-06
            c vs a  dc/c = 1.20e-07,  dGamma/Gamma = 1.9e-05
            d2 vs d1  dc/c = 0.0,  dGamma/Gamma = 0.0      (BITWISE)
            f2 vs f1  dc/c = 0.0,  dGamma/Gamma = 0.0      (BITWISE)
```

* **THE TRANSVERSE COUPLING IS NOT THE BUG.** The wave rotated onto x2 travels at the
  same speed and damps at the same rate as along x1, to 1.2e-7 and 1.9e-5 -- truncation
  error, not a scheme difference -- and for the two arms that share a mesh shape
  (explicit, and the marginally thick implicit pair) the agreement is **bitwise**.  The
  x2/x3 momentum deposit, the `bmom_half` rule at periodic faces, the derived
  `F_2 = face mean + A_2 E`, the work term `vbar.delta(rho v)` with transverse
  components, the upwinding of `A_2` and the velocity-dependent `E0` correction are all
  exercised here (the radiative force carries a large part of the wave's restoring
  force: without it the speed would be the gas value 1.29099).
* The implicit scheme reproduces the analytic mixture sound speed to **2.0e-4** and the
  45-degree wave to 3.6e-4 (at 64/sqrt(2) = 45 cells per wavelength, and with 1.4x the
  damping, as its resolution implies).  No arm grows.
* The EXPLICIT scheme, gated for the first time in multi-D with gas feedback, is
  isotropic to the last bit but **1.1 % slow and damps 8x too fast** at `tau_cell = 1.6e3`
  -- the ordinary non-asymptotic-preserving `O(c dx)` numerical diffusion.  It is a
  cross-check here, not a gate.
* At `tau_lambda = 1e2` the radiation is a fast heat bath and the wave becomes
  ISOTHERMAL: the measured 0.999975 against `sqrt(P_gas/rho)` = 1 exactly, in both
  directions.  That is a second, independent physical check of the coupling.

### What T10 does NOT cover

Stratification and gravity; the well-balanced effective potential and
`force_reference = wb_arad`, whose reference acceleration is a function of `z` alone and
is subtracted from the x1 momentum ONLY; an optically THIN region, where the
face-eliminated central transverse flux has no HLL dissipation and `c^2 dt/dx` ~ `c` x 7e3
turns a 1e-3 relative `E` asymmetry into a super-luminal face flux; the free surface and
the imposed-flux bottom; the seed itself.  Every one of those is present in the He slab
and absent here.

### The seeded He slab, split by coupling piece

`<rad_m1>/dbg_gas_force` and `dbg_gas_heat` used to be hard-wired to `true` inside
`ImplicitSolve`; phase E honours them there (both default `true`, so every earlier
configuration is bitwise unchanged) and adds `dbg_gas_force_trans`, which keeps the x1
radiative force -- and with it the well-balanced reference -- but hands the gas no
x2/x3 momentum and no transverse work.  That is the one piece of the multi-D coupling
with no hydrostatic reference to be measured against.

Four 200 s arms of the seeded slab, `implicit_offdiag = operator` +
`implicit_closure_lag = step` (growth rates are `d ln KE/dt` over 1 < t < 20 s; the
convective rate is `v_MLT/H_p` = 2.1e-3 /s):

```
quantity                 s_base      s_notrans    s_noheat     s_noforce
                        (production) (no x2 force)(force only) (heat only)
d ln KE_1/dt [1/s]        0.1807       0.1904       0.1621       0.1599
d ln KE_2/dt [1/s]        0.1568       0.3236       0.0981       0.2042
max|v1| / v_MLT           178.76       191.14       108.43       100.12
max|v2| / v_MLT           739.92       107.71       452.47        45.19
max|v| / c_s               2.584        6.207        0.721         0.809
F1top/F_in                0.39389      0.09002      0.29015       0.12985
min E, last dump           199.0        0.0317       1.0e-4        15.65
```

* `s_base` reproduces `runs_3b5/seed_op` **bitwise** (history, Picard 12.75090 / 9
  non-converged, inner 38.62696, breakdowns 151, residual 5.680650e-11), which is the
  proof that honouring the debug switches here is inert at the defaults.
* **`KE_1` grows at 0.16-0.19 /s in every arm**, 80-90x the convective rate, including
  the arm with no radiative force at all.  The growth rate is not carried by any single
  piece of the gas coupling.
* The transverse force is an **amplifier, not the driver**: removing it cuts `max|v2|`
  6.9x and `KE_2` 14x, but raises the x2 growth rate (0.157 -> 0.324 /s) and makes
  everything else worse (`max|v|/c_s` 2.6 -> 6.2, `F1top/Fin` 0.39 -> 0.09, `min E`
  199 -> 0.032) -- as it must, since the momentum the radiation loses is then given to
  nothing.  `s_noforce` is likewise not a clean control under `wb_arad` (the
  well-balanced gravity source subtracts `rho a_rad_ref` expecting the radiation to put
  it back).

### What is NOT done

The seeded slab is still unphysical and the cause is **not** the transverse gas
coupling.  In the order they should be attacked: (1) the optically thin top, where the
face-eliminated CENTRAL transverse flux has no HLL dissipation and `c^2 dt/dx ~ c x 7e3`
turns a 1e-3 relative `E` asymmetry into a super-luminal face flux (the positivity
fallback fires on 1388 of 1393 steps) -- an asymptotic-preserving transverse face flux
is the next build, and T10 at `tau_lambda = 1e2` is where it can be gated; (2) a
horizontal well-balance gate for `force_reference = wb_arad`, whose reference is a
function of `z` alone; (3) the imposed-flux bottom and the Marshak free surface in 2-D
(T10 is periodic in every direction); (4) whether the seed excites the same modes here
as in the two-stream production runs.  Also still open from phase D: the 9 non-converged
steps, a positivity-preserving 9-/19-point operator, transverse Marshak / imposed-flux
boundaries, SMR/AMR, GPU, MPI (reasoned, not measured), a coarse preconditioner space.

---

## 15. The TRANSVERSE realizability limiter (`implicit_trans_limit`, default off)

### What it is

Two new `<rad_m1>` parameters, both inert at their defaults:

| parameter | default | meaning |
| --- | --- | --- |
| `implicit_trans_limit` | `none` | `none` \| `lp`.  `none` takes the OLD expression for the transverse face `theta` term for term, so every earlier configuration is bitwise unchanged (gate G-a). |
| `implicit_trans_fmax` | `1.0` | the reduced-flux cap `fmax` of the `lp` limiter. |

### Why

The face-eliminated transverse (x2/x3) flux of sect. 11,

```
F_f' = theta_f [ F_f^n - c^2 dt G_f - c dt v_f g0_f ],   G_f = gr_f + off_f,
theta_f = 1/(1 + c dt kt_f)
```

has **no free-streaming bound**.  Its steady state is `F_f = -c G_f / kt_f`, the
unlimited diffusive flux `-c grad P/(rho kappa)`, which at `c dt/dx ~ 7e3` is
super-luminal wherever `rho kappa` is tiny.  In the optically thin top of the seeded 2-D
He slab (sect. 14) that is exactly what happens: `|F_2|/(c E)` saturates at the
post-solve clip of 1 in the top ~12 cells, `E` varies horizontally by 2-5 decades, the
x1 transport collapses (`F1top/Fin` 0.07-0.66) and the positivity fallback fires on
nearly every step.

### The formula

Under `lp` every transverse face carries a LAGGED "limiter opacity"

```
klim_f = |G_f| / (phi_f E_f),
E_f    = (E_L + E_R)/2 of the lagged iterate (floored),
phi_f  = fmax sqrt(max(1 - f1_f^2, 0.01)),
f1_f   = face mean of the lagged cell x1 reduced flux F1/(c E), clipped to [-1,1],
theta_f = 1/(1 + c dt (kt_f + klim_f)).
```

Then in steady state `|F_f| = c|G_f| / (kt_f + |G_f|/(phi_f E_f)) <= c phi_f E_f`: the
transverse flux can never exceed the room the x1 flux leaves in the realizability cone.
Where `R = |G_f|/(kt_f E_f) << 1` -- every optically thick face -- the change is `O(R)`
and the diffusion limit is untouched.  This is a Levermore-Pomraning flux limiter
written as an opacity, so that it enters the already-linear face elimination.

`klim` only ever makes `theta_f` SMALLER, and `theta_f` multiplies the transverse
off-diagonals of the 7-point row together with their own diagonal contribution, so the
row keeps its signs and its diagonal dominance: the M-matrix property of sect. 11 is
preserved.

`theta_f` is computed ONCE per evaluation, in `RadiationM1::ImplicitTransTheta`, into
the face arrays `thx2`/`thx3`, and the face-flux kernels, the cell-terms kernel
(`TDIA`, `CJM..CKP`), `ImplicitOffDiagOp` (hence the BiCGStab operator) and the
post-solve reconstruction all READ that one number, so they cannot drift apart.  `klim`
follows `implicit_closure_lag`: frozen at the entry state for the whole step under
`step`, recomputed each Picard pass under `pass`.  `M1_IW_F1` was added to the
transverse halo (`M1_NHALO_T` 13 -> 14) so that both blocks of a shared face build
`f1_f` from bit-identical numbers; nothing else reads it.

### Gates (serial CPU, `tests_m1/runs_3b8/RESULTS.txt`)

* **G-a, default-off inertness**: the seeded slab to `tlim = 3` gives `m1slab.hydro.hst`
  and `m1slab.user.hst` BYTE-IDENTICAL to the pre-change binary.
* **G-b(i), the static slab** (`vpert = 0`, `tlim = 20`): limiter on = limiter off,
  bitwise, in both history files -- with no horizontal structure `G_f = 0`, so
  `klim = 0`.
* **G-b(ii), T10** along x2 and at 45 degrees with the limiter on: still PASS, with the
  phase speed unchanged in the digits T10 reports.
* **G-c, the cure** (seeded slab, `tlim = 30`): at `fmax = 0.5`, `KE_2(30 s)` falls from
  `1.7e29` to `6.8e26` (2.4 decades), `F1top/Fin` from 0.07-0.66 to 0.95-1.00, the
  Picard count from 12.5 mean / 200 max with 1 non-converged step to 7.6 / 10 with
  none, and `dt` stops eroding (constant 0.1613 s against 0.033 s at 30 s).  What the
  limiter does NOT fix: the positivity fallback still fires on 172 of 186 steps and the
  per-cell `|F_2|/(c E)` still reaches 1 in the topmost cells -- the bound is on the
  face-MEAN `E`, and with `E_min/E_mean ~ 0.03` across a row the thin cell can still sit
  on the clip.  A positivity-preserving operator and/or an AP transverse face flux are
  still needed.

## 16. ANDERSON acceleration of the Picard map (`implicit_accel`, default off)

### What it is for

`transport = implicit` solves the face-eliminated scalar-`E` system with the closure
`(chi, n)` LAGGED, either per Picard pass (`implicit_closure_lag = pass`) or frozen over
the step (`= step`).  On the seeded 2-D He slab at `c dt/dx ~ 7e3` neither works: per pass
the Picard iteration DIVERGES in the optically thin top (under-relaxation to 0.3 diverges
too), and per step each step converges but the run is unstable step to step (sect. 15).
The mechanism is that terms like `c^2 dt d_1(P_12)`, `P_12 = (3 chi - 1)/2 n_1 n_2 E`, act
as an advection of the transverse flux treated EXPLICITLY at that CFL, so the fixed-point
map has gain >> 1.  A CONVERGED nonlinear (backward-Euler) closure should be stable, and
the cheapest way to converge a divergent fixed-point map without building a Jacobian is
Anderson acceleration (Walker & Ni 2011, SIAM J. Numer. Anal. **49**, 1715), which is
equivalent to a multisecant quasi-Newton method on `g(x) = G(x) - x`.

### The fixed-point vector

One Picard pass is the map `x -> G(x)` with

    x = ( E , F_1 , F_2 , F_3 )   per cell,

the four quantities the TOP of a pass rebuilds the closure from (`rf = |F|/(cE)` gives
`chi`, `n = F/|F|`, hence the whole Eddington tensor, `D_ab`, the enthalpy-flux
coefficients `a_d` and `de0`).  Nothing else a pass reads is independent state: the face
fluxes `f0x1/f0x2/f0x3` are recomputed from `E` at the end of every pass, the velocities
and (by default) the opacities are frozen over the step, and `T'` is slaved to `E'` by the
scalar root find.  Accelerating `E` alone would be inconsistent -- the closure would still
be built from the unaccelerated `F`.

**Scaling.**  `E` and `F` are not commensurate (`F ~ c E`) and `E` spans ~10 decades
between the base and the thin top, which is exactly where the iteration needs the help.
Every cell is therefore divided by its OWN energy scale `S = max(E^n, e_floor)`, frozen
over the step so that the least-squares problem stays linear across passes, and the flux
components by `c S`:

    x = ( E/S , F_1/(c S) , F_2/(c S) , F_3/(c S) ) ,

i.e. the Anderson residual is the RELATIVE fixed-point residual cell by cell.  On an
x1-only solve `iw` has no `F_2/F_3` slot and the vector is `(E/S, F_1/(cS))`.

### The algorithm

With `g_k = G(x_k) - x_k` and the histories `dX_t = x_{t+1} - x_t`, `dG_t = g_{t+1} - g_t`
of the last `m` passes,

    gamma   = argmin || g_k - dG gamma ||_2 ,
    x_{k+1} = x_k + beta g_k - (dX + beta dG) gamma .

`beta = 1` is the undamped form; `gamma = 0` (empty history) is plain Picard.  The least
squares is solved on the HOST from the normal equations `(dG^T dG + lambda I) gamma =
dG^T g_k` with `lambda = 1e-10 tr(dG^T dG)/m` (Tikhonov, so a nearly dependent history is
harmless), `m <= 10`.  Every inner product is a GLOBAL sum: one `MPI_Allreduce` of `m+1`
doubles per history row (`M1GlobalSumArr`), so the answer does not depend on the
decomposition.  The histories are device Views in a ring buffer, allocated only when the
acceleration is on.

Safeguards:

* the accelerated iterate has the realizability limits re-applied (`E >= e_floor`,
  `|F| <= c E`) before it becomes the next lagged state;
* if `||g_k||` grows by more than 10x over the previous pass the history is DROPPED and
  that pass falls back to plain Picard (counted as a "history restart");
* the last pass of a step is never accelerated, so the state a step ends on is the one the
  stored face fluxes were built from;
* the convergence test is UNCHANGED (`implicit_tol` on `|dE|/E` plus the true 7-point
  linear residual against `implicit_lin_tol`), measured on the UNACCELERATED map, which is
  the fixed-point residual itself.

### Parameters (all inert at their defaults)

    implicit_accel          = none | anderson   (default none = bitwise the 3b/3c code)
    implicit_anderson_m     = 5                 (history depth, 1..10)
    implicit_anderson_beta  = 1.0               (mixing, (0,1])
    implicit_anderson_start = 1                 (first 0-based pass that is accelerated)

`ImplicitReport` prints `accelerated passes` (total and per solve) and `history restarts`
next to the existing Picard mean/max/NON-CONVERGED line.

### Gates (serial CPU, `tests_m1/runs_3e_newton/RESULTS.txt`)

* **N-a, default-off inertness**: the seeded slab to `tlim = 3` with
  `implicit_closure_lag = step` gives `m1slab.hydro.hst` and `m1slab.user.hst`
  BYTE-IDENTICAL to the stored pre-change output, and every solver statistic agrees
  (19 solves, Picard mean 8.157895, max 11, 14 positivity fallbacks, min E -5.957282e6).
* **N-b, the static slab** (`vpert = 0`, `tlim = 20`, `lag = pass`): `m = 5`, `beta = 1`
  cuts the Picard count from **80.35 mean / 129 max to 22.74 / 50** (3.5x / 2.6x) with
  0 non-converged steps in both and 16 history restarts, and the BiCGStab reductions
  fall 155120 -> 75688.  The answer is unchanged: mass agrees to 1.1e-14, total energy
  to 6.3e-12, the transverse momenta and kinetic energies exactly; the SIGNED,
  zero-crossing `1-mom` and `V1mid` agree to 2.0e-8 of their own amplitude at t = 1.1,
  drifting to 1.6e-7 by t = 20 -- the accumulation of two iterations stopping at
  different points inside the same `implicit_tol = 1e-8`, not a defect.
* **N-c, the seeded slab: FAIL.**  With `lag = pass` and the closure per pass, *every*
  step of *every* arm exhausts `implicit_maxit = 200` and is NON-CONVERGED -- no
  acceleration (72/72 steps), `m = 5, beta = 1` (56/56), `m = 10, beta = 0.5` (42/42),
  and `m = 5` plus `implicit_trans_limit = lp, fmax = 1` (19/19).  The residual stalls
  at a MEDIAN `|dE|/E ~ 1.9e2` with excursions to 1e8-1e10, always in the optically
  thin top: the worst cell's depth index is confined to `i = 75..86` of the active
  range `3..86` (the top 12 cells, tau < ~0.1), never below.  The positivity fallback
  fires on 100 % of the steps in all four arms and the linear solve returns
  `min E = -3e6 .. -3e7`, i.e. the closure is handed a non-realizable state every pass.
  What the acceleration does buy is a monotone improvement in the final 7-point linear
  residual (mean 1.0e3 with none, 1.5e2 at `m = 10, beta = 0.5`, 7.2e1 with the lp
  limiter) and, with the limiter, a `dt` that stays at its initial 0.1613 s.
  **Conclusion:** the stall is not the RATE of the outer iteration, so accelerating this
  map further is not the way forward; the non-M-matrix 9-point off-diagonal operator and
  the unbounded transverse face flux are.  A positivity/asymptotic-preserving transverse
  flux, or a true Newton step on the coupled `(E,F)` system with the closure Jacobian,
  is what the evidence asks for.

`ImplicitSolve` also gained a diagnostic (unconditional, but printed only for a step
that fails to converge): one line with the pass count, the Picard and linear residuals
and the `(m,k,j,i)` of the worst cell.  It changes no number.
