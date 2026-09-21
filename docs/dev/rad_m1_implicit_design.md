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
