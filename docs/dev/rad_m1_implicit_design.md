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
