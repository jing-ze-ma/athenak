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
