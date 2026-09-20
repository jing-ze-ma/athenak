# Grey M1 radiation transport (`<rad_m1>`): design note

Status: DESIGN ONLY, no code yet (2026-09-21, branch `rt-integration`).

Today every radiative setup in this tree (box_convection, red_giant, deep_hot_jupiter) uses

```
radial:      grey / correlated-k two-stream column (src/utils/two_stream_rt.hpp),
             mode-3 implicit column, blended with radiative conduction by tau
transverse:  radiative conduction, STS or ADI (src/diffusion/conduction_transverse.cpp)
force:       (1-w) rho kappa F/c + P_rad grad w, tied to the EOS radiation taper
```

That is three operators stitched together, with the radiation energy and pressure living
inside the EOS.  This note specifies a replacement in which the radiation field is
evolved: a non-relativistic, grey, photon two-moment (M1) module.

**Scope of this note**: stage 1 and stage 2 of the plan in section 9, i.e. Cartesian,
explicit transport, implicit local matter coupling, uniform mesh and SMR/AMR through the
existing cell-centred machinery.  **Out of scope, stated so that nothing here forecloses
it**: implicit transport (stage 3), spherical-polar and cubed-sphere geometry (stage 4),
multigroup / correlated-k, MHD coupling, scattering beyond an isotropic `kappa_s`.

The upstream GR neutrino module (`origin/project/ccsn:src/radiation_m1`, ~6000 lines)
was read for reuse.  Verdict: take the task list, the cell-centred boundary wiring and the
module skeleton; write the physics fresh (55-60 % of it is metric and neutrino
machinery; it has no hydro coupling outside `<adm>`, no restart support, and
`A_jp12` uses `dx1` in all three directions).  Its root finders look like GSL ports
(GPL) and are not needed in flat space: do not copy them.

---

## 1. Equations

Mixed frame: lab-frame moments `E`, `F_i`, comoving-frame opacities, all in code units.
`beta = v/c`.  The reduced speed of light `chat` multiplies the radiation time derivatives
only; every source term and the reduced flux `f = |F|/(cE)` keep the true `c`.

```
dE/dt   + (chat/c) div F     = - chat   G0
dF_i/dt + chat c   d_j P_ij  = - chat c G_i
d(E_gas)/dt = ... + c G0            d(rho v_i)/dt = ... + G_i
```

so what is conserved is `E_gas + (c/chat) E` and `rho v + F/(chat c)`; at `chat = c`
these are the true totals.

**Source terms.**  They are written by evaluating the four-force in the comoving frame and
transforming it, not as a term-by-term expansion:

```
E0     = (1 + beta^2) E - 2 beta.F/c + beta.P.beta        comoving energy density
F0_i/c = F_i/c - beta_i E - beta_j P_ij                   comoving flux
g0     = rho (kappa_E E0 - kappa_P a T^4)                 comoving energy exchange
g_i    = rho (kappa_F + kappa_s) F0_i/c                   comoving momentum exchange
G0     = g0 + beta.g            G_i = g_i + beta_i g0
```

`kappa_P` Planck (emission), `kappa_E` energy (absorption) mean, `kappa_F` flux
(Rosseland) mean.  Expanded, `G0` carries `(kappa_F - 2 kappa_E) beta.F/c`, which is
Krumholz, Klein & McKee (2007) `(kappa_0R - 2 kappa_0P) v.F/c^2`, and
`-(kappa_F - kappa_E)(beta^2 E + beta.P.beta)`; this is AREPO-IDORT (Ma et al.,
arXiv:2503.16627) eqs. (56)-(57) up to the `(1 + beta^2/2)` Lorentz factor on the
emission term, which is dropped here (relative size `beta^2 ~ 1e-9`, never promoted).

Why this form and not the usual O(v/c) truncation (Skinner & Ostriker 2013,
Rosdahl & Teyssier 2015 main text, AREPO-RT, QUOKKA I):

* In dynamic diffusion (`tau beta > 1`, which is every deep cell of ours) the lab flux is
  `F ~ (4/3) v E + F0`, so `beta.F/c` and `beta^2 E` are the SAME order and cancel to leave
  the physical work term `beta.F0_diff`.  Keeping `beta.F/c` without `beta^2 E` leaves a
  spurious heating `~ (4/3) c rho kappa beta^2 E`.  Writing `g_i` with `F0` makes the
  cancellation exact by construction.  This is the Krumholz et al. argument that terms
  formally O(beta^2) are leading order at `beta tau > 1`.
* A single opacity (`kappa_P = kappa_E = kappa_F`) gives the `beta.F/c` coefficient `-1`
  instead of `kappa_F - 2 kappa_E`; with a stellar table the means differ by orders of
  magnitude near the photosphere.
* What falls out of the form above without further terms: advection of trapped radiation
  enthalpy (the stiff `g_i` relaxes `F` to `(4/3) v E + F0`), the `v.grad P_rad` work term
  (`beta.g` with `F0 = -c grad E/(3 rho kappa_F)`), and the radiation force `G_i`.

Not kept: the `(v.F) v/c^3` and `(1/2) beta^2 kappa F/c` momentum terms of Krumholz et al.
(they are O(beta^2) corrections to an O(1) term and are never promoted).

**Closure** (Levermore 1984), identical in every M1 code surveyed:

```
chi = (3 + 4 f^2) / (5 + 2 sqrt(4 - 3 f^2)),   f = |F|/(c E),   n = F/|F|
P_ij = E [ (1-chi)/2 delta_ij + (3 chi - 1)/2 n_i n_j ]
```

Algebraic in flat space (no root find).  The closure uses the LAB-frame `f`; the comoving
correction to `chi` is O(beta) on a term that is already isotropic where beta tau matters.

**The EOS becomes gas-only.**  `<hydro>/eos_radiation = false` and none of
`eos_rad_rho_hi/lo`, `eos_rad_t_hi/lo`.  The taper is not an alternative (it removes
radiation only in the thin wing), and `problem/rt_rad_force` is fatally tied to the taper
(box_convection.cpp:2193), so that path is refused together with `<rad_m1>`.  Consequence
for initial conditions: a hydrostatic profile built with radiation pressure in the EOS
must be rebuilt as gas pressure + `E/3`, with `E = a T^4` and `F` from the diffusion limit.

---

## 2. Reduced speed of light: where it is and is not allowed

Consensus condition (Skinner & Ostriker 2013 eq. 35; Melon Fuksman et al. 2021):

```
chat >> v_max * max(1, tau_max)
```

and RSLA is wrong at O(1) in dynamic diffusion: with `chat` on the time derivative the
trapped radiation is advected at `(chat/c)(4/3) v`, and the diffusion time is long by
`c/chat`.  For a He-star envelope (`v ~ 1e6`, `tau ~ 1e4-1e6`) and for the deep hot
Jupiter the condition gives `chat ~ c`.  Therefore:

* `<rad_m1>/chat_over_c` exists (default 1) and is a TEST and shallow-box tool.
* At startup the module evaluates `v_max tau_max / chat` over the domain and prints it;
  above `<rad_m1>/rsla_warn` (default 0.1) it warns, above 1 it is fatal unless
  `rsla_force = true`.
* Explicit M1 is a validation vehicle for the closure, flux, sources and coupling.
  Production on our stars needs stage 3 (implicit transport), or a hybrid that hands
  `tau > tau_hi` to the existing implicit conduction; that decision is taken after stage 2.

---

## 3. Transport flux

Reconstruct `(E, f_i)` with the existing PLM (`src/reconstruct/plm.hpp`) and rebuild
`F = c E f`, so `|f| <= 1` survives reconstruction (QUOKKA).  Wave speeds: the closed form
of Skinner & Ostriker (2013, eq. 41a), no table:

```
lam_pm/c = { mu f +- sqrt[ (2/3)(4 - 3f^2 - s) + 2 mu^2 (2 - f^2 - s) ] } / s,
s = sqrt(4 - 3 f^2),   mu = cos(angle between F and the face normal)
```

HLL with `b_L = min(0, lam_-)`, `b_R = max(0, lam_+)` taken over both states.

**The thick limit is the main trap.**  Plain HLL has `D_num = c dx/2` against the physical
`c/(3 rho kappa)`: wrong as soon as `tau_cell = rho kappa_F dx > 1`; ours reach 1e3-1e6.
Two corrections are implemented behind `<rad_m1>/thick_flux`, because they share
everything else and the gates decide:

`ap_hll` (default; Berthon, Charrier & Dubroca 2007):

```
alpha = (b_R - b_L) / (b_R - b_L + c tau_face),   tau_face = (rho kappa_F)_face dx
F_E   = alpha * [ (b_R F_L - b_L F_R) + b_R b_L (E_R - E_L) ] / (b_R - b_L)
```

In the thick limit `alpha -> (b_R-b_L)/(c tau)`, `b = +-c/sqrt3`, and the dissipation term
becomes `-(c/(3 rho kappa dx))(E_R - E_L)`: the compact two-point physical diffusion
flux.  The `P`-flux and the source remainder follow Berthon's eqs. (UNVERIFIED here: read
the paper before coding; the survey verified the `alpha` formula and the limit only).

`blend` (Radice et al. 2022; upstream AthenaK): `F = F_c - A (1-phi)(F_c - F_LLF)`,
`A = min(1, 1/tau_face)`, minmod ratio `phi`, plus their sawtooth switch.  Cheaper, but the
thick-limit diffusion operator is centred on a `2 dx` stencil, i.e. it decouples odd and
even cells; the sawtooth switch is there for that reason.  We have been bitten by
checkerboards before (B-star photosphere), hence not the default.

`none`: plain HLL, for the beam and shadow tests and as the failing control of T3.

**Moving fluid.**  Bloch et al. (2021) state that the AP correction is not AP in a moving
fluid: the advective enthalpy flux `v E + v.P` sits inside `F` and is multiplied by
`alpha -> 0`.  Fix, to be proven by T4: split the face flux,

```
F_E = upwind_v[ v_n E + (v.P)_n ]  +  ap_hll[ F - v E - v.P ]
```

with the hydro face velocity; the same split for the `P`-flux is NOT obviously right and is
the open design question of this note (section 8, R2).

Face opacity: harmonic mean of `rho kappa_F` (the diffusion-limit flux across an opacity
jump, e.g. the Fe bump, is set by the harmonic mean).

---

## 4. Matter coupling: the implicit local solve

Per cell, after the explicit transport update (`*` values), opacities frozen at `T^n`
(one optional Picard pass, `<rad_m1>/opac_picard`).  Three steps, none iterative except a
scalar root find:

**(a) energy exchange**, the stiff part `g0` only.  Backward Euler on `E`:

```
E' = [ E* + chat dt rho kappa_P a T'^4 ] / (1 + chat dt rho kappa_E)
y(T') = rho e(rho,T') + (c/chat) E'(T') - [ rho e* + (c/chat) E* ] = 0
```

`y` is strictly increasing in `T'`; bracketed Newton (bisection fallback, tol 1e-12) using
`EOS_Data::ThermoAt` for `e` and `c_v` (eos.hpp:604).  `T^4` is NOT linearised: that is
what makes the tabulated EOS converge (AREPO-IDORT; Skinner & Ostriker's quartic is the
ideal-gas special case).  Energy conservation is algebraic, not tolerance-limited: the gas
energy is SET to `total - (c/chat) E'` (QUOKKA 2024b).

**(b) flux**, linear backward Euler with `E'`, `P' = P(E', f*)`:

```
F' = [ F* + chat dt rho (kappa_F+kappa_s)(v E' + v.P') - chat c dt beta g0' ]
     / (1 + chat dt rho (kappa_F + kappa_s))
delta(rho v) = - (F' - F*) / (chat c)
```

then re-apply `|F'| <= c E'`.

**(c) work term**, explicit (non-stiff: it is `v.grad P_rad`): `W = vbar.delta(rho v)` with
`vbar` the mean of old and new velocity (exactly the kinetic energy change); `u(IEN) += W`
and `E' -= (chat/c) W`.  The `E0 - E` velocity corrections inside `g0` are applied
explicitly with `*` values (they are O(beta) of a term that vanishes in equilibrium).

Floors: `E >= <rad_m1>/e_floor`; flux limiting scales `F` by `c E/|F|` (upstream scales by
the square of that ratio: a bug, not copied).  Both are counted in `EventCounters`.

**Opacities.**  `kappa_F` = the existing Rosseland table (`RosselandTable`,
conduction.hpp:750).  A Planck mean does not exist in this tree (correlated_k.hpp:33 flags
it).  Stage 1 tests use constant/power-law opacities; stage 2 starts with
`kappa_P = kappa_E = kappa_F` behind `<rad_m1>/planck_from_rosseland = true` and a second
table `<rad_m1>/planck_table` is added before any result is quoted.

---

## 5. Time integration and task lists

Radiation is sub-cycled inside the operator-split lists, modelled on the resistive RKG
loop (driver.cpp:430-445): new lists `m1_before_stagen / m1_stagen / m1_after_stagen`,
`N_sub = ceil(dt_hydro / dt_rad)`, `dt_rad = cfl_rad * min(dx)/chat`, placed after
`after_timeintegrator`; with `<rad_m1>/strang = true` half the substeps run in
`before_timeintegrator` (the Strang finding for the two-stream, that the unsplit form
seeds the closure before the solve, applies here too).  `dtnew` is NOT folded into
`Mesh::NewTimeStep` unless `subcycle = false`.

Each substep is IMEX: target PD-ARS (Chu et al. 2019 as used by QUOKKA, arXiv:2404.08247;
AP, implicit stage cell-local; exact tableau UNVERIFIED, read before coding).  First
implementation: two-stage SSP-RK2 transport with the section 4 solve after each stage.

`m1_stagen` chain, the upstream order: `CopyCons -> Closure -> Opacity -> Fluxes ->
SendFlux -> RecvFlux -> Update -> Coupling -> RestrictU -> SendU -> RecvU -> PhysicalBCs
-> Prolongate`.  The coupling writes hydro `u0`, so one hydro `ConToPrim` closes each
substep block (not each substep: `T` for the next substep comes from the solve itself).

---

## 6. Plumbing (file list)

New: `src/rad_m1/rad_m1.{hpp,cpp}`, `_tasks.cpp`, `_fluxes.cpp`, `_update.cpp`,
`_coupling.cpp`, `_closure.hpp`, `_newdt.cpp`; `src/bvals/physics/rad_m1_bcs.cpp`;
`src/pgen/tests/rad_m1_tests.cpp`; all added by hand to `src/CMakeLists.txt`.

Touched: `meshblock_pack.{hpp,cpp}` (`<rad_m1>` branch after :188, fatal with
`<radiation>` and with `eos_radiation`), `driver.cpp` (sub-cycle loop), `mesh.cpp:890`,
`mesh_refinement.cpp:517,554,591` and `load_balance.cpp` (~7 sites),
`basetype_output.cpp` (`m1_e`, `m1_f1..3`, derived `m1_f`, `m1_chi`, `m1_tr`),
`restart.cpp:73-85,173-179,306-331` with the mirror in `pgen.cpp:155-165,323-342,577+`
(same order in writer and reader; upstream has no restart support at all).

Variables: `u0(m, n, k, j, i)`, `n = 0..3` = `E, F1, F2, F3`, one
`MeshBoundaryValuesCC` with `InitializeBuffers(4)`, `coarse_u0` when multilevel.

Boundary states (ghost fills read the conserved radiation variables, never `w0`: the
restart-bitwise rule of box_convection.cpp:3510): top `vacuum` = free streaming,
`F_n = c E` outward, zero incoming, the analogue of `rt_top_vacuum`; bottom `flux` =
imposed `F_n = <hydro>/rad_flux_inner` with `E` from the diffusion limit,
`E_g = E_1 + 3 rho kappa_F dx F_n/c`; plus `reflect`, `outflow`, `periodic`, `user`.

---

## 7. Bitwise requirement

With no `<rad_m1>` block nothing may move a bit: the module pointer is null, the new task
lists are empty, the driver loop is skipped.  Proof obligation: the box G1 gate
(`tests_r8/g1.sh`, modes 3 and 0) and `red_giant_cs` 100 cycles, against `e688efd4`.

---

## 8. Test plan

All in `rad_m1_tests` (`problem/m1_test = ...`), inputs under `inputs/tests/rad_m1_*`,
regression versions under `tst/test_suite/rad/test_rad_m1_*_cpu.py`.

**(T1) Beam.**  128^2, 45 degrees, `kappa = 0`, `thick_flux = none`.  FWHM 24 cells with
the computed eigenvalues (30 with `lam = +-c`; Gonzalez et al. 2007).  Pass: <= 25.

**(T2) Shadow** (Hayes & Norman; HERACLES numbers).  1 x 0.12 cm, clump at (0.5, 0) of
semi-axes (0.1, 0.06), `rho_0 = 1`, clump 1000x, `sigma = 0.1 (T/T0)^-3.5 (rho/rho0)^2`
per cm, `T0 = 290 K`, source `T_r = 1740 K`, 280 x 80.  Pass: light-crossing 3.33e-11 s
within 5 %, shadow `E < 1e-3` of the lit value behind the clump.

**(T3) Static thick pulse: the AP gate.**  1-D Gaussian in `E`, scattering-dominated,
`tau_cell` = 0.1, 10, 1e3, 1e6, no hydro.  Measure the variance growth rate against
`2 c/(3 rho kappa)`.  Pass: within 2 % at every `tau_cell` for `ap_hll`; second order
between two resolutions at `tau_cell = 10`; control: `none` must FAIL at `tau_cell >= 10`
(if it passes, the test is not testing).  Add a Nyquist-mode seed at `tau_cell = 1e3`:
pass = monotone decay, no growth, for `ap_hll`; record what `blend` does.

**(T4) Advected thick pulse: dynamic diffusion.**  QUOKKA AP-paper parameters: static
`tau = 2.9e3, beta = 3.3e-5`; dynamic `tau = 1.4e4, beta = 1e-3` (`beta tau = 14`),
periodic, prescribed `v`.  Pass: pulse centre at `x0 + v t` within one cell, width equal
to the static case within 2 %, no net gas heating beyond 1e-10 relative (this is the test
that the `beta^2` terms cancel).  Control: the O(v/c)-truncated source must show the
spurious heating.  Run at `chat = c` and once at `chat = c/10` to exhibit the
`(chat/c)(4/3) v` advection error of section 2.

**(T5) Equilibration.**  Single cell, HERACLES numbers (`E_r = 1e12`, `sigma = 4e-8`,
`e_0 = 1e2` and `1e10`), `dt` from 0.1 to 1e5 coupling times.  Pass: correct equilibrium
for every `dt`, no sign oscillation, `E_gas + (c/chat) E` conserved to round-off; repeat
with the tabulated He EOS at (`rho = 1e-8`, `T = 2e5 K`).

**(T6) Marshak wave / Su-Olson.**  Bloch et al. numbers: pass <= 2 % L1 (their AP scheme
1.1 %, uncorrected HLL 84 %).

**(T7) Radiative shocks** (Lowrie & Edwards 2008), Mach 2 subcritical and Mach 5
supercritical, coupled to hydro.  Pass: L1 < 2 % in `rho`, `T_gas`, `T_rad` at `chat = c`.
Note Skinner & Ostriker's result that fixed `lam = +-c/sqrt3` beat computed eigenvalues
here (0.4 % vs 2-8 %): record both.

**(T8) Conservation and restart.**  Any coupled test, 200 cycles: totals to round-off;
restart at cycle 100 bitwise against the uninterrupted run; 1 rank vs 4 ranks bitwise.

**(T9) Radiation-supported atmosphere.**  Plane-parallel hydrostatic column with
`F = const` and gas + radiation pressure, `tau_cell ~ 100`: velocities stay below 1e-6 of
the sound speed for 100 sound crossings.  This is the precursor of the box.

Open risks, to be closed by the gates above and not before:

* **R1** thick-limit flux: `ap_hll` needs Berthon's `P`-flux and source remainder read
  from the paper; T3 decides between `ap_hll` and `blend`.
* **R2** the advective split of section 3 in the `F` equation; T4 decides.
* **R3** PD-ARS tableau and its diffusion-limit order; T3/T6 at radiation CFL > 1.
* **R4** Planck-mean table for stellar mixtures: source to be chosen (OPAL/OP do not
  ship one; the correlated-k route exists for the hot Jupiter only).

---

## 9. Staging

1. Stage 1: module + T1-T6, T8 (radiation alone or with prescribed `v`); section 7 gate.
2. Stage 2: hydro coupling, T7, T9; then a shallow box_convection photosphere with
   `chat` inside the section 2 condition, compared with the two-stream at equal setup:
   emergent flux map, `T(tau)`, radiation force.
3. Stage 3 (separate note): implicit transport, Eddington tensor lagged, matrix-free
   Krylov, preconditioned by the block-tridiagonal radial column solve and the existing
   ADI line partition.
4. Stage 4 (separate note): spherical-polar, then cubed sphere (geometric sources for
   `F`, seam basis transform as for the momentum).

## References

Audit et al. 2002 (astro-ph/0206281); Berthon, Charrier & Dubroca 2007 (J. Sci. Comput.
31, 347); Gonzalez, Audit & Huynh 2007 (A&A 464, 429); Krumholz, Klein & McKee 2007
(astro-ph/0611003); Lowrie, Morel & Hittinger 1999 (ApJ 521, 432); Jiang, Stone & Davis
2012 (1201.2223); Skinner & Ostriker 2013 (1306.0010); Sadowski et al. 2013 (1212.5050);
Rosdahl & Teyssier 2015 (1411.6440); Foucart et al. 2015 (1502.04146); Kannan et al.
2019 (1804.01987); Melon Fuksman & Mignone 2019 (1903.10456), Melon Fuksman et al. 2021
(2005.01785); Jiang 2021 (2102.02212); Bloch et al. 2021 (2011.13926); Wibking &
Krumholz 2022 (2110.01792); Menon et al. 2022 (2202.08778); Radice et al. 2022
(2111.14858); He, Wibking & Krumholz 2024 (2404.08247, 2407.18304); Ma, Pakmor, Justham
& de Mink 2025 (2503.16627).
