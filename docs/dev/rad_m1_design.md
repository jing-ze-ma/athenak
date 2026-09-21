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
The literature has three families; read at formula level (Berthon & Turpault 2011 in
full, Bloch et al. 2021 incl. appendices, Jiang 2021 App. A, He, Wibking & Krumholz 2024,
Rosdahl & Teyssier 2015; NOT read: Berthon, Charrier & Dubroca 2007, Jiang, Stone & Davis
2014, Gonzalez et al. 2007; VETTAM and AREPO-IDORT only through summaries).

`<rad_m1>/thick_flux = ap_hll` (default).  The asymptotic-preserving factor multiplies the
E-flux only; the pressure flux is plain HLL (the Bloch et al. 2021 variant, eq. 12/C.1):

```
alpha = 1 / [ 1 - 3 tau_face (1 - f^2) lam+ lam- / (c (lam+ - lam-)) ]     (Bloch eq. 25)
F_E   = alpha * F_E^HLL(reconstructed L, R)  +  (1 - alpha) * F_diff  + A_upwind
F_diff = - c [ E(i+1) - E(i) ] / (3 tau_face)            CELL-CENTRE values
F_F   = c^2 [ (b_R P_L - b_L P_R) ] / (b_R - b_L) + b_R b_L (F_R - F_L) / (b_R - b_L)
tau_face = (1/2) [ (rho kappa_F)_L + (rho kappa_F)_R ] dx
```

**Second order.**  Berthon's scheme is first order, and no second-order extension is in
print.  Its original form, `alpha * [HLL flux with the dissipation term
b_R b_L (E_R - E_L)]`, recovers the diffusion flux only because `E_L, E_R` are the two
cell values.  With PLM face states the jump is O(dx^2), the term vanishes, and the scheme
would under-diffuse grossly (this was wrong in the first version of this note, caught
before 1b was coded).  The blend above is algebraically identical to Berthon's flux for
piecewise-constant states at the diffusion value of `F` (checked: `alpha * dissipation =
(1 - alpha) F_diff` exactly in the isotropic limit), and is second order in both limits:
reconstructed HLL where thin, central compact diffusion flux where thick.
`<rad_m1>/reconstruct = dc | plm` (later `ppm4 | wenoz` from `src/reconstruct/`, which
need `nghost >= 3`) selects the face states; `dc` reproduces the original scheme.

with `f` the face mean, and `(1 - f^2)` Bloch's guard that keeps `f < 1` near free
streaming.  Isotropic limit: `lam = +-c/sqrt3`, `alpha -> 2/(sqrt3 tau_face)`, and the
dissipation term becomes `-(c/(3 tau_face))(E_R - E_L)`: the compact two-point physical
diffusion flux.  Checked here: when the cell fluxes sit at their diffusion value,
`alpha F + (dissipation) = F` EXACTLY at every `tau_face`, i.e. the scheme is the blend
`alpha F_HLL + (1 - alpha) F_diff` of Foucart et al. (2015) with a COMPACT `F_diff`.  This
is what distinguishes it from the Radice/upstream-AthenaK blend, whose thick limit is the
centred average of cell fluxes (a `2 dx` stencil with odd-even decoupling, hence their
sawtooth switch); that variant is not implemented.

Face opacity is the ARITHMETIC mean of `rho kappa_F` (= harmonic mean of the diffusivity,
the correct flux across an opacity jump such as the Fe bump; Bloch eq. 19).  An earlier
draft of this note said harmonic mean of the opacity: that was wrong.

Why not Berthon & Turpault's original (alpha on the WHOLE flux vector, their Sect. 4.3):
it is the variant with the proof (E > 0 and |F| <= cE for CFL <= 1/2, Thm 3.1) and it is
fully explicit, because its interface-based source has the bounded weight
`c sigma/(2 + sigma~ dx)` plus the remainder `(alpha_{i+1/2} - alpha_{i-1/2})/dx f(w_i)`.
But with alpha on `c^2 P` the discrete radiation momentum balance is
`w_eff F = -alpha c^2 grad P`, so a force taken as minus the discrete F-source is
`alpha grad P_rad -> 0` in thick cells.  For radiation hydrodynamics the gas must feel
`-grad P_rad`; with the pressure flux unmodified and the true stiff `rho kappa_F` in the
implicit source (section 4b) the force is conservative and correct.  The price: no
admissibility proof (flux limiting of section 4 stays), and Bloch report spurious flux
oscillations with their well-balanced source near free streaming (we use the cell-centred
implicit source instead, which they recommend for the moving case anyway).

`thick_flux = scaled`.  The modified HLLE of Jiang (2021, App. A) and AREPO-IDORT carried
over to moments.  Their intensity flux has signal speeds
`S = c |mu| sqrt[(1 - exp(-tau_c^2))/tau_c^2]` (downwind side `exp(-tau_c^4)`), with
`tau_c = alpha_J (rho_L + rho_R)(kappa_L + kappa_R) dx`, `alpha_J = 5` (= 20 tau_cell;
IDORT: `5 rho kappa dR`).  Its angular moments are an HLL for `(E, F)` with BOTH wave
speeds scaled by `eps(tau) ~ 1/tau_c`: the VETTAM / KORAL (`a -> min(a, 4/(3 tau))`)
family.  What it is: the physical diffusion comes from the centred average of the slaved
cell fluxes (`2 dx` stencil) and the scaled HLLE term adds `D_num = S dx/2` on top, a
FIXED fraction of the physical diffusivity, independent of `tau`:
`D_num/D_phys = 3<|mu|>/(2 * 20) = 0.04` for Jiang's prefactor, `~0.15-0.25` for IDORT's.
It does not vanish as `tau -> inf` (not AP in the strict sense; neither paper claims a
proof), but it is compact, so unlike a pure centred flux it damps the odd-even mode, at
4-20 % of the physical rate.  What does not carry over: the per-angle upwinding and the
`tau^2`/`tau^4` asymmetry, and the positive/negative coefficient split and first-order
reconstruction (those serve the implicit matrix only).  Moment form, prefactor exposed as
`<rad_m1>/scaled_prefactor` (default 20):

```
eps = sqrt[ (1 - exp(-tau_c^2)) / tau_c^2 ],   b_L,R -> eps * b_L,R   in both fluxes
```

`thick_flux = none`: plain HLL; beam and shadow tests, and the failing control of T3.
Note QUOKKA runs exactly this (PPM + uncorrected HLL) and obtains its AP property from the
time integrator alone (section 5).  That removes the `c^2 dt` error mode of the time
discretisation, but the HLL dissipation `(c/(2 sqrt3))(E_R - E_L)` is still in the E-flux:
it is small only where the reconstruction makes the face jump small.  At a limiter-clipped
extremum the jump is O(dx) and the numerical flux is `~tau_cell` times the physical one.
Their AP demonstration at `tau_cell = 1e5` is a UNIFORM advected medium.  For turbulent
convection at `tau_cell = 1e3-1e6` this is not acceptable as the default (own assessment,
to be exhibited by T3c).

**Moving fluid.**  Bloch et al. (Sect. 6.2, App. A): the thick limit in a moving fluid is
`F -> -(c/(3 sigma)) grad E + (4/3) E v`, and since `alpha -> 0` switches off the
transport flux that carries `(4/3) E v`, "this scheme does not capture the asymptotic
regime in a moving fluid" (they say the same of Gonzalez et al. 2007 and Berthon &
Turpault 2011).  Hence the split written above: `F0 = F - A`, `A = v E + v.P` per cell;
`alpha` acts on the HLL flux of `F0`; `A_upwind` is the FULL enthalpy flux upwinded with
the hydro face velocity from the reconstructed states.  Published relatives: Jiang (2021,
eqs. 16-17) splits `c n I = (c n - f v) I + f v I`, `f = 1 - exp(-tau_c^2)`, and upwinds
the second piece (moments: `f v E` and `f v (x) F`, i.e. `v E`, not the enthalpy flux:
with a scaled HLLE nothing is switched off, so the split there only protects the advected
piece from the HLLE dissipation); Rosdahl & Teyssier (2015) move trapped photons into a
gas-advected variable (rejected here: AREPO-RT's objection that ~60 % of the flux is
already trapped at `tau_c ~ 1` and the radiation pressure becomes isotropic); QUOKKA does
not split (`F = (4/3) v E` is a fixed point of its implicit source); AREPO-IDORT rejects
Jiang's split ("unstable behavior with local time-stepping") and lets the moving mesh
carry the advection.  We sub-cycle with a global step, not local time-stepping, but T4
runs with and without sub-cycling for that reason.  No split is applied in the F
equation: `F` is slaved by the stiff source to `A - c grad.P/(rho kappa_F)` (section 4b
is written in terms of `F0`), and advection of radiation momentum is O(beta^2).  With the
light-speed CFL the advective Courant number is `~beta`, so even a centred `A` would be
stable in the explicit scheme; upwinding is chosen because stage 3 (implicit) needs it.

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

Each substep is IMEX PD-ARS (Chu et al. 2019) as used by QUOKKA (He, Wibking & Krumholz
2024, eqs. 27-29, read), free parameter `eps = 0`; `L` = explicit transport, `S` = the
section 4 solve, which is cell-local, so the implicit stages need no communication:

```
U1      = U^n + dt L(U^n)                     + dt S(U1)
U^{n+1} = U^n + dt/2 [ L(U^n) + L(U1) ]       + dt/2 [ S(U1) + S(U^{n+1}) ]
```

Their analysis (eqs. 56-73) is the reason: any scheme with a stage that contains transport
but not the source (SSP-RK2 with the source split off) leaves `F = F^n - (c dt/3) grad E`
in that stage and a spurious diffusion `c^2 dt/3`, dominant as soon as
`dt > 1/(c rho kappa)`, i.e. always for us (their Marshak front ran 2x too fast).  So the
operator-split "transport, then coupling" of upstream AthenaK is NOT an option, and
section 4 is called inside both stages.  PD-ARS is second order streaming, first order in
the diffusion limit, and was run to radiation CFL ~ 10 in thick problems; with `ap_hll`
the explicit limit stays the light-speed CFL (Berthon & Turpault: 1/2).

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
`2 c/(3 rho kappa)`.  Pass: within 2 % at every `tau_cell` for `ap_hll`; for `scaled`
record the excess against the
predicted `1 + 3<|mu|>/(2 prefactor)` (4 % at the default); second order
between two resolutions at `tau_cell = 10`; control: `none` must FAIL at `tau_cell >= 10`
(if it passes, the test is not testing).  Add a Nyquist-mode seed at `tau_cell = 1e3`:
pass = decay at the physical rate of the Nyquist mode for `ap_hll`; record `scaled`
(prediction: 4 % of that rate).  **(T3b) opacity jump** (Bloch Sect. 5.2): constant-flux
steady state across a 1e3 jump in `rho kappa_F`; pass: flux uniform to 1e-3, no peak at
the jump.  **(T3c) clipped extremum**: a top-hat in `E` at `tau_cell = 1e4` with
`thick_flux = none` and PLM, to exhibit the `~tau_cell` excess flux at limiter-clipped
cells that rules out the uncorrected solver; `ap_hll` must stay within 5 %.

**(T4) Advected thick pulse: dynamic diffusion.**  QUOKKA AP-paper parameters (read):
static `T0 = 1e7, T1 = 2e7, rho = 1.2, w = 24 cm, kappa = 100` (`tau = 2.9e3`,
`beta = 3.3e-5`), 512 cells, their advected-vs-static difference < 0.03 %; dynamic
`v = 3e7 cm/s, kappa = 500` (`tau = 1.4e4`, `beta = 1e-3`, `beta tau = 14`), where they
needed 1024 cells "to reduce the magnitude of odd-even decoupling instability" (< 0.06 %);
we run 512 AND 1024 and require no odd-even mode at 512 with `ap_hll`.  Periodic,
prescribed `v`.  **(T4b)** their Sect. 5.5: uniform medium with `F = (4/3) v E`,
`v = 0.01 c`, `tau_cell = 1e5` (`beta tau = 1e3`); pass: gas temperature drift < 1e-12
(they get 4e-5 with O(v/c) sources, < 1e-15 with O(v^2/c^2)); this is the direct test of
the section 1 source form and of the advective split.  Run T4 with and without
sub-cycling (AREPO-IDORT found Jiang's split unstable under local time-stepping).
Pass: pulse centre at `x0 + v t` within one cell, width equal
to the static case within 2 %, no net gas heating beyond 1e-10 relative (this is the test
that the `beta^2` terms cancel).  Control: the O(v/c)-truncated source must show the
spurious heating.  Run at `chat = c` and once at `chat = c/10` to exhibit the
`(chat/c)(4/3) v` advection error of section 2.

**(T5) Equilibration.**  Single cell, HERACLES numbers (`E_r = 1e12`, `sigma = 4e-8`,
`e_0 = 1e2` and `1e10`), `dt` from 0.1 to 1e5 coupling times.  Pass: correct equilibrium
for every `dt`, no sign oscillation, `E_gas + (c/chat) E` conserved to round-off; repeat
with the tabulated He EOS at (`rho = 1e-8`, `T = 2e5 K`).

**(T6) Marshak wave / Su-Olson.**  Bloch et al. numbers: pass <= 2 % L1 (their AP scheme
1.1 %, uncorrected HLL 84 %).  Also QUOKKA's nonlinear Marshak
(`chi = 300 (kT/keV)^-3` per cm, 60 cells on 0.66 cm, cell optical depth 3 to 3e9):
they get 4.5 % at CFL 0.9; pass: <= 4.5 %.

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

* **R1** (narrowed) thick-limit flux: formulas now read, default chosen (section 3); what
  remains unproven is admissibility of the E-only `alpha` variant (no theorem; Bloch saw
  `f > 1` near free streaming) and its second-order extension (none in print): T1-T3.
* **R2** (narrowed) moving fluid: the enthalpy-flux split in the E equation is our own
  construction (relatives: Jiang 2021, rejected by AREPO-IDORT under local
  time-stepping); T4/T4b decide, with and without sub-cycling.
* **R3** (closed) PD-ARS tableau read; residual: first order in the diffusion limit.
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

## 10. Findings from the implementation (supersede the text above where they differ)

Milestone 1a `6b9a14f9`, 1b `f271d57d`.

* **Section 3, `ap_hll`.**  The identity "`alpha F_HLL` = physical flux" holds for
  piecewise-constant states, so with `reconstruct = dc` the flux is Berthon's
  `alpha F_HLL` and NOTHING is added (adding `(1-alpha) F_diff` double counts: measured
  1.8966 at `tau_cell = 10`, predicted 1.8965).  With `plm` the flux is
  `alpha F_HLL + (1-alpha) F_diff`, and the reconstructed dissipation term inside
  `F_HLL` carries `alpha^2`, not `alpha`: at a kink the limiter returns the cell values,
  the dissipation is then as large as `F_diff` itself, and a weight `alpha` left 65 % of
  the flux as excess at an opacity jump and +2 % in the smooth pulse.  Measured diffusion
  rate / exact: dc 1.000102 / 1.000001 / 1.000000, plm 1.00079 / 1.000015 / 1.000000 at
  `tau_cell` = 10 / 1e3 / 1e6; plm second order (error ratio 5.1 per doubling); Nyquist
  mode decays at 0.99996 of the physical rate; clipped top-hat 1.00005 (plain HLL: 78.7x;
  smooth pulse at `tau_cell` 1e6: 21365x).  To try in 1c: the single form
  `alpha F_HLL + (1-alpha) F_diff * (1 - dE_recon/dE_cell)`, which reduces to both cases
  without a switch.
* **`scaled`** (Jiang / AREPO-IDORT form in moments): +4.3 % with dc at every `tau_cell`,
  as predicted; -1 % with plm; but it damps the Nyquist mode at only 1.5 % of the
  physical rate.  Kept as an option, not the default.
* **Section 5, stability.**  With the source inside both PD-ARS stages the thick regime is
  stable and answer-unchanged to `cfl_rad = 8`; the thin limit degrades at 1.0 and is
  clean to 0.6-0.8.  Sub-cycling `nsub = 10` reproduces `nsub = 1` to 6 digits.
* **Section 4.**  The root find must be `rtsafe` (Newton only when it shrinks the bracket
  faster than bisection): on the tabulated EOS `c_v` is a separate interpolated surface,
  not the slope of its own `e(T)`, and plain Newton stalls.  `EOS_Data::ThermoAt` serves
  both the ideal and the tabulated branch.  Code temperature of the tabulated EOS is not
  kelvin: `arad = a_cgs temp_cgs^4 / pres_cgs`.  The T5 coupling time that probes the
  solver is the thermal one, `rho c_v/(4 c rho kappa a T^3)`, not `1/(c rho kappa)`.
* **Coupling and hydro ghosts.**  After the coupling writes hydro `u0`, a `ConToPrim` is
  not enough: the next hydro fluxes read ghosts exchanged earlier.  The module refreshes
  the hydro ghosts before the inversion (without it: spurious momentum, 2e-10
  conservation break).  Cost item for stage 2: once per substep block may suffice.
* **T3b (opacity jump) as specified fails, and the specification was wrong.**  In steady
  state the FACE energy flux is uniform by conservation and `E(x)` is correct (one 8e-5
  bump); what deviates (0.88 at a 1e3 jump) is the CELL-centred `F` of the two cells
  straddling the jump, because the implicit source divides a centred `2 dx` pressure
  gradient by the cell's own opacity.  The radiation force, `rho kappa F/c`, is then
  exactly the centred `-grad P_rad`, the same discretisation as the gas pressure force:
  correct and conservative.  The cell `F` also enters `alpha avg(F)`, which is negligible
  when the cells are thick and matters only for a jump at `tau_cell ~ 1`.  T3b is
  redefined: `E(x)` against the analytic two-slope solution, for jumps at `tau_cell`
  = 1e3 and ~1; Bloch's interface source stays rejected unless the second case fails.
* `<hydro>/evolution = static` skips the driver loop; tests use a cold dynamic gas or
  `<rad_m1>/gas_feedback = false`.

## 11. Findings from milestone 1c (moving media, and the remaining radiation-only gates)

Milestone 1c: the advective enthalpy-flux split, the unified `ap_hll` form, the O(v/c)
source control, the T4/T4b/T6 problem generators and the redefined T3b.

* **The unified `ap_hll` form was tried and REJECTED as the default, contrary to the
  guess in section 10.**  `F_E = alpha F_HLL + (1-alpha) F_diff (1 - dE_recon/dE_cell)`
  (ratio clamped to [0,1], = 1 where `dE_cell` vanishes) does reduce to Berthon's
  `alpha F_HLL` for `dc` **exactly** (measured: bit-identical T3 rates) and to the blend
  for smooth `plm`, so it removes the reconstruct-dependent switch.  But it is worse
  wherever the two forms differ: T3 `plm` rate error 1.23e-3 vs 7.9e-4 at
  `tau_cell = 10`, resolution-error ratio 3.2 vs 5.1 per doubling (i.e. 1.7th vs
  2.3rd order), Nyquist decay 1.019 vs 0.99996 of the physical rate, T6 L1 0.48 % vs
  0.19 %.  It is better only on the clipped top-hat (1.00000 vs 1.00005) and at
  `tau_cell = 1e3`.  `<rad_m1>/ap_form = alpha2` (the 1b pair) is the default;
  `unified` is kept as an option.
* **The advective split is what makes the thick moving limit work, by a factor 1000.**
  T4 dynamic (`beta tau = 14`, 512 cells): with `advect_split = true` the pulse centre
  lands 0.07 cells from `x0 + v t`; with it off, 69.9 cells short of it (half the
  advection distance) and 17 % too wide.  With `thick_flux = none` the split changes
  nothing (0.053 vs 0.086 cells), which is the direct confirmation of Bloch et al.'s
  mechanism: it is `alpha -> 0` switching off the transport flux that carries
  `(4/3) E v`, not the HLL dissipation.  `scaled` reproduces `ap_hll` with the split.
* **The O(v/c) control.**  With a SINGLE opacity (`kappa_E = kappa_F`) the `beta^2 E`
  and `beta.P.beta` terms cancel analytically between `E0` and `beta.g`, but NOT in the
  discretisation: `E0` enters the implicit backward-Euler energy solve (bounded), while
  `beta.g` is the explicit work term (unbounded).  `source_form = ovc` therefore leaves
  exactly the design's spurious heating `(4/3) c rho kappa beta^2 E`.  T4b measures
  `d ln T_gas` = 0.735 in ONE radiation substep against the predicted 0.634 (16 %,
  the excess from the two source applications of PD-ARS); the full form gives 0.0
  exactly.  At `beta tau_cell = 1e3` the O(v/c) truncation is not a small error.
* **T4b cannot see the split**: the medium is uniform, so `div F = 0` whichever way the
  flux is assembled, and split and no-split are identical to the last digit.  T4b tests
  the SOURCE form; T4 tests the split.
* **T4b must be initialised at the exact fixed point** (`F0_i = 0` and `E0 = arad T^4`,
  i.e. `f = beta(1 + chi(f))` solved for and `E = arad T^4/[(1+beta^2) - 2 beta f +
  beta^2 chi]`, which is `arad T^4 (1 + (4/3) beta^2)` to leading order).  Setting
  `E = arad T^4` with `F = (4/3) v E` instead leaves a REAL O(beta^2) relaxation that
  has nothing to do with the scheme.
* **T3b, redefined, does not pass at 1 %.**  Slope RATIO across the jump (which isolates
  the jump treatment from the overall flux the Dirichlet ends settle on): 2.2e-5
  (`tau_cell` 1 -> 1e3) and 5.4e-5 (1e-3 -> 1) -- excellent.  But the ABSOLUTE slopes are
  1.40 % (first case) and 0.33 % (second) low, because the one-cell `E` bump on the thin
  side of the jump acts as an extra resistance and depresses the steady flux by the same
  amount.  The face flux reconstructed from `E` is uniform to 4-11 % away from the jump
  and has an O(20) spike on the face next to it.  The 1b claim "E(x) is correct, one
  8e-5 bump" was measured at `tau_cell` 0.01 -> 10, where the thin side carries 100x more
  optical depth per unit length; at `tau_cell ~ 1` the same bump costs 1.4 % of the flux.
  `thick_flux = none` fails outright (slope ratio 0.81).  Bloch's interface source is
  therefore NOT clearly refused any more: this is the case that was supposed to decide.
* **T6 is the constant-c_v Marshak wave, not Su-Olson.**  `c_v = alpha T^3` needs
  `e ~ T^4`, which no EOS here has; `t6_marshak.py --marshak` runs its own S_N reference
  with the same constant `c_v` and an incident isotropic bath.  L1(E) at `tau_cell` 2.3:
  0.19 % (`alpha2`), 0.48 % (`unified`), 0.22 % (`scaled`), 1.76 % (`none`); at
  `tau_cell` 4.7: 1.5 / 1.9 / 0.67 / 5.1 %; at 9.4: 4.3 / 4.8 / 2.8 / 16 %.  `none` is
  the failing control from `tau_cell` ~ 4 up.  The reference's own self-error is 0.15 %.
* **The bin writer is single precision** and several of these gates need more: the T3b
  thin side changes `E` by 3e-6 per cell and the T4b drift target is 1e-12.  Those gates
  read `file_type = tab` with `data_format = %26.17e`.  The apparent `|f| = 1 + 7e-8`
  admissibility violation of the beam dump is float32 rounding: in double precision the
  same run gives `max |f| = 0.9999992`.

## References

Audit et al. 2002 (astro-ph/0206281); Berthon, Charrier & Dubroca 2007 (J. Sci. Comput.
31, 347); Berthon & Turpault 2011 (Numer. Methods PDE 27, 1396); Chu et al. 2019
(PD-ARS); Gonzalez, Audit & Huynh 2007 (A&A 464, 429); Krumholz, Klein & McKee 2007
(astro-ph/0611003); Lowrie, Morel & Hittinger 1999 (ApJ 521, 432); Jiang, Stone & Davis
2012 (1201.2223); Skinner & Ostriker 2013 (1306.0010); Sadowski et al. 2013 (1212.5050);
Rosdahl & Teyssier 2015 (1411.6440); Foucart et al. 2015 (1502.04146); Kannan et al.
2019 (1804.01987); Melon Fuksman & Mignone 2019 (1903.10456), Melon Fuksman et al. 2021
(2005.01785); Jiang 2021 (2102.02212); Bloch et al. 2021 (2011.13926); Wibking &
Krumholz 2022 (2110.01792); Menon et al. 2022 (2202.08778); Radice et al. 2022
(2111.14858); He, Wibking & Krumholz 2024 (2404.08247, 2407.18304); Ma, Pakmor, Justham
& de Mink 2025 (2503.16627).
