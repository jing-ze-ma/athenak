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
`<rad_m1>/reconstruct = dc | plm | ppm4 | ppmx | wenoz` (the last three from
`src/reconstruct/`, and they need `nghost >= 3`) selects the face states; `dc` reproduces
the original scheme, and only `dc` selects Berthon's form of the E-flux -- every
reconstruction that puts a second-order or better polynomial on the face makes the face
jump O(dx^2) and needs the `alpha2` blend, exactly as plm does (milestone 1c-B).

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
`restart.cpp` and the mirror in `pgen.cpp` (same order in writer and reader; upstream
has no restart support at all).  DONE in milestone 1c-B: the four moments are written as
one more cell-centred block, immediately after the `<radiation>` intensities and before
the turbulence forcing, in both size sums and both loops, for the single-file and the
one-file-per-rank layouts, and behind `pradm1 != nullptr` throughout, so a restart file
written by a run without `<rad_m1>` is byte-identical to what it was before (verified).
Nothing else of the module is written: the PD-ARS stage-1 gas increment `ugas1` is
intra-step state and the `EventCounters`-style solve counters are diagnostics.

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

## 12. Findings from milestone 1c-B (restart, higher-order reconstruction, T8)

Milestone 1c-B: restart support for `<rad_m1>`, the four high-order reconstructions,
the reconstructed velocity in the advective split, the `advect_shear` test and T8.

* **Restart.**  The two fatals are gone and the four moments travel in the restart file
  (section 6).  A `.rst` written by a run WITHOUT `<rad_m1>` is byte-identical to one
  written by the pre-change binary (box_convection, general EOS, 4 files).  T8a is
  bitwise for T4 dynamic (coupled, moving, sub-cycled: both `tab` dumps and the hydro
  `hst` identical to the last digit) and for the 2-D beam.
* **Test problem generators must run on a restart.**  `RadiationM1Tests` used to return
  at once when `restart` was set, so a test with user BCs or a user source term died in
  `pgen.cpp`'s "not enrolled during restart" check: the FIELDS come back from the file,
  the function POINTERS do not.  Each branch now does its parameter reads and its
  enrolment and returns just before writing the arrays.
* **A DEFAULTED Real is serialised into the restart file with six digits.**  The beam
  restart was not bitwise (3e-6 in E) until `beam_angle` was written into the input file
  at full precision: `GetOrAddReal` records its default in the embedded parameter dump
  as `0.785398`, and the restarted run therefore builds a slightly different beam
  direction.  This is a code-wide property of `ParameterInput`, not a `<rad_m1>` bug,
  but any pgen whose state is DERIVED from a defaulted Real is not bitwise restartable.
* **Higher-order reconstruction: worth having, not worth defaulting to.**  pulse1d
  (free streaming, the only exactly solvable M1 case) L1(E) at 512 cells and the orders
  over N = 32..512: plm 2.81e-3 (1.42/1.63/1.76/1.88), ppm4 7.25e-3
  (1.19/1.42/1.10/0.83), ppmx 4.33e-4 (2.28/1.80/1.95/1.99), wenoz 4.34e-4
  (1.69/2.00/2.00/2.00).  ppmx and wenoz are 6.5x more accurate than plm at 512 and
  reach the second order that the RK2 time integration caps them at; **ppm4 is worse
  than plm and loses order** -- its Colella-Woodward limiter clips the reconstructed
  `(E, f_i)` pair hard enough to destroy the convergence, and it should not be used.
  The beam FWHM is 11.38 / 11.32 / 11.38 / 11.31 cells (plm / ppm4 / ppmx / wenoz)
  against the 11.3-cell exact width, `min(E) > 0` in every case, and the apparent
  `max|f| = 1 + 4e-8 .. 8e-8` is the float32 bin dump as in section 11.
* **...but the AP gate does NOT improve, and at high resolution it degrades.**  T3
  `d(sigma^2)/dt` / exact at `tau_cell` = 10 / 1e3 / 1e6: plm 1.000791 / 1.000015 /
  1.000000, ppm4 0.999623 / 0.999999 / 1.000000, ppmx 0.999454 / 0.999997 / 1.000000,
  wenoz 0.999456 / 0.999997 / 1.000000.  All pass the 2 % gate, but the resolution
  sequence at fixed `tau_cell` = 10 (N = 64..512, error `|ratio - 1|`) is plm
  3.45e-3 / 7.91e-4 / 1.56e-4 / 1.4e-5 (order 2.1-3.5) against wenoz
  1.13e-3 / 5.44e-4 / 1.97e-4 / 7.7e-5 (order ~1.4); the high-order methods start
  better and END WORSE, converging to a floor that plm passes through.  The Nyquist
  decay rate at `tau_cell` = 1e3 is 0.99996 (plm) against 0.99973 (ppm4), 0.99971
  (ppmx), 0.99948 (wenoz).  The reason is that the `alpha2` thick-limit blend is built
  around the assumption that the reconstructed face jump is O(dx^2): a 5-cell stencil
  changes the residual dissipation term, not the physical `F_diff`.
* **Cost per cell-update relative to plm** (serial, 256 cells): pure transport
  (`thick_flux = none`) ppm4 1.64x, ppmx 1.71x, wenoz 1.24x; the production
  configuration (`ap_hll` + the implicit coupling) ppm4 1.39x, ppmx 1.36x, wenoz 1.17x.
* **DECISION: plm stays the default.**  wenoz is the one to reach for when free
  streaming accuracy is what matters (the beam, shadows, a thin atmosphere) and costs
  17-24 %; it must not be used where the thick limit is the point, and it needs
  `nghost = 3`.  ppm4 is a trap.  Two risks that the gates do not close: near `f -> 1`
  the 5-cell kernels reconstruct `f_i` components whose rebuilt `|f|` is then clipped by
  `M1Rebuild`, i.e. the limiting happens AFTER the interpolation and the effective
  reconstruction of the direction is not monotone; and at a thick/thin transition the
  face state jump is O(dx) for any method, so the `alpha2` dissipation weight is doing
  the work and the extra order buys nothing (T3b is unchanged by the reconstruction).
* **The velocity in the advective split: `split_vel = recon` is the new default and it
  IMPROVES T4.**  It reconstructs `v` to the face with the same method as `(E, f_i)`
  and upwinds on the mean of the two reconstructed face-normal velocities; `cell` is the
  1c-A behaviour and is bit-identical to the 1c-A binary (checked: the whole flux-kernel
  rewrite is bit-preserving on the `dc`/`plm` path).  On T4 dynamic at 512 cells the two
  differ by L1 1.8e-3 of E -- NOT round-off, because hydro gives the pulse a real
  velocity structure -- and `recon` is the better one on every T4 metric.  At 512 cells:
  pulse centre 0.059 vs 0.074 cells from `x0 + v t`, advected-vs-static difference
  2.356e-3 vs 3.785e-3, width error 1.014e-2 vs 1.015e-2.  At 1024 cells: 0.049 vs
  0.064 cells, 9.404e-4 vs 1.665e-3, width 1.000e-2 both.  The QUOKKA metric is 1.6-1.8x
  smaller with `recon` at both resolutions and converges slightly faster (order 1.33 vs
  1.18), so this is a genuine accuracy gain, not a reshuffling.
* **`advect_shear` (T4c) exists but is not the discriminating test it was meant to be.**
  The T4 dynamic pulse in `v(x) = v0[1 + 0.2 sin(2 pi (x - x1min)/L)]`, prescribed and
  re-imposed every hydro stage by a user source term with `gas_feedback = false`,
  converges against a 2048-cell reference with L1(E) 4.408e-3 / 1.059e-3 / 3.289e-4 /
  1.221e-4 at N = 64/128/256/512 (orders 2.06/1.69/1.43, the tail limited by the
  reference's own error), and `recon` and `cell` agree to three digits in that L1; the two
  differ by only 7.7e-6 / 2.2e-6 / 5.8e-7 / 1.5e-7 of L1(E) at N = 64..512, i.e. second
  order small.  A SMOOTH prescribed shear is the case in which the two velocities differ
  by a clean O(dx) that largely cancels between `F0 = F - A` and the `A` added back
  (the residual weight is `-b_L/(b_R - b_L)`); what T4 has and T4c does not is a
  velocity field with structure on the scale of the pulse itself.  T4c is kept as the
  convergence test of the split in a non-uniform flow; T4 stays the gate that decides
  the option.
* **The stale-hydro-ghost gap does not exist.**  T4 dynamic at 512 cells on 1, 2 and 4
  MeshBlocks is BITWISE identical in the `tab` dumps with `gas_feedback = true`, with
  `gas_feedback = false` and with `coupling = false`.  When the module does not write
  hydro `u0`, hydro's own stage chain has already exchanged its ghosts; when it does,
  `HydroConToPrim` refreshes them (the 1b finding).  No code change was needed.
* **T8b.**  1 rank vs 4 ranks bitwise for T4 dynamic on 4 MeshBlocks (`tab`, 17 digits)
  and for the 2-D beam on 4 MeshBlocks (`bin`).  The ROMIO limitation stands: the output
  directories must pre-exist and must not be on `/tmp`.
* **T8c.**  Over exactly 200 cycles of T4 dynamic (512 cells), `max |d(E_gas +
  (c/chat) E)/E| = 1.8e-13` and `max |d(rho v + F/(chat c))|/|.| = 1.1e-13`, i.e.
  ~1e-15 per cycle: round-off, accumulated as a drift.  Over the full 4.8e-6 s run
  (770 cycles) the two are 3.5e-13 and 3.2e-13.

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

## 13. Findings from milestone 1d (T2 shadow, T7 radiative shocks, T3b closed)

Milestone 1d closes the three items left open at the end of stage 1: the shadow test
(T2), the coupled radiative shocks (T7) and the T3b flux deficit.  New code:
`src/pgen/tests/rad_m1_tests2.cpp` (the `shadow` and `radshock` generators, dispatched
from the existing `rad_m1_tests.cpp` entry), `<rad_m1>/f_source = cell | wb`, and
`tests_m1/t2_shadow.py`.  The default path is BITWISE unchanged: 68/68 dump files of
T3 (`tau_cell` 1e3), T3c (top-hat), T4 dynamic 512, T4b, T6 (128) and T3b case A are
byte-identical to the pristine `rt-integration` HEAD binary run on the same inputs.

### 13.1 T3b: the flux deficit was mostly the TEST, and `wb` fixes the rest

Section 11 attributed the "absolute slopes 1.40 % (case A) / 0.33 % (case B) low" to a
one-cell `E` bump at the jump acting as an extra resistance.  **That attribution is
retracted.**  The decisive control is the SAME run with `jump_ratio = 1`, i.e. a uniform
medium at `tau_cell = 1` with the same Dirichlet ends: it shows a flux deficit of
0.676 % at 64 cells and 0.282 % at 128, first order in `dx` and with no jump anywhere.
The cause is in `RadM1FixedEBC`: it filled every ghost cell with the end value `E(x1min)`,
i.e. it imposed the boundary FACE value at the ghost cell CENTRE.  That half-cell offset
removes `dE_total/nx1` of driving from a problem whose flux is fixed by the total `dE`
across the domain.  `<problem>/exact_ghost = true` (default false, so nothing moves
without it) continues the analytic slope into the ghosts; with it the uniform case is
EXACT to 1e-11 in every metric -- the scheme is already well balanced on a constant
`rho kappa_F` -- and what is left in the jump runs is the jump error alone.

With `exact_ghost` on, and the jump window widened to 6 cells so that the script's
reconstructed face flux is not sampled inside the kink (steady state, `t = 1e5` for
case A and `2000` for case B, 64 cells):

| | A `cell` | A `wb` | B `cell` | B `wb` |
| --- | --- | --- | --- | --- |
| absolute flux error | 1.346e-3 | **1.780e-5** | 3.467e-3 | **4.892e-4** |
| `max｜E/E_exact-1｜` | 1.288e-4 | **3.00e-6** | 3.094e-4 | 2.539e-4 |
| slope ratio error | 2.37e-7 | 6.89e-8 | 2.30e-6 | 4.31e-6 |
| cell-centred `F` error | 43.25 | **0.656** | 0.1287 | 0.0971 |
| face-flux peak at the jump | 18.03 | **0.881** | 18.15 | 77.1 |

So the < 0.2 % gate is met by `wb` in both cases (0.0018 % and 0.049 %) and by the
default only in case A (0.135 %; case B is 0.347 %).

**The diagnosis.**  Mechanism (i) of the brief is the real one and `wb` removes it.  In a
constant-flux steady state the implicit `F` source balances the discrete centred
`2 dx` pressure gradient against the cell's OWN `rho kappa_F`, and
`-(c/3)(E_{i+1}-E_{i-1})/(2dx) = (1/2)[sigma_{i-1/2} + sigma_{i+1/2}] F0`, so the cell
`F` comes out `F0 (sigma_{i-1} + 2 sigma_i + sigma_{i+1})/(4 sigma_i)` -- 250x too large
in the last thin cell of case A (measured 43x after the neighbouring faces have reacted).
That wrong `F` enters the E-flux through `alpha*avg(F0)` on the neighbouring face, where
`tau_face ~ 1` and `alpha ~ 0.5`.  `f_source = wb` divides by the Bloch et al. (2021)
eq. 18 trapezoidal interface value instead, `(1/2)[sigma_{i-1/2} + sigma_{i+1/2}]` with
the same ARITHMETIC face mean the thick-limit flux uses, i.e. the 1-2-1 average
`(sigma_{i-1} + 2 sigma_i + sigma_{i+1})/4`, per direction.  The analytic two-slope
solution is then an exact discrete steady state of BOTH the `E` and the `F` equation,
which is why case A collapses to 1.8e-5.  Mechanism (ii) (the arithmetic face mean in
`alpha`) is NOT a defect: for the two-slope solution the arithmetic mean of `rho kappa`
is exactly the right face resistance, and `F_diff` built with it reproduces the analytic
`dE` at every face including the jump face to round-off.

What survives in case B is a third thing, neither (i) nor (ii): at `tau_cell = 1e-3 -> 1`
the jump face has `tau_face ~ 0.5` and `alpha^2 ~ 0.25`, and the PLM limiter CLIPS the
reconstruction on the thick side, leaving a face jump `E_R - E_L` of order half the cell
step.  The residual HLL dissipation `alpha^2 b_+ b_- (E_R - E_L)/(b_+ - b_-)` is then a
real O(10 %) flux error at that one face, which the steady state pays for with a slightly
depressed global flux.  `alpha^2` only kills it where `alpha` is small.  No remedy for
that is shipped: it needs a reconstruction that knows about the opacity kink, which is
outside the scope of stage 1.

**The force on the gas.**  `wb` does not touch the momentum exchange, which is
`delta(rho v) = -(F' - F*)/(chat c)` by construction and therefore conservative to
round-off whatever opacity the `F` update uses.  In steady state the force per unit time
is `-dP_rad/dx` with the same centred discretisation as the gas pressure force, exactly as
with `cell` (the opacity cancels: `F sigma_eff = -c dP/dx` gives
`delta(rho v)/dt = -dP_rad/dx`).  Out of steady state the two cells straddling a jump
exchange momentum at the interface-averaged opacity rather than at their own, which is
the same averaging the flux already uses.

**Regression (default vs `wb`, plm, otherwise stage-1 defaults).**

| gate | `cell` | `wb` |
| --- | --- | --- |
| T3 `d(sigma^2)/dt`/exact, `tau_cell` 10 | 1.000791 | 1.000791 |
| ... `tau_cell` 1e3 | 1.000015 | 1.000015 |
| ... `tau_cell` 1e6 | 1.000000 | 1.000000 |
| T3 Nyquist decay / physical | 0.999961 | 0.999961 |
| T3c clipped top-hat | 1.000054 | 1.000054 |
| T4 dynamic 512, `dcentre` (cells) | 0.054 | 0.054 |
| T4 dynamic 512, Nyquist power | 8.323e-14 | 8.321e-14 |
| T4b `dT_gas/T` drift | 0.0 | 0.0 |
| T6 Marshak L1(E), 64 / 128 | 0.0148 / 0.0019 | 0.0148 / 0.0019 |

Every uniform-opacity gate is bit-identical (a 1-2-1 average of a constant IS the
constant); T4, whose `rho kappa_F` varies smoothly across the pulse, differs by 2e-4
relative in the Nyquist power and by nothing else at the printed precision.

**Should `wb` be the default?**  Recommended, but left OFF in this branch.  It is free
where the opacity is smooth (bitwise), it is the only setting that meets the T3b gate in
both cases, and stage 2 runs a tabulated stellar opacity whose Fe bump is exactly the
`tau_cell ~ 1` jump case A models.  Cost: four extra loads of `opac` per cell in the
coupling kernel.  The switch should be flipped together with the stage-2 opacity table,
with the box G1 gate re-run.

### 13.2 T2, the shadow test

`<problem>/m1_test = shadow` (`inputs/tests/rad_m1_shadow.athinput`), the design note's
HERACLES numbers: 1 cm x 0.12 cm on 280 x 80, ambient `rho0 = 1 g/cc` (the design note's
value; the classic Hayes & Norman setup uses 1e-3, but only `rho/rho0` enters the opacity
law so the two are the same problem), an elliptical clump at `(0.5, 0)` with semi-axes
`(0.1, 0.06)` and 1000x the ambient density with a Fermi edge of width `delta = 10`,
`sigma = 0.1 (T/T0)^-3.5 (rho/rho0)^2` per cm as the module's `powerlaw` law
(`kappa0 = 0.1, opac_a = 1, opac_b = -3.5, rho_ref = 1, t_ref = 290`), absorption only
(`kappa_P = kappa_E = kappa_F`, no scattering), `T0 = 290 K` for gas and radiation, and a
source at `T_r = 1740 K` entering the inner-x1 face with `f = 1 - 1e-6`.  Outer-x1 and
outer-x2 are `vacuum`; the inner-x2 face is `reflect`, i.e. the upper half of the
symmetric problem is modelled.  Units are cgs with the gas temperature carried in kelvin
(`e = rho T/(gamma-1)`, `arad = 7.5657e-15`); `gas_feedback = false` and a user source
term re-impose the gas every hydro stage, so the medium is exactly static and the opacity
is a pure function of position -- the gas heat capacity in these units exceeds `a T0^4` by
ten orders of magnitude, so it is a thermostat in any case.  `subcycle = false`, so the
mesh timestep IS the radiation one.

### 13.3 T7, the radiative shocks

`<problem>/m1_test = radshock`, Lowrie & Edwards Mach 2 subcritical and Mach 5
supercritical with the parameter sets printed by `tests_m1/t7_radshock.py --athinput`
(`P0 = 1e-4`, `sigma_a = 1e6`, `kappa = 1`; cgs `rho0 = 5.69`, `T0 = 2.18e6 K`,
`sigma = 577.16` per cm -- 577.35 per `L~ = 1.0003 cm`).  The runs are INITIALISED FROM
THE SEMI-ANALYTIC PROFILE (`--write-ref`, read by the generator and shifted by
`m1_shock_xs`), not from a step: at the light-speed CFL a step start would need several
flow-through times (4.4e-10 s at `dt = 4e-16 s`) to build the precursor, which is
~5e6 substeps on one core.  Both x1 faces are Dirichlet at the far-field states, so the
shock is held in its own frame and the run measures whether the scheme KEEPS the
solution.  Opacity: `powerlaw` with `opac_a = -1`, which stores a constant `rho kappa`
independently of `rho`, as the problem specifies.  Units are honest cgs with the code
temperature `T_code = [k/(mu m_H)] T_kelvin` (`m1_shock_tunit = 8.2499208e7`) and
`arad = a_cgs/tunit^4 = 1.6332399e-46`, which reproduces the `c_v = 1.2374881e8`
the analysis script expects.


**T2 result.**  The front arrives at `x = 0.98` at `3.1893e-11 s` against `x/c =
3.2689e-11 s`, i.e. `-2.44 %`, and a linear fit of the 1 % contour gives a front speed of
`1.0053 c`: PASS on the 5 % criterion.  Taken at the HALF maximum instead, the arrival is
`+8.54 %` late and the fitted speed `0.868 c`; that is the PLM + HLL smearing of the
front over tens of cells, not a wave-speed error, which is why the gate is written on the
leading edge and both numbers are reported.  The shadow depth
`E(0.8, 0)/E(0.8, 0.115)` deepens to `2.1e-3` after one crossing and then FILLS to a
steady plateau of `7.4-7.8e-2` from four crossings on; the final value is `7.727e-2`.
That does NOT meet the design note's `< 1e-3`, and it is not the ambient re-emission
either (`a T0^4` is `8.5e-4` of the lit value).  The `closure = eddington` control gives
`8.858e-1`, i.e. the umbra reaches 89 % of the lit value and the shadow is destroyed:
M1 is 11.5x sharper, which is the qualitative statement the test exists to make, but the
7.7 % residual is a real M1-plus-PLM artefact at this resolution and is recorded as such.
One deviation from the brief: the lit reference is taken at `y = 0.115`, not at the
`y = 0.1` the design note names, because the Fermi edge puts `rho = 2.3 rho0` at
`y = 0.1` and `7.7 rho0` at `y = 0.09`, so `y = 0.1` sits in the halo's penumbra (`E`
there is 0.45 of the free-streaming value) and a ratio formed with it is 0.974.

**T7 result** (`t = 1e-10 s`, L1 gate 2 %, `rho / T_gas / T_rad`):

| run | L1 `rho` | L1 `T_gas` | L1 `T_rad` |
| --- | --- | --- | --- |
| M0=2 M1, 512 | 4.03e-3 | 6.26e-3 | 9.99e-3 |
| M0=2 eddington, 512 | 4.25e-4 | 4.94e-4 | 8.84e-5 |
| M0=2 `chat/c` 0.1 / 0.01 | 4.09e-3 / 4.58e-3 | 6.26e-3 / 6.31e-3 | 9.99e-3 / 1.00e-2 |
| M0=5 M1, 2048 | 2.48e-3 | 5.68e-3 | 7.54e-3 |
| M0=5 M1, 1024 | 2.99e-3 | 5.60e-3 | 7.21e-3 |
| M0=5 eddington, 1024 | 4.05e-4 | 5.21e-4 | 1.31e-4 |
| M0=5 `chat/c` 0.1 / 0.01 | 3.15e-3 / 3.64e-3 | 5.62e-3 / 6.22e-3 | 7.23e-3 / 7.79e-3 |

All PASS.  Zel'dovich spike and precursor length against the reference: M0=2 M1
`-0.25 %` and `+11.4 %`; M0=5 M1 `-2.96 %` (2048) / `-5.01 %` (1024) and `+4.1 %`; the
spike error halves between 1024 and 2048 because the M0=5 relaxation region is only 6
cells wide at 512 and 23 at 2048.  `closure = eddington` is 8-100x more accurate on
every L1 norm, which is the Skinner & Ostriker ordering (they report 1.7/6.1/7.8 % for
M1 with computed eigenvalues against 0.42/0.49/0.42 % for `lambda = +-c/sqrt3`); our M1
numbers are 4-10x better than theirs, but the ordering stands, and it is a diffusion-limit
comparison, not a test of the closure.

**RSLA.**  `chat/c = 0.1` changes nothing; `chat/c = 0.01` grows the errors by 14 %
(M0=2) and 22 % (M0=5) RELATIVE, i.e. 0.40 -> 0.46 % and 0.30 -> 0.36 % in `rho`.  This
is far smaller than PLUTO's `> 40 %` at `c/1000`, and the reason is worth recording: at
`d/dt = 0` the `chat` on the transport term and the `chat` on the source CANCEL, so the
steady state of the moment system is chat-independent.  A shock held in its own frame and
started from the semi-analytic profile can therefore only expose RSLA through the
transient and through terms of order `v/chat`; PLUTO's number comes from `c/1000`, where
the supercritical upstream `v = 8.7e7 cm/s` EXCEEDS `chat = 3e7 cm/s` and the scheme is
not even hyperbolically sensible.  That is a stronger statement of section 2's condition,
not a weaker one: RSLA is safe for a steady state and wrong for everything time
dependent.

**Sub-cycling.**  With `subcycle = false` (the input's default for this test) the mesh
step IS the radiation step.  With `subcycle = true` the hydro step is `2.28e-13 s` against
the radiation `3.99e-16 s`, so `N_sub = 571`, and the three L1 norms agree with the
un-sub-cycled run to 1.6 % of the error itself.
