---
name: general-eos-stage3-table
description: "Stage 3 of the general EOS: the tabulated analytic EOS (H2 + Saha + radiation) — implemented, verified and committed as 03febafa"
metadata: 
  node_type: memory
  type: project
  originSessionId: 9b922d54-af85-46a0-8a62-d8857d1141fd
  modified: 2026-08-12T15:33:16.721Z
---

Done 2026-08-12, committed as **`03febafa`** on `general-eos` (see [[general-eos-project]]).
Fills in all ten `TODO(stage3)` branches in `src/eos/eos.hpp`. ~1050 lines, 3 new files.

## The three decisions the user made (do not relitigate)

1. **Tabulated** `(log10 rho, log10 T)`, not Saha-per-cell and not fitted fractions.
2. **Radiation behind a runtime switch** (`eos_radiation`), default off.
3. **Keep a gamma-law mode** so the existing `general == ideal` regression tests survive.

## Files

- `src/eos/eos_composition.hpp` — the analytic (rho,T) physics, **HOST ONLY**, evaluated
  only at table build. H2/H/H+/He/He+/He++/electrons + inert metals, closed by charge
  neutrality as one 1D bisection on n_e. e = translational + H2 rot/vib (partition
  functions, ortho:para 3:1, Mulholland expansion above 10*theta_rot) + chemical, zero
  point = ground-state H2 + neutral He. NOT included: Coulomb, degeneracy, pressure
  ionization, excited states.
- `src/eos/eos_table.{hpp,cpp}` — the table + device interpolation + the three inversions.

## Three design points that matter

- **Only the GAS is tabulated**; radiation is added analytically in `Eval()`. Makes the
  switch exact rather than a second table, and keeps the surfaces smooth (aT^4 spans 16
  decades across the grid and would dominate the error everywhere else).
- **Bicubic Hermite from node values AND node derivatives**, the derivatives central
  differenced from the ANALYTIC MODEL (h1=1e-5 first, h2=1e-3 cross), not from the table.
  So `chi_rho`/`chi_T`/`c_v` are the analytic derivatives OF THE INTERPOLANT and
  `Gamma_1 = chi_rho + p chi_T^2/(rho T c_v)` holds on the surface the solver actually
  uses. This is the property the WB scheme and Riemann solvers need; matching the true EOS
  is a separate, weaker requirement.
- **Anisotropic resolution, 5x finer in T** (`eos_dlogd`=0.05, `eos_dlogt`=0.01). Not a
  tuning knob: an ionization front is exp(-chi/kT), width in ln T is kT/chi ~ 0.01 dex at
  low density; nothing varies fast in density. Uniform-in-both wastes nodes AND still
  fails. Default table 321 x 651 = **19 MB, ~5 s to build**, single threaded.

## Two bugs found the hard way — remember these

1. **Newton limit cycle.** A bracket-only safeguard is NOT enough. Across an ionization
   transition `dln e/dln T` rises then falls, and Newton from either side overshoots to
   the other while STAYING INSIDE the bracket — an endless ping-pong (observed: y bouncing
   between 2.96 and 4.31 forever, root at 3.545, T wrong by 5x). Fixed with the classical
   `rtsafe` test: bisect if the step leaves the bracket **or** fails to halve the interval
   (`fabs(2*g) > fabs(dzold*dg)`). All three inversions share `SolveLog()`.
2. **`std::max(x, NaN)` returns x**, silently. My Gamma_1 test reported a perfect
   `0.000e+00` over 3000 samples while every sample was NaN. Always guard `if (!(r >= 0))`.
   The underlying test bug: `de = (p/rho^2) drho` is for the **SPECIFIC** energy `u = e/rho`;
   applying it to the volumetric `e` gives garbage.

## Verification (all passing)

Standalone harness lives in **`/orion/u/jinma/ATHENAK/eostest/tabletest/`** (`make &&
./test_table`, prints ALL CHECKS PASSED). Its `stub/athena.hpp` + `stub/parameter_input.hpp`
let `src/eos/eos_table.cpp` compile VERBATIM on the host without linking Kokkos, so the
tests exercise the REAL interpolation and root finds rather than a reimplementation. Also
there: `test_comp` (physics only), `map_err` (error map over the plane, takes a dlog
argument), `diag_g1` (Gamma_1 vs adiabat), `analyze` (turns "x rho e" triples into
T/p/Gamma_1/mu — how the shock tube output was checked). Results:
- interpolation vs analytic model: worst 2.2e-5 on the ionization ridge, 1e-7..1e-12 typical
- derivatives consistent with the interpolant to FD truncation (chi_rho 2.8e-10)
- **Gamma_1 identity vs an independent numerical adiabatic compression: 3.6e-8**
- all three inversions round trip to ~1e-14 from a COLD start
- gamma-law limit (`eos_h2=false`, `eos_ionization=false`): exact to 1.8e-10 — this is the
  strongest check, the whole pipeline against a closed form
- radiation: Gamma_1 -> 5/3 gas dominated, -> 4/3 radiation dominated, exactly

In-code: **all 16 eostest cases still reproduce `eos_compare_baseline.txt` byte for byte**
and `dhj` is bitwise identical, because `general_eos` defaults to `gamma`.

End-to-end run: `eostest/tabsod_gen.athinput` (new), a cgs shock tube through the H
ionization zone. Across the contact rho/e/T all jump but **p is continuous to 5e-4**,
Gamma_1 varies 1.24..1.53, mu 1.22..1.26. Physically correct.
Physics spot checks: mu -> 2.33077 molecular / 0.60548 fully ionized (both match the
analytic limits exactly); Gamma_1 dips to **1.14** in the H2 dissociation zone and 1.27 in
the H ionization zone — the classic structure.

## Input parameters (all in the `<hydro>`/`<mhd>` block)

`general_eos` = gamma|table; `eos_xh`, `eos_yhe`, `eos_a_metal`, `eos_h2`,
`eos_ionization`, `eos_radiation`; `eos_logd_min/max`, `eos_logt_min/max`, `eos_dlog`,
`eos_dlogd`, `eos_dlogt`.

## State of the tree

Working tree clean, `03febafa` is the branch tip and is **NOT yet pushed** to
`origin/general-eos` (everything up to `7bcedea7` is). Build dir is configured
`PROBLEM=built_in_pgens`.

## What is NOT done yet

1. **Built-in test pgens still do `e = p/(gamma-1)` unconditionally**
   (`src/pgen/tests/shock_tube.cpp:94,111,184,204`, and linear_wave likewise). Under
   `general_eos=table` that means `<problem>/pl` is an ENERGY proxy, not the pressure.
   Should be routed through `EintFromP` in `pgen_eos_utils.hpp` (which keeps the ideal
   path bitwise).
2. **No regression test for table mode.** The existing ones only cover `gamma`.
3. Remaining `TODO(stage3)`: the entropy floor in `general_c2p_{hyd,mhd}.hpp`, the
   `EnergyFromPressure` re-solve, and `wb_background.hpp:253`.
4. Table build is ~5 s serial; could be parallelised if it ever matters.
5. No wiki documentation yet (required before merge).

## Validity domain: envelope/atmosphere only (checked 2026-08-15)

Evaluated `EOSCompositionModel` directly against analytic limits. It hits both exactly:
`mu = 0.60548` fully ionized and `mu = 2.33077` all-H2. But:

- **Metals are permanently INERT** (`a_metal = 16`, one particle per nucleus, never an
  electron). Fully ionized they should give 1 + A/2 = 9, so the model's hot-limit mu is
  **0.40% high** at every temperature — 0.60548 against the correct 0.60304. It also means
  no metal electrons, and in the real photosphere the low-chi metals (Na, Mg, Al, Fe, Ca)
  supply ~8x MORE free electrons than hydrogen does.
- **No pressure ionization, no Coulomb/degeneracy corrections.** At high density Saha's
  n_e denominator drives runaway recombination. Error in mu against the model's own ideal
  limit, along a solar track: CZ base (rho 0.2, T 2.2e6) +0.7%; 0.5 Rsun +2.6%;
  0.3 Rsun **+40%**; centre (rho 150, T 1.57e7) **+129%**, where it claims hydrogen is 85%
  neutral and 83% of it **molecular**. T needed for mu within 1%: 6.3e4 K at rho 1e-6,
  1.2e6 K at 1e-1, 5.4e6 K at rho 1.

**DECISION (user, 2026-08-15): leave both as they are.** Correct for the intended use —
the solar runs span rho 1.7e-5 down to 8e-11, nine decades below where Saha breaks, and
in the envelope metals are only singly ionized in reality anyway, so the model's error is
one particle per metal nucleus out of a 8.4e-4 nucleon budget, ~0.1% in mu. That is below
the single-precision `.bin` floor. Revisit ONLY if (a) the grid is reused for rho > ~0.1,
or (b) `get_kapr` is replaced by a real n_e-based H- opacity, which needs the metal
electrons. See [[solar-convection-general-eos]].
