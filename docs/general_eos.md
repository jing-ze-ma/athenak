# General equation of state (non-relativistic Hydro and MHD)

Replaces the ideal-gas assumption in non-relativistic `hydro` and `mhd` with a general
`p(rho, e)`. The target physics is an **optically thick interior with radiation locked to
the gas (LTE), plus H2 dissociation and H/He ionization** — gas-giant and stellar
envelopes, where the adiabatic exponent `Gamma_1` drops well below 5/3 through the
partial-ionization zones and an ideal gas gets the sound speed and the convective
stability boundary wrong.

The general path is **added alongside** the ideal path, not in place of it. Runs with
`eos = ideal` allocate nothing extra, take the original code path, and are bitwise
unchanged.

---

## Turning it on

```
<hydro>                      # or <mhd>
eos         = general
general_eos = table          # gamma (default) | table

<units>                      # REQUIRED -- see below
length_cgs = 1.0
mass_cgs   = 1.0
time_cgs   = 1.0
mu         = 1.0
```

### `<units>` is required

An ideal gas is scale free; a general EOS is not. Ionization and dissociation depend on
the **absolute** density and temperature, so the run must declare its physical scale.
`<hydro|mhd>/eos = general` aborts without a `<units>` block.

Running in cgs (`length_cgs = mass_cgs = time_cgs = 1.0`) is the recommended setup: the
conversion factors are then all 1 and there is nothing to get wrong.

### The two modes

| `general_eos` | what the interface evaluates |
|---|---|
| `gamma` (default) | the same ideal gamma law as `eos = ideal` |
| `table` | the tabulated EOS: H2 dissociation, H/He Saha ionization, optional radiation |

`gamma` is not really a physics option. It is a **control**: it makes the general code
path — the reconstruction of `p` and `Gamma_1`, the general Riemann solvers, the general
well-balanced background, refinement — exercisable against an answer that is known
exactly. The regression tests use it to assert that the general machinery reproduces the
ideal path run for run. Selecting `eos = general` without `general_eos = table` therefore
buys nothing but cost.

---

## Composition convention: `mu_ref = 1`

Under a general EOS **composition lives in the EOS**, not in the unit system. The code
works internally with `p = rho T / mu(rho,T)` and a reference mean molecular weight of 1.

Consequences, all of which surprise people:

* `<units>/mu` **does not set the composition**. It is divided back out of the temperature
  scale. Set it to 1.0 and forget it.
* A problem generator's `<problem>/Rgas` (or equivalent hardcoded gas constant) is
  **ignored** on the general path. `deep_hot_jupiter_rt` warns if the value it was given
  disagrees with the units block; `solar_convection` and `cooling_convection` hardcode
  `Rgas = 1.38e8` in several places, and under `general_eos = table` the atmosphere they
  build will differ from the ideal-gas one. That is correct behaviour, not a bug.
* Code temperature times `EOS_Data::temp_cgs` gives kelvin. Do **not** use
  `Units::temperature_cgs()`, which folds in the fixed `<units>/mu`; the two agree only
  when `mu = 1`.

---

## Parameter reference

All of these live in the `<hydro>` or `<mhd>` block. They are read only when
`general_eos = table`.

### Composition

| parameter | default | meaning |
|---|---|---|
| `eos_xh` | `0.7381` | hydrogen mass fraction X |
| `eos_yhe` | `0.2485` | helium mass fraction Y |
| `eos_a_metal` | `16.0` | atomic weight of the "metal" component (Z = 1 − X − Y). It contributes particles, and never dissociates. |
| `eos_h2` | `true` | include H2 and its dissociation equilibrium, and its rotational/vibrational energy |
| `eos_ionization` | `true` | include Saha ionization of H, He and He+ |
| `eos_metal_ionization` | `false` | let Na, K, Ca, Al, Mg and Fe ionize — see below |
| `eos_metal_mh` | `<problem>/met`, else `0.0` | [M/H] in dex for those donors |
| `eos_metal_condensation` | `false` | rain each donor out below its own T_cond(p, [M/H]) |
| `eos_metal_tcond` | `0.0` | crude alternative: one temperature for all species; overrides the above when positive |
| `eos_radiation` | `false` | add `a T^4` energy and `a T^4/3` pressure on top of the tabulated gas |

`eos_xh` and `eos_yhe` must be non-negative and sum to at most 1; the run aborts
otherwise. Radiation is added **analytically on top of** the tabulated gas rather than
being tabulated, which makes the switch exact and keeps the tabulated surfaces smooth.

### Metal ionization, and why it exists

`eos_metal_ionization` is about **n_e, not about p or e**. Below ~4000 K hydrogen supplies
essentially no free electrons — 13.6 eV is far too high — and the electrons come from
Na, K, Ca, Al, Mg and Fe at 4–8 eV. Anything that needs an ionization degree, above all
the Ohmic resistivity, is wrong by orders of magnitude in a cool atmosphere without them.

The thermodynamic cost of switching it on is negligible, by construction: the tracked
species are already inside the metal lump and contribute one heavy particle each whether
ionized or not, so the only change to `n_tot` is their electrons. Measured, μ moves by
~0.004% in a molecular atmosphere and ~0.4% in the fully ionized limit. Their ionization
energy is carried too — far too small to matter for p or e, but it is what keeps `c_v`
consistent across the metal ionization zone.

**Metallicity is set in one place.** `eos_metal_mh` defaults to `<problem>/met` when that
parameter exists, because a problem generator with a metallicity uses it for its opacity,
and an atmosphere opaque at one metallicity while conducting at another is not a physical
configuration — it is a parameter that got updated in one place only. Set `eos_metal_mh`
explicitly to override; `deep_hot_jupiter_rt` then aborts if it disagrees with its `met`.

**Condensation.** Below roughly 1000–1800 K the donors condense out and gas-phase
abundances fall by orders of magnitude, taking n_e with them; Saha with full solar
abundances would make those regions far too conductive. `eos_metal_condensation` removes
each species below its own condensation curve, in the standard published form

    10^4 / T_cond = a - b log10(p_T/bar) + d (log10 p_T)^2 - c [M/H]

with p_T the total pressure. The coefficients are taken from the literature, not fitted
here, and give

Every coefficient is from the literature. Per species, so each row can be checked against
its own source:

| element | condensate | published fit | source |
|---|---|---|---|
| Na | Na2S | 10⁴/T = 10.045 − 0.72 log p − 1.08 [Fe/H] | [Morley+2012](https://doi.org/10.1088/0004-637X/756/2/172), Table 1 |
| K | KCl | 10⁴/T = 12.479 − 0.879 log p − 0.879 [Fe/H] | [Morley+2012](https://doi.org/10.1088/0004-637X/756/2/172), Table 1 |
| Ca | CaAl₄O₇ (grossite) | 10⁴/T = 4.990 − 0.2394 log p + 1.398e−3 (log p)² − 0.595 [Fe/H] | [Wakeford+2017](https://doi.org/10.1093/mnras/stw2639) |
| Al | Al₂O₃ (corundum) | 10⁴/T = 5.014 − 0.2179 log p + 2.264e−3 (log p)² − 0.580 [Fe/H] | [Wakeford+2017](https://doi.org/10.1093/mnras/stw2639) |
| Mg | Mg₂SiO₄ (forsterite) | 10⁴/T = 5.89 − 0.37 log p − 0.73 [Fe/H] | [Visscher+2010](https://doi.org/10.1088/0004-637X/716/2/1060), eq 18 |
| Fe | Fe metal | 10⁴/T = 5.44 − 0.48 log p − 0.48 [Fe/H] | [Visscher+2010](https://doi.org/10.1088/0004-637X/716/2/1060), eq 2 |

Pressure is the total pressure in bar. Visscher+2010 state their expressions are valid for
800–2500 K and [Fe/H] ≤ +0.5; the Wakeford+2017 fits hold where that phase is the
highest-temperature Al- or Ca-bearing condensate. Behind Morley+2012 sit the
saturation curves log p'_Na = 8.550 − 13889/T − 0.5[Fe/H] (from
[Visscher+2006](https://doi.org/10.1086/506245)) and log p'_KCl = 7.611 − 11382/T.

Which gives:

| | Ca | Al | Fe | Mg | Na | K |
|---|---|---|---|---|---|---|
| T_cond at 1e-3 bar | 1748 K | 1758 K | 1453 K | 1429 K | 819 K | 662 K |
| T_cond at 1 bar | 2004 K | 1994 K | 1838 K | 1698 K | 996 K | 801 K |
| T_cond at 10 bar | 2104 K | 2084 K | 2016 K | 1812 K | 1072 K | 862 K |

`log10 p_T` is clamped to [−8, +3] before use: the refractory fits are quadratic in it,
and the table grid extends far outside any real atmosphere.

These are the analytic curves the field uses — virga, PICASO, and the Ackerman & Marley
cloud lineage all draw on them. They are not the newest: [virga v1](https://arxiv.org/abs/2508.15102)
moved Fe, Mg and Al onto Morley et al. 2024 saturation pressures, shifting their T_cond by
1–2%, while KCl and Na2S are unchanged from Morley+2012. Since K and Na dominate n_e
wherever the refractories are condensing, that revision does not move the electron budget.
Full equilibrium condensation with rainout ([FastChem Cond](https://arxiv.org/abs/2309.02337),
GGchem) is the accuracy ceiling, but it is a Gibbs minimisation per node rather than a
curve, which is far more machinery than deciding when a donor leaves the gas phase needs.

Two approximations remain, both removing an element slightly too early: Ca is given the
grossite curve, though grossite and hibonite consume only a few tenths of the Ca inventory
with the rest following into Ca-silicates over the next few hundred kelvin; and Mg is
given forsterite, ignoring the enstatite behind it. Neither matters much for the electron
budget, because below ~1700 K — where all four refractories have gone — n_e is dominated
by Na and K, which stay in the gas down to ~800–1000 K.

The transition is a tanh over 5% in temperature rather than a step, because the result is
differenced to build the interpolation table and a discontinuity there would give node
derivatives that make the Hermite patch ring. Only the electron donation is suppressed:
the condensed material stays in the mass and particle budget, since whether it rains out
of the column or remains as cloud is a transport question this local model cannot answer,
and the mass involved is ~1e-4 of the gas.

The net effect on the electron budget is a cliff rather than a gentle taper, because the
species that dominate n_e at low temperature — K, then Na — are also the last to condense.
Measured at ρ = 1e-5 g/cm³, n_e is unchanged above 900 K, falls to 0.17× at 700 K, 0.013×
at 600 K and effectively zero by 400 K. `eos_metal_tcond` remains as a cruder alternative,
a single temperature applied to every species.

Condensation costs a second pass through the composition model, because it needs the
pressure and the pressure is an output. That is host-side table-build code, so it costs
nothing at run time, and it is exact enough because the donors are ~1e-4 of the particles
and p is insensitive to whether they have condensed.

Ground-state statistical weights are used throughout. That is good for the alkalis, whose
first excited states lie ~2 eV up, and poorer for Fe, whose low-lying levels make the true
neutral partition function ~1.5–2× the ground-state 25 near 3000–5000 K — worth tens of
percent in n_e where iron dominates.

### Table grid

The table is in `(log10 rho_cgs, log10 T_K)` and is built once, on the host, at setup.

| parameter | default | meaning |
|---|---|---|
| `eos_logd_min` | `-14.0` | lower bound in log10 density (cgs) |
| `eos_logd_max` | `2.0` | upper bound in log10 density (cgs) |
| `eos_logt_min` | `1.5` | lower bound in log10 temperature (K) |
| `eos_logt_max` | `8.0` | upper bound in log10 temperature (K) |
| `eos_dlog` | `0.05` | default spacing, used for whichever of the two below is not set |
| `eos_dlogd` | `eos_dlog` | node spacing in log10 density |
| `eos_dlogt` | `0.2*eos_dlog` | node spacing in log10 temperature |

**The temperature grid is five times finer by default, and that is physics rather than a
tuning knob.** An ionization or dissociation front is a factor `exp(-chi/kT)`, so its
width in `ln T` is `kT/chi` — about a hundredth of a decade at the low densities where
hydrogen ionizes near 3000 K — while nothing in the EOS varies faster than smoothly in
density. A grid uniform in both logs wastes most of its nodes and still fails to resolve
the only feature that is hard.

The defaults give a 321 x 651 table (about 19 MB). At startup the code prints the grid,
the composition, and `mu` at the grid centre — check that line, it is the cheapest way to
confirm the run is using the EOS you think it is.

Interpolation is bicubic Hermite, with node derivatives obtained by central differencing
the **analytic model** rather than the table, which makes the interpolant fourth-order
accurate rather than second.

### Floors

The usual `dfloor`, `pfloor` and `tfloor` all work. `sfloor` does not — see below.

#### `tfloor` changes meaning, and `tfloor_kelvin` avoids it

`tfloor` is a floor on the **code** temperature, and which quantity that is depends on the
EOS:

| | code temperature | kelvin per unit |
|---|---|---|
| `eos = ideal` | `p/rho` | `Units::temperature_cgs()`, which **includes** `<units>/mu` |
| `eos = general` | the EOS's own `T` | the same factor with that `mu` divided back out, because the general EOS works at `mu_ref = 1` and keeps composition inside the EOS |

The two differ by exactly `<units>/mu`. So a `tfloor` copied unchanged from an ideal-gas
input into a general-EOS one is silently wrong by that factor — and if the ideal input had
no `<units>` block at all, its `tfloor` was in units of the problem's own `Rgas` and the
number means nothing here. This is a real trap: `inputs/mhd/deep_hot_jupiter_rt.athinput`
carried `tfloor = 400.0`, which with `Rgas = 4.593e7` is about `1e-5` K, i.e. no floor.

`<hydro|mhd>/tfloor_kelvin` states the same floor in **kelvin** and converts it at startup,
so an ideal and a general input can carry the identical number:

```
tfloor_kelvin = 200.0     # instead of tfloor = 1.663e10 under a cgs general EOS
```

The conversion is reported at startup. Notes:

- Setting both `tfloor` and `tfloor_kelvin` is a fatal error — they are one floor in two
  units, so there is no sensible precedence.
- It needs a `<units>` block, which the general EOS requires anyway. For an **ideal** run
  the conversion goes through `<units>/mu`, so that must be the mean molecular weight the
  problem actually assumes (for a pgen carrying `Rgas`, `mu = k_B/(m_u Rgas)`).
- Supported for the non-relativistic `ideal` and `general` EOSs only. It is refused under
  `isothermal` (whose C2P has no temperature floor) and under SR/GR (where `tfloor` is not
  this quantity).
- `tfloor` itself is completely unchanged, so no existing input behaves differently.

---

## What is supported, and what is refused

**Supported:** all reconstruction methods (`dc`/`plm`/`ppm4`/`ppmx`/`wenoz`); the
`llf`, `hlle`, `hllc`, `hllclm`, `lhllc` and `ausmpup` hydro solvers and the `llf`,
`hlle` and `hlld` MHD solvers; the CFL timestep; first-order flux correction (FOFC);
static and adaptive mesh refinement; the well-balanced (Kappeli & Mishra) scheme in all
its options; thermal conduction; and the gravity and cooling source terms.

**Refused, loudly:**

| what | why |
|---|---|
| `<hydro>/rsolver = roe` | The Roe average is only defined for an ideal gas. Substituting a local `Gamma_1` would not give a valid Roe average — the linearized matrix would stop satisfying `F_R - F_L = A (U_R - U_L)`, and shock speeds would be silently wrong. A correct general-EOS Roe solver needs Vinokur–Montagne averaging, which is not implemented. (MHD's Roe solver is commented out upstream, so there is nothing to refuse there.) |
| `<hydro\|mhd>/sfloor` with `general_eos = table` | `sfloor` floors the ideal-gas entropy variable `p/rho^gamma`, which is not an invariant of a tabulated EOS. Use `pfloor` or `tfloor`. It works normally under `general_eos = gamma`. |
| `eos = general` with SR/GR | The general EOS is non-relativistic only. Relativistic runs use `eos/primitive-solver/`. |

**Known accuracy caveat:** `hlle` differs from the ideal path at truncation level *by
design*. Its ideal-gas wave speeds use a Roe average with no general-EOS analogue, so the
general path falls back to a Davis/Einfeldt bound. This also shows up at a **polar
boundary in spherical coordinates**, which forces `hlld -> hlle` at the pole.

---

## Cost model

Under a general EOS the physics is naturally a function of `(rho, T)`, but the code stores
`(rho, e)`. Recovering `T` from `(rho, e)` is a root find and is **the** expensive
operation; everything else is cheap once `T` is known.

The code is structured around that:

* `ConsToPrim` solves for `T` **exactly once per cell**, warm-started from the cached `T`
  of the previous stage, and caches it in `Hydro/MHD::wtemp`.
* It also evaluates `p` and `Gamma_1` once per cell into `wder`, which are then
  *reconstructed* to interfaces. **The Riemann solvers never call the EOS.**
* Every EOS accessor comes in two forms: a primary `(d,e,T)` form that costs nothing extra
  when the caller already has `T`, and a convenience `(d,e)` form that solves for `T`
  itself. **Kernels on the hot path must use the `(d,e,T)` form.** Calling the `(d,e)`
  form N times on one cell means N root finds.

If you add code that touches the EOS, that last point is the one to remember.

---

## Validity of the physics model

The composition model covers H2, H, H+, He, He+, He++, free electrons, a metal component
that contributes particles, and — with `eos_metal_ionization` — Na, K, Ca, Al, Mg and Fe
as singly ionizing electron donors. It is closed by charge neutrality.

**Not included:** Coulomb (non-ideal) corrections, electron degeneracy and pressure
ionization, excited bound states beyond ground-state statistical weights, condensation of
the metals, and any second ionization of them.

It is therefore valid for the **weakly coupled, non-degenerate** regime — gas-giant and
stellar envelopes — and not for deep stellar interiors or degenerate objects. Measured
against the analytic fully-ionized limit, the error is +0.7% at the base of the solar
convection zone and +2.6% at 0.5 R☉, but +40% at 0.3 R☉ and +129% at the centre, where
the absence of pressure ionization makes Saha recombine hydrogen that should be ionized.
The temperature at which μ comes within 1% of the ionized limit is 6.3e4 K at
ρ = 1e-6 g/cm³ and 5.4e6 K at ρ = 1, so a table grid extending to high density contains
regions the model cannot describe.

---

## The electron fraction, and Ohmic resistivity

`EOS_Data::ElectronFraction(d, e, T)` returns n_e/n_tot. It reads a fourth tabulated
surface, `log10(n_e/n_tot)`, which is deliberately **not** part of `EOSTable::Eval()`:
only non-ideal MHD wants it, and charging every EOS call for another Hermite patch would
be a pure loss. An ideal gas returns 0, since it carries no composition — check
`IsGeneral()` first.

`<mhd>/ohmic_resistivity = eos` builds the diffusivity from it as

    eta = 230 sqrt(T)/x_e  +  5.2e11 lnL/T^1.5      [cm^2/s]

the electron-neutral (Blaes & Balbus) and Spitzer terms, added because the collision
frequencies add. The alternative `perna` carries only the first term and fits potassium
alone: it is good to a factor of ~3 over 1500–5000 K, but it has no saturation ceiling, so
it over-predicts x_e once the alkalis are fully ionized, and no hydrogen term, so it
under-predicts above ~6000 K, where the Spitzer term is also the larger of the two.
Selecting `eos` without a general EOS is a fatal error rather than a silent floor.

Cost is not a reason to avoid it: a table lookup is ~6× **cheaper** than evaluating the
`perna` fit, which is transcendental-bound (`log10`, `pow`, two `sqrt`). Solving Saha
per cell instead would be ~17× more expensive — which is the point of tabulating it.

---

## Restarts

The table is rebuilt from the input parameters when a run restarts. The restart file
carries the input, so a plain restart reproduces the same EOS — but overriding any
`eos_*` parameter on a restart command line **changes the EOS mid-run**, silently.

---

## Tests

| test | what it covers |
|---|---|
| `tst/scripts/hydro/hydro_general_eos.py` | `general_eos = gamma` must reproduce `eos = ideal` run for run, over plm/wenoz x llf/hllc x uniform/SMR; also asserts that `roe` aborts |
| `tst/scripts/mhd/mhd_general_eos.py` | the same for MHD, with hlld |
| `tst/scripts/hydro/hydro_general_eos_table.py` | `general_eos = table`: second-order convergence of a sound wave whose background sits in the hydrogen partial-ionization zone (`Gamma_1 ~ 1.30`), plus a contrast run in gamma mode so the test cannot pass on a silent fallback |
| `tst/scripts/mhd/mhd_general_eos_table.py` | the same for a fast magnetosonic wave |

`hlle` is exempt from the equality checks for the reason given above, and is only required
to stay within a factor of two of the ideal path.

---

## See also

[`docs/solar_convection.md`](solar_convection.md) — a worked example: what has to change
in a stratified problem generator when the ideal-gas assumption is dropped, and what the
halved pressure scale height does to the atmosphere above the photosphere.
