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
| `eos_a_metal` | `16.0` | atomic weight of the inert "metal" component (Z = 1 − X − Y). It contributes particles but neither dissociates nor ionizes. |
| `eos_h2` | `true` | include H2 and its dissociation equilibrium, and its rotational/vibrational energy |
| `eos_ionization` | `true` | include Saha ionization of H, He and He+ |
| `eos_radiation` | `false` | add `a T^4` energy and `a T^4/3` pressure on top of the tabulated gas |

`eos_xh` and `eos_yhe` must be non-negative and sum to at most 1; the run aborts
otherwise. Radiation is added **analytically on top of** the tabulated gas rather than
being tabulated, which makes the switch exact and keeps the tabulated surfaces smooth.

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

The composition model covers H2, H, H+, He, He+, He++, free electrons and an inert metal
component, closed by charge neutrality.

**Not included:** Coulomb (non-ideal) corrections, electron degeneracy and pressure
ionization, excited bound states beyond ground-state statistical weights, and any
ionization of the metal component.

It is therefore valid for the **weakly coupled, non-degenerate** regime — gas-giant and
stellar envelopes — and not for deep stellar interiors or degenerate objects.

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
