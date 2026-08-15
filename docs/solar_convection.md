# `solar_convection`: running it with a general EOS

`src/pgen/solar_convection.cpp` is a user problem generator (`-DPROBLEM=solar_convection`)
for a box of the solar surface layers: a convection zone below, a radiative atmosphere
above, coupled by a two-stream grey radiative transfer solver, with a CO5BOLD-style open
bottom that sets the inflow entropy. A reference input is
`inputs/hydro/solar_convection.athinput`.

This page covers what changes when it is run with [`eos = general`](general_eos.md), and
the two runtime parameters that were added for it.

---

## Why the pgen needed changing

The pgen was written for an ideal gas and assumed it in nine places. The worst was the
initial condition: the hydrostatic integration returned an EOS-consistent `(rho, p)`, but
the state written to memory used `e = p/(gamma - 1)`. Under the tabulated EOS the internal
energy also carries H2 dissociation and ionization latent heat, so `(gamma - 1) e`
overestimates the pressure by a factor of ~4 and the IC starts 75% out of hydrostatic
balance.

All nine sites now go through `pgen_eos_utils.hpp` and are EOS-agnostic. The diagnostic
worth remembering: **evaluate `(dp/dz + rho g)/(rho g)` on the `t = 0` snapshot before
trusting any new stratified IC.** It was 0.75 before the fix and 4.3e-3 after (the ideal
control gives 9.7e-4), and a large *uniform* residual points straight at a `p`-to-`e`
conversion rather than at the integration.

---

## `<problem>/T_top_fix`

```
<problem>
T_top_fix = 3500.0        # K, default 3500
```

Temperature of the fixed-`T` ghost lid at the top boundary. It was previously hard-coded.
3500 K was tuned for the ideal-gas run, whose atmosphere settles at 3525 K and therefore
sits consistently with it; the general-EOS atmosphere is far more rarefied and settles
colder, so the value is worth varying.

Doing so turns out **not** to matter for the atmospheric dynamics — see below — but the
parameter is useful for testing that, and the default reproduces the historical runs
exactly.

---

## `<problem>/sponge`: absorbing layer below the top boundary

```
<problem>
sponge      = true        # default false
sponge_zbot = 0.8         # bottom of the layer, as a fraction of the box height
sponge_c    = 0.1         # damping rate, in units of cs/dz
```

Damps the velocity toward zero over the top `(1 - sponge_zbot)` of the box, at rate
`sponge_c * cs/dz` — a fraction per cell-crossing time, the same resolution-aware idiom
the CO5BOLD bottom relaxation uses, so the layer behaves the same at every resolution. The
ramp is quadratic from the bottom of the layer to the top so that it opens smoothly and
does not itself reflect, and the update is implicit, `v /= (1 + rate*dt)`, hence
unconditionally stable.

**The kinetic energy removed is discarded, not thermalised.** The layer stands in for the
escape an unbounded atmosphere would provide; returning the energy as heat would go on
inflating the atmosphere, which is the problem being solved. Total energy is therefore
deliberately not conserved inside the layer and the sink is visible in the history file.
In practice it is small — in the runs below the total-energy drift *improved*, from +1.88%
without the sponge to +1.37% with it, because much of that "drift" was the accumulating
wave energy.

Off by default, so existing runs are bitwise unchanged.

### Why it is needed

Convection drives an acoustic flux of ~5-10e9 erg/cm^2/s up through the photosphere, and
it is the same under either EOS. What differs is the stratification it climbs through:

| | ideal gas (`mu = 0.602`) | tabulated EOS |
|---|---|---|
| `mu` at the photosphere | 0.602 | 1.26 |
| pressure scale height | 275 km | **126 km** (the solar value) |
| scale heights above the photosphere | ~5 | **~10** |
| wave amplitude gain, `rho^-1/2` | ~11x | **~134x** |

A wave amplitude grows as `rho^-1/2`, so the general-EOS atmosphere drives the same flux
to Mach 1 before it reaches the lid, and above `z/zmax ~ 0.57` the atmosphere is
shock-dominated rather than hydrostatic. That much is physical: real chromospheres shock.

What is not physical is that it never saturates. Measured on the dumps, radiative damping
cannot remove the energy (`t_rad/t_ac = 11-17`, so a wave crosses a scale height
adiabatically) and the top boundary reflects part of what arrives, so wave energy
accumulates indefinitely — 28% of a 192^2 box supersonic by `t = 20000`, still climbing.
Varying `T_top_fix` over 2500-4860 K changes nothing, which rules out the lid temperature
as the driver.

### What it does

96 x 64^2, tabulated EOS, identical except for the sponge:

| | no sponge | `zbot 0.8, c 0.1` | `zbot 0.7, c 0.3` |
|---|---|---|---|
| fraction of the box with Mach > 1 at `t = 20000` | 4.80% | ~0.1% | ~0.1% |
| wave energy above the photosphere, `t = 11000 -> 20000` | **x5.0** | x0.80 | x1.00 |
| `dE/dt` / acoustic luminosity, after `t = 12000` | +0.002..+0.009 | ~0 | ~0 |
| mid-atmosphere Mach, `t = 10000 -> 20000` | 0.109 -> **0.336** | 0.085 -> 0.066 | 0.070 -> 0.061 |

The gentler setting is enough, and is the recommended default: it damps only the top 20%
and leaves more of the atmosphere physical.

### Caveat

The sponge is a numerical device, not a physical model. It produces *a* steady state, but
that state is not proven to match what a genuinely taller box with a real escape boundary
would give — the sponged runs' mid-atmosphere is ~40% less dense than the unsponged
control's. A taller-box run remains the real validation.

---

## Offline analysis

The tabulated EOS cannot be called from host code or Python, and AthenaK has no
temperature output for non-relativistic hydro, so the diagnostics above are computed by
`tools/solar_convection/`, which dumps the code's own composition model to a
`(log rho, log T)` grid and interpolates it. See that directory's `README.md`.
