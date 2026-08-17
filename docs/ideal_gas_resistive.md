# Resistive MHD with an ideal gas and a tabulated electron fraction

Recommended settings for Ohmic-resistive hot-Jupiter runs that keep ideal-gas
thermodynamics. Reference input: [`inputs/mhd/deep_hot_jupiter_rt_ideal_xe.athinput`](../inputs/mhd/deep_hot_jupiter_rt_ideal_xe.athinput).

Every number below was measured on that problem — a 64 x 56 x 128 spherical-polar
ultra-hot Jupiter at `bbot = 10 G`, run for one rotation period (2pi/omega = 3.05e5 s) on
16 MPI ranks x 7 OpenMP threads — unless it says otherwise. For the general-EOS version of
the same problem see [`general_eos.md`](general_eos.md).

---

## The short version

```
<mhd>
eos               = ideal
ohmic_resistivity = eos
max_eta           = 1.0e13
use_rkg_sts       = false
dfloor            = 5.44e-12
pfloor            = 1.0e-1
tfloor            = 9.186e9      # = 200 K * Rgas; NOT tfloor_kelvin, see below
```

plus a `<units>` block and the `eos_*` composition keys, which the table build needs even
though the thermodynamics stay ideal.

---

## Why `ohmic_resistivity = eos` under an ideal gas

`eos` takes the electron fraction from the same Saha composition model the general EOS
uses — H, He, and the metal donors Na, K, Ca, Al, Mg, Fe, with condensation curves — and
forms

```
eta = 230 sqrt(T) / x_e          electron-neutral   (Blaes & Balbus 1994)
    + 5.2e11 * lnL / T^1.5       electron-ion       (Spitzer, lnL = 20)
```

The alternative, `perna`, is potassium-only sqrt-Saha at a fixed abundance: no saturation
ceiling (so it over-predicts x_e ~3x above 3500 K at low pressure), no hydrogen term (so it
under-predicts above ~6000 K), and no condensation or metallicity.

**It is free.** With `eos = ideal` the composition table is built for its x_e surface
alone: `tbl.active` stays false, every thermodynamic accessor keeps the ideal-gas
expression, and only `ElectronFractionKelvin()` reads the table. One bicubic Hermite patch
costs no more than `perna`'s `log10 + pow(10,x) + 2 sqrt`. Measured at **-0.02 ms/cycle**,
i.e. 7% faster than `perna` as well as more accurate. Startup pays a one-off ~5 s table
build.

`eos_metal_ionization = true` is **not optional**. Below ~4000 K hydrogen's 13.6 eV
contributes essentially no free electrons and the metals supply all of them; without them
x_e collapses and the configuration is pointless.

### What it does not buy

An ideal gas carries a fixed mu. `Rgas = 4.593e7` implies mu = 1.810, while the
composition model gives ~1.26 in the same gas, so an ideal run's temperature is ~1.44x the
general EOS's there. With `d ln x_e / d ln T ~ 13` in the alkali range, that is ~2 orders
of magnitude in x_e — far more than the factor ~3 that separates `perna` from full Saha.

So this gives the best electron fraction **consistent with the ideal run's own (rho, T)**.
It does not make an ideal run resemble a general-EOS one. Use `eos = general` when the
magnetic coupling itself is the science; the general EOS costs **4.3x per cycle / 3.2x per
simulated second** (measured, 200 cycles, constant eta, against the identical ideal setup).

---

## `max_eta = 1e13` and `use_rkg_sts = false`

These two go together and were chosen by measurement, not by argument.

`max_eta` is applied twice, and the first application is the important one: x_e is floored
at `230 sqrt(T)/max_eta` before the Spitzer term is added, so **the cap is really an
electron-fraction floor**.

A cap is unavoidable. Over a real snapshot, the uncapped Saha eta has a median of 1.4e9
cm^2/s but a 99th percentile of 3.2e28 and a maximum of 5e59, because eta ~ exp(+chi/2kT)
once the alkalis recombine. What each choice truncates:

| `max_eta` | % cells capped | % mass in them | % magnetic energy in them |
|---|---|---|---|
| 1e12 | 9.5 | 9.9e-3 | 2.68 |
| **1e13** | 6.1 | 4.9e-3 | **0.87** |
| 1e14 | 3.3 | 1.4e-3 | 0.26 |
| 1e15 | 2.3 | 1.5e-5 | 0.20 |

(Computed without condensation, which would remove donors and raise eta further, so these
are lower bounds on the capped fraction.)

### The cost, measured

Wall seconds per simulated second, matched time windows, same problem:

| window | 1e14 + STS | 1e14, no STS | **1e13, no STS** | mean dt |
|---|---|---|---|---|
| 6e4 - 9e4 | 0.01720 | 0.01929 | **0.01064** | 1.66 |
| 1.2e5 - 1.5e5 | 0.02043 | 0.01981 | **0.01374** | 1.39 |
| 1.8e5 - 2.1e5 | 0.02414 | 0.01994 | **0.01647** | 1.18 |
| 2.4e5 - 2.6e5 | 0.02723 | 0.02030 | **0.01823** | 1.05 |

`1e13 + no STS` is **1.49x faster** than the `1e14 + STS` configuration at late time, at
essentially the same timestep (1.036 vs 1.045 s). At 1e13 the resistive timestep lifts to
dt_diff ~ 9.5 s, far above the hydrodynamic CFL value, so diffusion never binds and
super-stepping has nothing to do.

Whole-rotation totals for the same four runs (16 x 7, to t = 3.05e5):

| configuration | wall [s] | cycles | vs baseline |
|---|---|---|---|
| `max_eta = 1e14`, STS on (baseline) | 6347.7 | 221210 | 1.000 |
| `max_eta = 1e14`, STS on, `pfloor` 1e-3x | 6257.5 | 219351 | 0.986 |
| `max_eta = 1e14`, STS off | 5877.8 | 316614 | 0.926 |
| **`max_eta = 1e13`, STS off** | **4129.4** | 219988 | **0.651** |

All four completed the rotation with finite mass and total energy agreeing to 5-6 digits.
Note the STS-off run takes 43% more cycles yet still finishes 7% sooner.

### The accuracy cost, measured

1e13 against 1e14 at matched simulated time (t = 1.22e5):

| radial shell | rho relative difference | shell magnetic energy ratio |
|---|---|---|
| 0 - 32 (r/R_p <= 1.17) | 2.4e-7 .. 2.0e-5 | 1.0000 - 1.0004 |
| 40 | 3.2e-3 | 1.011 |
| 48 | 2.0e-2 | 1.038 |
| 56 - 63 | 1 - 2e-2 | 0.954 - 0.960 |

The deep atmosphere is unchanged to 4-5 digits. Only the outer eight radial shells move,
by ~4%, and those are the shells where the floors are also active.

### When super-stepping IS worth it

RKG stage count is `s = ceil(-a + sqrt((a+1)^2 + tau(3+2a)))` with a = 0.5 and
tau = dt_hydro/dt_diff, **floored at s = 3**. Break-even is tau ~ 1.5-2, and one sub-stage
costs ~0.14 of a full RK2 step, so `cost_on/cost_off = (1 + 0.14 s)/max(1, tau)`.

On this problem with `max_eta = 1e14` the sign flipped **during the run** as dt_hydro
decayed: STS was 1.12x ahead at t = 7.5e4 (dt = 1.66, tau ~ 1.8), level by t = 1.3e5
(dt = 1.39), and 1.34x behind by t = 2.5e5 (dt = 1.05, tau ~ 1.1). So "turn STS on at high
max_eta" is not a standing rule — **watch dt against dt_diff**, both of which are printed
every `ndiag` cycles. Turn it on only when dt sits well below the CFL value.

Note also that the super-step runs *after* the complete RK2 step (`driver.cpp`), i.e. Lie
splitting, where with STS off the resistive fluxes sit inside the RK2 stages.

---

## Floors

```
dfloor = 5.44e-12
pfloor = 1.0e-1
tfloor = 9.186e9        # 200 K
```

**State `tfloor` in code units, not `tfloor_kelvin`.** `tfloor_kelvin` converts with
`<units>/mu`, which is 1 in this input, while the run's `Rgas = 4.593e7` implies mu = 1.797
— the floor would land 1.8x off. Under `eos = general` the same key IS correct, because
there the scale has mu divided back out. The code-unit value is `T_K * Rgas`.

These floors are **not** inert. Late-run rates, as a fraction of C2P calls: dfloor 13.9%,
efloor 8.0%, tfloor 0.44%. But they are confined to the outer ~12% of the domain in radius
— no cell below shell 48 (r/R_p = 1.25) touches any of them, and below shell 48 there are
three to nine orders of magnitude of headroom.

### dfloor can be relaxed 1000x; it costs 5.9x

At 10 G, dfloor is what bounds the Alfven speed: v_A = B/sqrt(4 pi rho) is 1.2e6 cm/s at
5.44e-12 against a sound speed of only 4.1e5, so the floor already sits where v_A ~ 2 c_s.

| dfloor | outcome | wall / simulated second |
|---|---|---|
| 5.44e-12 | healthy | 0.0208 |
| 5.44e-15 | healthy, physics unchanged | 0.1225 (**5.9x**) |
| 5.44e-18 | **NaN at t ~ 7.4e4** | — |

At 5.44e-15 the deep atmosphere is unchanged (density agrees to 1e-5 .. 3e-4, shell
magnetic energy to 0.3% below r/R_p = 1.26; mass and total energy to 6 digits) and only the
outer eight shells move, by ~6%. Use it if upper-atmosphere field structure is the target;
otherwise it buys a 6% correction for a 6x cost.

5.44e-18 is fatal, and not gently: by then the floor is *below* the density the gas
actually reaches (rho_min ~ 2e-14), so it fires only 3.6e-6 of the time — it has stopped
being a bulk floor and become a rare-event catcher. Those rare events are the whole point,
because one cell can then reach v_A ~ 1e9 cm/s within a single step. Note that equal mean
dt does **not** imply equal robustness: the 1e-3x and 1e-6x runs had the same timestep for
the first 20 minutes and only one of them survived.

### Do NOT relax pfloor or tfloor

Cooler gas is more resistive here — `d ln eta / d ln T` reaches -13 around 1800 K and -26
at 1000 K — so a lower temperature floor drives eta to `max_eta` over a widening region,
dt_diff collapses, and the run either stalls (measured: 3053 s of simulated time in 7
hours, super-stepping stage count exploding to 258) or goes NaN outright (measured: t =
1.2e4, reached in 95 s of wall time). Both failure modes were reproduced at 1e-3x these
values, with dfloor relaxed and with dfloor held at baseline.

`pfloor` is in any case nearly redundant: `dfloor * Rgas * tfloor = 0.050` already bounds p
from below, and pfloor is 0.1 — it has teeth only in that factor-2 window, which is why the
energy floor fires at 8% and the temperature floor at 0.44%.

Relaxing `pfloor` **alone** by 1e-3x confirms this directly: the run completes the rotation
in 0.986x the wall time with the same mass and total energy, and the flooring simply moves
from one floor to the other — energy-floor rate 7.97% -> 2.12%, temperature-floor rate
0.44% -> 6.95%. So it is harmless, and pointless. `tfloor` is the one that must not move;
it is also the one the resistivity couples to, through `d ln x_e / d ln T`.

---

## Where this model stops being physical

**It is Ohmic only.** `resistivity.hpp` says so — ambipolar diffusion and the Hall effect
are not implemented. The ratios `eta_Hall/eta_Ohm` and `eta_AD/eta_Ohm` both have n_e cancel
exactly, so they depend only on rho, T and B, not on the composition model. At 10 G:

| rho | eta_Hall/eta_Ohm | eta_AD/eta_Ohm | dominant |
|---|---|---|---|
| 1e-5 | 1.7e-3 | 1.5e-8 | Ohmic |
| 1e-7 | 1.7e-1 | 1.5e-4 | Ohmic |
| 1e-8 | 1.7 | 1.5e-2 | Hall |
| 5.4e-12 | 3.1e3 | 5.2e4 | ambipolar |

Hall overtakes Ohmic at rho ~ 1.7e-8, ambipolar at ~1.2e-9. On this problem that puts the
Ohmic-only approximation in the clear through shell 49 (r/R_p <= 1.26), holding **99.999%
of the mass and 99.1% of the magnetic energy** — and wrong in the outer 14 shells, which
are also where the floors and the cap are active. That region is numerically controlled
rather than physically resolved and should not be interpreted either way.

**LTE is assumed** — no stellar XUV photoionization, which in a real ultra-hot Jupiter
upper atmosphere maintains x_e well above LTE Saha. Same region as the Hall/ambipolar
failure, and both errors point the same way.

**`max_eta` is not physics.** In the capped cells Saha wants eta ~ 1e29, i.e. an insulator
whose field should relax freely; capping makes it orders of magnitude more conducting, so
the field stays frozen where it should not.

---

## The temperature sensitivity, and why it drives everything above

For trace Saha, `d ln x_e / d ln T = 3/4 + chi/(2kT)` — the `25188/T` in `perna` is exactly
chi_K/2k for potassium's 4.34 eV. Evaluated at each shell's own (rho, T):

| r/R_p | median T | eta [cm^2/s] | d ln eta / d ln T | |
|---|---|---|---|---|
| 1.003 | 12177 K | 9.3e6 | -2.4 | weak (Spitzer) |
| 1.085 | 7576 K | 6.4e7 | -8.0 | moderate |
| 1.167 | 3794 K | 4.2e9 | -5.7 | moderate |
| 1.208 | 2448 K | 5.0e10 | -9.6 | strong |
| 1.249 | 1765 K | 4.0e11 | **-13.2** | strong |
| 1.326 | 3547 K | 2.6e8 | -1.3 | weak (Spitzer) |

Three regimes: above ~8000 K hydrogen is ionized and the Spitzer term carries eta, so the
slope is just -1.5; between 3000 and 8000 K metals saturate while hydrogen switches on;
below ~2500 K alkali ionization gives -10 to -13, steepening to -26 at 1000 K. There is no
flat plateau in between — hydrogen's exponential tail switches on before the metals can
flatten the curve.

The resistivity is therefore a strong function of T only in the cool band around the
temperature minimum, here r/R_p ~ 1.21-1.29. That band is exactly what `tfloor` and
`max_eta` act on, which is why both have so much leverage on cost and stability.

---

## Diagnostics

Add the event-counter output:

```
<output2>
file_type = log
dt        = 3.05e3
```

It writes `eos_dfloor eos_efloor eos_tfloor eos_vceil eos_fail c2p_it fofc`. Two things to
know: the counters **reset after every write**, so each row is per-interval and not
cumulative; and a row is written only when some counter is non-zero, so with a large `dt`
the only row you get is the last one.

`ndiag` also prints `dt` and, with super-stepping on, `dt_diff` and the stage count — which
is what tells you whether `max_eta` is binding. A `dt_diff` that stays constant to within a
few percent while the flow evolves means eta is pinned at the cap.
