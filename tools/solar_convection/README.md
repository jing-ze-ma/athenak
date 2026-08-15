# solar_convection: offline analysis for the ideal-vs-general-EOS comparison

Analysis for paired `solar_convection` runs that differ only in their equation of state —
an ideal-gas control against `eos = general`, `general_eos = table` (H2 dissociation +
H/He Saha, radiation off). Everything here reads AthenaK `.bin` dumps directly and
re-evaluates each run's own EOS in Python, because:

* **AthenaK cannot output temperature for Newtonian hydro** — there is no `hydro_t` in
  `var_choice` (`src/outputs/outputs.hpp`) — and
* the tabulated EOS lives in a `DvceArray`, so it cannot be called from host code or from
  Python at all.

`eos_dump.cpp` closes that gap: a standalone program that `#include`s the code's *own*
`src/eos/eos_composition.hpp` and dumps a (log10 rho, log10 T) grid, so the offline physics
cannot drift from what the solver did. No Kokkos needed:

```bash
g++ -O2 -I<athenak>/src -o eos_dump eos_dump.cpp
./eos_dump eos_grid.bin              # writes the grid eoslib.py interpolates
```

## Layout and configuration

The scripts expect a run tree

```
$SOLEOS_ROOT/            default /orion/u/jinma/ATHENAK/soleos
  <run>/bin/sun.hydro_w.*.bin
  <run>/sun.hydro.hst
  tools/eos_grid.bin     the table dumped by eos_dump
  plots/                 everything is written here
```

Set `SOLEOS_ROOT` to point them somewhere else. `bin_convert` is imported from the
repository's own `vis/python/`, located relative to this directory, so a checkout
anywhere works. `slices.py` also picks up `run/mesa2.mplstyle` if it is present
(override with `MPLSTYLE`).

## The tools

| script | what it does |
|---|---|
| `dumps.py` | Finds a run's dumps and rejects truncated ones. A GPFS quota can silently truncate a dump mid-write — the run continues, the file is short, nothing is logged. Any dump more than 64 KB below the run's largest is dropped, with the slot kept as `None` so frame numbers keep matching file numbers. |
| `eoslib.py` | `IdealEOS` and `TableEOS` behind one interface. **Every accessor takes `(rho, e_specific)`** — not `(rho, T)`, not energy density. `TableEOS.temperature` inverts by vectorised bisection on log10 T. Derivatives (chi_rho, chi_T, cv, Gamma_1, grad_ad) come from `np.gradient` on the log-log grid. `get_kapr` mirrors the pgen's opacity so tau matches. |
| `compare.py` | `compare.py [idx] [ideal_run] [table_run]` — the headline comparison: horizontally averaged structure, photospheric granulation, conservation history, and a text summary. Outputs are tagged `_<ideal>_<table>` unless the pair is the default `ideal`/`table`. |
| `slices.py` | `slices.py <run> [first] [last]` — vertical x3-x1 slices of vz, T, log10 dtau and Mach through mid-x2. |
| `tau_surface.py` | `tau_surface.py <run> [first] [last]` — granulation on the *corrugated* tau=2/3 surface rather than a flat geometric height, linearly interpolated in tau. `TAG` appends to the output directory so one run can be rendered at two sets of colour limits. |
| `lidtest.py` | Sensitivity of the atmosphere to `problem/T_top_fix`, across runs that differ only in the lid temperature. |
| `acoustic.py` | Acoustic energy budget: up/down-going characteristic fluxes, the reflection coefficient `R = F-/F+`, wave energy stored above the photosphere, and the radiative relaxation time against the acoustic crossing time. Distinguishes a trapping boundary from an undamped atmosphere. |

Colour limits in `slices.py` and `tau_surface.py` default to the values the older
`run/plotsun.py` and `run/plot_tau_surface.py` used, so figures stay comparable with
earlier work; all are overridable by environment variable (`VZLIM`, `TMIN`, `TMAX`,
`VMIN`, `VMAX`, `MACHMAX`). Note that the `vz` defaults **saturate** for the tabulated
runs, whose atmospheric velocities reach tens of km/s.

Two correctness notes carried over from those originals: they computed T with the
ideal-gas formula (wrong for `general_eos = table`) and `plotsun.py` used a stale opacity
(C=1e2, no floor, no Kramers). Both are fixed here.

## Caveats

* `.bin` output is **single precision**, so agreement beyond ~1e-6 is not meaningful.
* `compare.py`'s `_tau3` and the tau surfaces integrate the pgen's grey opacity, not the
  two-stream solver's actual optical depth; they agree in the mean but not cell by cell.
* `acoustic.py`'s characteristic decomposition is linear, so it is quantitative only while
  the atmosphere is subsonic.
