---
name: resistivity-xe-table-for-ideal
description: "DONE (5d1c3436): ohmic_resistivity=eos now works under eos=ideal via a table built for x_e alone. Measured 7% FASTER than perna. Two latent bugs found on the way."
metadata:
  type: project
---

Asked 2026-08-16: with `eos = ideal` (to avoid the general EOS's ~4.4x cost), can the
resistivity still beat the potassium-only `perna` fit? **Yes, and at no runtime cost.**

## Why it is free — MEASURED, not argued

From the cost decomposition in [[general-eos-table-cost]] (dhj, 500 cycles, 16x7):
`general(table) + perna` = 88.17 ms/cycle, `general(table) + eos` = 88.15. The x_e Hermite
lookup is **-0.02 ms/cycle against perna**, i.e. free. A bicubic patch plus two log10 costs
no more than perna's `log10 + pow(10,x) + 2 sqrt`. This confirms the older standalone
benchmark in [[resistivity-perna-uhj]] (bilinear lookup 13.5 ns vs perna 82 ns).

The 61 ms/cycle the general EOS costs is the `T(rho,e)` root find in ConsToPrim, which an
IDEAL gas never performs. So an ideal run can have the table's x_e without paying for the
table's thermodynamics.

## Accuracy gained over perna

perna is K-only sqrt-Saha with a fixed a_K = 1e-7: no saturation ceiling, no hydrogen, no
Na/Ca/Al/Mg/Fe individually, no condensation, no metallicity. Ratio perna/full-Saha
(from [[resistivity-perna-uhj]]):

| T | ~0.01 bar | ~1-5 bar |
|---|---|---|
| 1500-3000 K | 0.7-2.1 | 0.7-0.8 |
| 3500-5000 K | 2.8-3.3 | 1.0-1.2 |
| 7000-8000 K | 0.7-0.3 | 0.5-0.3 |

So ~factor 3 over most of the range, plus the qualitatively missing behaviour: saturation
once the alkalis are fully ionized, hydrogen taking over above ~6000 K, condensation, [M/H],
and the Spitzer electron-ion term that `ResistivityEOS` already adds.

## THE CATCH, which must be stated whenever this is offered

An ideal gas has a FIXED mu, so its temperature is `T = p/(Rgas rho)` with that mu baked in.
For the dhj setup `Rgas = 4.593e7` means mu = 1.810, while the table reports mu = 1.258 at
the grid centre -- the ideal run's T is ~1.44x the general EOS's there. Since
`d ln x_e/d ln T ~ 13`, that is a factor ~100 in x_e.

So: this gives the BEST x_e consistent with the ideal run's own (rho,T), which is a real and
free improvement over perna. It does NOT make an ideal run's magnetic coupling resemble a
general-EOS run's -- the temperature mismatch across H2 dissociation dominates and no
resistivity model can repair it. Do not oversell it as "general-EOS accuracy for free".

## Scope (contained, ~half a day + a test)

`BuildGeneralEOS` is currently called ONLY from the `General{Hydro,MHD}` constructors, and
`IdealMHD`'s ctor sets neither the cgs scales nor a table. Needed:

1. `IdealMHD`/`IdealHydro` ctor: when asked, set `dens_cgs/pres_cgs/temp_cgs` from `punit`
   and build the table. Building the full 16-surface table and using only ITXE is wasteful
   (~16 MB, ~5 s at startup) but needs no new build code; an x_e-only variant is tidier.
2. Relax the fatal guard in `Resistivity::SetResistivity` from `IsGeneral()` to
   `tbl.active`.
3. Add a Kelvin-taking entry point, e.g. `ElectronFractionKelvin(rho_cgs, T_K)` ->
   `tbl.ElectronFractionCgs`. **Do not route the ideal path through
   `ElectronFraction(d,e,t)`**: it multiplies by `eos_data.temp_cgs`, which the ideal EOS
   never sets (stays 1.0). This is the same mu trap as `tfloor_kelvin`; the resistivity
   already has T in Kelvin for both branches, so pass that.
4. Requires `<units>` and the `eos_*` composition keys in an otherwise ideal-gas input.
5. Regression test: extend `tst/scripts/mhd/mhd_eos_electrons.py`, which already has the
   machinery for banner-vs-fluid x_e comparisons.


## DONE 2026-08-16, commit `5d1c3436`

`IdealMHD` builds the composition table when, and only when, `ohmic_resistivity = eos` is
selected, and **leaves `tbl.active = FALSE`** -- `active` is what every thermodynamic
accessor in `EOS_Data` tests to choose table-vs-ideal, and such a run wants the ideal gas
for all of them. Only x_e comes from the table, via a new `ElectronFractionKelvin()` gated
on the new `EOS_Data::xe_from_table`. That entry point takes KELVIN deliberately: the
`(d,e,T)` form scales by `EOS_Data::temp_cgs`, which an ideal gas never sets. The
resistivity holds kelvin on both branches so both now use it; the general path is **bitwise
unchanged, verified over 22 dumps**.

**MEASURED 20.35 ms/cycle against perna's 21.94 -- 7% FASTER as well as more accurate.**
Startup pays a one-off ~4.7 s table build.

**Regression test:** `mhd_eos_electrons.py`'s old "eos = ideal must abort" assertion was
exactly the behaviour being added, so it is replaced by a stronger positive check -- with
`<problem>/Rgas = 3.5765e7` (mu = 2.3247) the ideal gas sits at the same 1985 K the table
gives that background, so an ideal run reading the same x_e must damp the Alfven wave
identically. **It agrees to 0.03%.** `eos = isothermal` is asserted to abort instead.

## Two latent bugs found on the way, both fixed in the same commit

1. **`HotJupiterParam` was uninitialised** and only filled when `<problem>/hot_jupiter` is
   true, yet `SetResistivity` divides by its `Rgas` on the ideal branch unconditionally.
   Every non-hot-Jupiter problem using `perna` was dividing by stack garbage. Now
   zero-initialised, and the resistivity refuses a non-positive Rgas.
2. **The probe for "does the resistivity want x_e" must not be `GetOrAddString`** -- that
   ADDS the parameter, and `Resistivity`'s constructor keys off it merely EXISTING, so a run
   with no resistivity acquired one of type `'none'` and aborted. Same class as the
   `GetOrAddReal`/`tfloor` trap in [[general-eos-project]]. Caught by `mhd_general_eos`.

**Harness note:** the fake test tree at `xetest/tstsim` APPENDS to `*-errs.dat` across
runs, so `athena_read.error_dat` eventually fails to reshape. `rm -f build/src/*-errs.dat
build/src/*.hst` between modules. `run_tests.py` does not have this problem -- it recreates
`build/`.


## 2026-08-17: an ideal run at max_eta = 1e14 falls into the DIFFUSIVE limit

Measured while comparing ideal vs general with `ohmic_resistivity = eos` on dhj at
bbot = 10 G. The ideal run's dt collapses partway through to exactly **1.1731 s**, which is
dt_diff at max_eta = 1e14 with STS off -- it becomes resistivity-limited. The general run
alongside it stays hydro-limited at 15.23 s the whole time.

The cause is the same fixed-mu problem this note is about, showing up as a TIMESTEP
pathology rather than an x_e error: the ideal run's log shows dfloor firing ~1e6 times
while the general run's event log is **completely empty**. mu = 1.810 drives it into cold
floored states the general EOS never reaches.

So an ideal run is only ~2.3x cheaper per simulated second than the general one WHILE it
stays hydro-limited, and once it drops to dt_diff that advantage largely evaporates. Use
`max_eta = 1e13`, which is what `docs/ideal_gas_resistive.md` recommends anyway -- it lifts
dt_diff to ~9.5 s so diffusion never binds. See [[general-eos-table-cost]].
