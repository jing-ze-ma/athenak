---
name: general-eos-stage4-rho-e-table
description: "SCOPE (not started). Baseline moved 2026-08-17: 3.07x ideal now, not 4.5x, and the no-inversion proxy re-measured at 2.00x. So Stage 4 is worth ~1.5x, not ~2.1x. Design, call-site audit, risks and verification plan below still stand."
metadata:
  type: project
---

Scoped 2026-08-16. Motivation and the measured ceiling are in [[general-eos-table-cost]]:
49% of the time loop is libm, only 14% is the fluid solver, and collapsing the root find
is worth a MEASURED 43% (4.52x -> 2.60x ideal, with a real table expected near 2.1x).

## RE-BASELINED 2026-08-17 -- the prize is smaller than scoped

Two things moved. The BASELINE fell to **3.07x ideal per cycle** (2.32x per simulated
second) after `2a154e7d` + `d1c289dd`, and the no-inversion ceiling was re-measured on the
current binary:

| configuration (500 cycles, 8x4, constant eta) | zone-cycles/cpu-s | x ideal |
|---|---|---|
| ideal | 9.319e6 | 1.00 |
| general, root find forced to ONE Newton step | 4.654e6 | **2.00** |
| general, current at the time | 2.128e6 | 4.28 |

The proxy (`logtol = 1e30`, a compile-time constant, so it needs a rebuild) is a genuine
lower bound and does not depend on any model of where the time goes. It says killing the
inversion was worth 53% THEN. Against today's 3.07x baseline the same endpoint implies
roughly **3.07x -> ~2.0x, i.e. ~1.5x, not the ~2.1x-from-4.5x originally scoped.**

**Confidence, stated honestly when the user asked (2026-08-17):** ~85% that it delivers at
least 40% off, because the ceiling is MEASURED rather than modelled and has now been
measured twice independently (43% and 53%); ~65% that it reaches ~2.0x ideal, since the
real table pays 2 log10 + 3-4 patches against the proxy's one residual plus final Eval --
close enough that landing NEAR the proxy, not below it, is the expectation. The reason to
trust this more than the estimate that missed (patch-counting predicted 2.5x, delivered
10%) is that the failed one was a MODEL of where time went, demolished by the profile
showing HermitePatch at 6%; this one measures the endpoint directly.

**Unverified risk:** the non-rectangular (rho,e) domain, still the main correctness hazard.
An attempt to check it with a throwaway Saha model failed -- the model was visibly wrong at
the cold end (full dissociation at 31 K) -- so it needs a host-side harness calling the real
`EOSCompositionModel`. That remains task 1.

**Do the stall check first: DONE.** It was not the EOS at all, see
[[eventlog-mpi-deadlock]]. Benchmarking Stage 4 on this problem was impossible until that
was fixed, and it would have credited the new table with a speedup that was really a
deadlock going away.

## The design: TWO tables, not a replacement

Keep the existing `(log rho, log T)` table exactly as it is, and ADD a second one on
`(log rho, log e)`. This is the whole trick that keeps the change small:

- **Hot path** reads the new (rho,e) table. `Temperature`, `Pressure`, `Gamma1`, `mu`,
  `x_e`, `c_v`, `chi_rho`, `chi_T` all become DIRECT LOOKUPS -- no iteration.
- **Everything else is untouched.** The floors need `e(rho,p)` and `e(rho,T)`, and the pgens
  need `DensityFromPressureTemperature`, `GradAd`, `EintFromP`; all of those are the (rho,T)
  direction and keep working against the old table with no edits at all.
- The new table is BUILT by inverting the old one once at startup, so the two cannot drift.
  ~127k root finds at build; the existing build already does ~1.1M Saha solves in ~5 s.

Replacing rather than adding was considered and rejected: the pgens call the (rho,T)
accessors from inside setup `par_for`s, so they need that direction ON DEVICE, and giving
it to them from a (rho,e) table means a 2D inversion.

## Call-site audit (done, `grep` over src/)

**Hot path — becomes a lookup:** `eos/general_c2p_{hyd,mhd}.hpp` (Temperature);
`diffusion/resistivity.cpp` (MeanMolecularWeight, ElectronFraction, per cell per stage);
`diffusion/conduction.cpp` (SpecificHeatCv in NewTimeStep); `srcterms/srcterms{,_newdt}.cpp`
(Temperature, MeanMolecularWeight).

**`utils/wb_background.hpp` is the awkward one** -- 25 EOS calls: 10 Pressure, 7
Temperature, 2 SpecificHeatCv, 2 ChiT, 1 each Enthalpy/ChiRho/EnergyFromTemperature/
EnergyFromPressure. Per cell per stage when well-balancing is on (solar_convection), and per
GHOST cell per stage from the dhj outer-x1 BC even when it is off. Most of those map onto
the new table directly; the two `EnergyFrom*` stay on the old one.

**Setup only, no change:** `pgen/*` via `pgen_eos_utils.hpp`
(DensityFromPressureTemperature x2, EnergyFromPressure x2, EnergyFromTemperature x2, ...),
`{hydro,mhd}_wellbalance.cpp` (`SetWbBackgroundPressure`, runs once at init and
early-returns unless `use_wellbalance_static`).

## Surfaces and memory

New table stores, at each (log rho, log e) node: log10 T, log10 p, Gamma_1, mu, log10 x_e,
c_v, chi_rho, chi_T = 8 surfaces x 4 (value + 2 first derivatives + cross) = 32 doubles.
At the dhj grid (281 x 451) that is ~32 MB against the current 16 MB. Fine on CPU and on a
16-64 GB GPU. `chi_rho`/`chi_T` (wb_background only) and `c_v` (conduction dt only) could be
made optional if memory ever matters.

## The one real subtlety: the domain is not a rectangle

A rectangle in (log rho, log T) does NOT map to a rectangle in (log rho, log e): at fixed
rho, e runs from e(T_min) to e(T_max), and both bounds move with rho. So the corners of the
new grid fall outside the sampled region and must be filled by continuation. `Interpolate`
already does linear continuation outside the table, but the BUILD has to decide what to put
at those nodes. Get this wrong and the error shows up only in rare cells.

First task, before any code: compute the actual e range over the intended (rho,T) box and
look at how badly non-rectangular it is.

## Verification plan

**It will NOT be bitwise** -- every value moves by interpolation error. So:
1. `tst/scripts/{hydro,mhd}/*_general_eos_table.py` (convergence) and
   `mhd/mhd_eos_electrons.py` (physics + x_e) must pass unchanged. Note the *equality*
   tests in `*_general_eos.py` use `general_eos = gamma` and are unaffected.
2. New round-trip check at build: `T(rho, e(rho,T)) == T` to interpolation accuracy over the
   whole grid, and the same for p.
3. Short dhj and solar_convection runs against the (rho,T) build: agreement should be at
   interpolation level, NOT round-off. Quantify it rather than eyeballing.
4. Re-run the cost benchmark; target ~2.1x ideal against 4.43x today.

## Effort

~2-4 focused days. `eos_table.cpp` build path ~150 lines; `eos_table.hpp` accessors ~100;
`eos.hpp` routing ~50 and the part needing most care (which accessor uses which table);
`general_c2p_*.hpp` get SIMPLER; wb_background is an audit rather than a rewrite.
Verification and benchmarking is a day on its own.

## Do the cheap wins first or not at all

If Stage 4 happens, `2dedcbdb`'s patch saving survives but the transcendental-elimination
idea in [[general-eos-table-cost]] becomes MOOT -- it optimises a root find that would no
longer exist. Do not spend effort there first.
