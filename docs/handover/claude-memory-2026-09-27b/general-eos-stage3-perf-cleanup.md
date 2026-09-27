---
name: general-eos-stage3-perf-cleanup
description: "Structural cleanup so the general EOS solves T once per cell instead of ~11-30 times — DONE, verified and committed as 7d1fed3f"
metadata:
  node_type: memory
  type: project
  originSessionId: f8810972-05f3-4467-b118-5f661b8a68a2
  modified: 2026-08-12T03:54:56.994Z
---

Done 2026-08-11, verified and committed 2026-08-12 as **`7d1fed3f`** on `general-eos`
(see [[general-eos-project]]). Before writing the analytic EOS (Stage 3), many call sites
asked the EOS for several quantities on the SAME `(d,e)`, each independently root-finding
`T` — ~11 root finds per cell per stage without well balancing, ~30 with. The cleanup makes
it one.

## What was changed (~630 lines across 23 files)

1. **`src/eos/eos.hpp`** — every accessor now has a primary `(d,e,T)` form plus a `(d,e)`
   convenience wrapper that solves `T` itself. `Temperature(d,e,tguess=-1.0)` is the single
   expensive entry point; `tguess<=0` means "no guess, bracket yourself". Hot-path kernels
   MUST use the `(d,e,T)` form. Added `BelowPressureFloor(d,e,T)`, `EnergyFloorBound(d)`
   (a cheap UPPER bound on `e(d,pfloor)`, exact for an ideal gas) and `ApplyEnergyFloor`.
2. **`general_c2p_{hyd,mhd}.hpp`** — ONE `T` solve per cell (was 3). Floors are tested on
   `p` and `T` instead of `e`, so inversions run only where a floor actually trips. `wtemp`
   is the solved value passed out, fed back as the next stage's warm start.
3. **`reconstruct/{wenoz,ppm}.hpp`** — six duplicated floor blocks collapsed into
   `eos.ApplyEnergyFloor(d, e)`, exact inversion behind the cheap bound.
4. **`utils/wb_background.hpp`** — one `T` per state. `WBEnergyFromEnthalpy` went from 4
   root finds per Newton iteration to 1.
5. **Static WB background precomputed** — `pwb`/`pfacewb` on Hydro and MHD, filled once by
   `SetWbBackgroundPressure()` from `Driver::InitBoundaryValuesAndPrimitives`. Biggest win
   (~19 calls/cell/stage on a background that never changes). Consumers:
   `WbStaticPiecewiseLinearDerX{1,2,3}` (which lost their now unused `eos` parameter),
   `RemoveWbFlux`, the spherical-polar coordinate source terms.
6. **Live-state pressure from the cache** — `coordinates.cpp` now reads `wder(IDPR)`.
   `SrcTermsGnomonicEquiangle` serves BOTH Hydro and MHD, so it takes `wder` as an argument.
7. **`diffusion/resistivity.cpp`** — reads `pmhd->wtemp`.

## Verification (complete)

`/orion/u/jinma/ATHENAK/eos test/eos_compare.py` runs each case twice from ONE input with only
`<hydro|mhd>/eos` overridden. (Path is `/orion/u/jinma/ATHENAK/eostest/`.)
- 16 built-in cases (lw/lwm/sod/ot/otf/smr/smrm/amr/cond/cool/coolm): worst-difference
  numbers **byte-identical** to the stashed baseline (`eostest/eos_compare_baseline.txt`).
- `wbs_plain`, `wbs_perturb`, `wb_iso`, `wb_isen` (`PROBLEM=hse_atm`): `wb_iso`/`wb_isen`
  bitwise identical baseline-vs-after; the two static-WB cases differ only in the velocity
  residual at ~1e-15 absolute (max |vely| ~1e-4), for BOTH the ideal and general paths —
  item 5 changes the ideal path too. No degradation: the residual is the same order before
  and after, so the item-6 site was kept.
- `dhj` (`PROBLEM=deep_hot_jupiter_rt`, spherical polar + dynamic WB + user srcs): **bitwise
  identical** baseline-vs-after in both `.bin` and `.hst`, ideal and general.

**Two harness traps fixed while doing this, both in `eos_compare.py`:**
- `tab_blocks` forced `data_format=%24.16e` on any tab/hst block, but AthenaK refuses to
  create a command-line parameter the input file does not define — `.hst` blocks have no
  `data_format`. It now only overrides blocks that already have the line.
- `dhj_gen.athinput` carries the PHYSICAL `Rgas = 4.593e7` (mu ~ 1.81), which the general
  branch ignores (composition lives in the EOS; the placeholder gamma law has mu = 1). The
  two runs then build different atmospheres and disagree at ORDER UNITY. The case now
  overrides `problem/Rgas=83144621.4563013`, and agrees exactly. Same trap as
  `solar_convection` — see [[eostest-input-inventory]].
- Earlier trap, still true: the inputs write `.tab` at `%12.5e`, hiding anything below ~1e-5
  relative. Always force `%24.16e`. The residual ideal-vs-general gap of ~1e-15 on linear
  waves is pre-existing and expected (general reconstructs `p` through PLM; ideal recomputes
  it from the reconstructed `e`).

## Left deliberately undone

`WBBackgroundStencil`'s five `Pressure()` calls could be handed back for free by `WBAdvance`,
but the `p -> e -> p` round trip is not bitwise neutral and well balancing lives on those
bits. There is a `TODO(stage3)` at the call site.

## State / next

- Build dir is currently `PROBLEM=deep_hot_jupiter_rt`.
- **The Stage 3 composition-model decision is still open** and still blocks the analytic EOS:
  Saha per cell (now ~1 solve/cell/stage, so ~2-5x total runtime), fitted ionization
  fractions, or a tabulated `(log rho, log T)` EOS (which removes the nested root find).
  Radiation pressure on/off and the fate of the `general == ideal` regression tests are also
  unanswered.
