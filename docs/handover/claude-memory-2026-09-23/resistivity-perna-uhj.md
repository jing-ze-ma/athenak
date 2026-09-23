---
name: resistivity-perna-uhj
description: "How good the Perna resistivity is for ultra-hot Jupiter atmospheres, the min_xe->max_eta input rename that breaks the dhj runs, and why a tabulated x_e would be FASTER not slower"
metadata: 
  node_type: memory
  type: project
  originSessionId: 6d2fc30e-834e-473c-b800-b2c8b589b054
  modified: 2026-08-15T16:16:13.699Z
---

Analysed 2026-08-15. AthenaK's `ohmic_resistivity = perna` (`src/diffusion/resistivity.hpp:105`,
Perna+2010 after Balbus & Hawley 2000) is thermal ionization of POTASSIUM ONLY, with
`eta = 230 sqrt(T)/x_e` (electron-neutral collisions) and
`x_e = 0.0563 (a_K/1e-7)^1/2 (T/1e3)^3/4 (2.4e15/n_n)^1/2 exp(-25188/T)`. The prefactor is
the code's `6.47e-13/1.15e-11`; 25188 K already contains the factor 2 of `chi_K/2k`.
`ResistivityKumar` (a polynomial fit) exists but is COMMENTED OUT.

**Accuracy vs a proper 6-species Saha** (Na, K, Ca, Al, Mg, Fe with real degeneracy ratios,
plus H from `eos_composition.hpp`), ratio x_e Perna / x_e full:

| T | ~0.01 bar | ~1-5 bar |
|---|---|---|
| 1500-3000 K | 0.7-2.1 | 0.7-0.8 |
| 3500-5000 K | 2.8-3.3 | 1.0-1.2 |
| 7000-8000 K | 0.7-0.3 | 0.5-0.3 |

**Good to a factor ~3 over 1500-5000 K, and to ~20% near 1 bar** — far better than a K-only
fit deserves, because the prefactor evidently absorbs the aggregate metal donation. Two
systematic failures in OPPOSITE directions: above ~3500 K at low pressure it over-predicts
x_e ~3x (eta too small) because the sqrt-Saha form never saturates, while the real metals
plateau at x_e = 7.1e-5 once fully ionized; above ~6000 K it under-predicts because there
is no hydrogen term. Also above ~7000 K the full eta drops BELOW Spitzer, so
electron-neutral scaling is the wrong physics there, not just miscalibrated.

**Already correct:** `SetResistivity` (`resistivity.cpp:130-137`) reads T from the `wtemp`
cache and asks the EOS for mu when `eos.IsGeneral()`, so under `general_eos = table` the
resistivity sees the mu swing 2.33 -> 1.26 across H2 dissociation. Note the resistivity
NEVER asks the EOS for n_e — it uses the Perna fit, which does carry the alkalis. So the
EOS's inert metals (see [[eos-selection-by-run-name-bug]] sibling note in
[[general-eos-stage3-table]]) do NOT cripple magnetic coupling; the gap is the factor-of-3,
not orders of magnitude.

**INPUT RENAME THAT BREAKS RUNS:** `min_xe = pin->GetReal(...)` is commented out at
`resistivity.cpp:60` and replaced by `max_eta = pin->GetReal("mhd","max_eta")`, which is
STRICT. The dhj athinputs still said `min_xe = 1.0e-9` and had no `max_eta`, so they abort
at startup and `min_xe` is silently ignored. Fixed 2026-08-15 in
`run/dhj_sph_rt_bd_test{,2,3}/deep_hot_jupiter.athinput` -> `max_eta = 1.0e13`. NOT an exact
restoration: `min_xe` floors the ionization fraction, `max_eta` caps the diffusivity, so the
effective x_e floor is `230 sqrt(T)/max_eta`, i.e. sqrt(T)-dependent. 1.0e13 matches the old
floor at T = 1890 K, -27%/+26% over 1000-3000 K. An exact fix would restore `min_xe` in the
CODE alongside `max_eta`. `dhj_sph_rt_hd_test` and `bd_test1` have the block commented out;
`run/ohmic_diffusion` uses `constant` + `eta_ohm_const` and is unaffected.

**Units: the dhj runs are in cgs code units** (mesh 9.44e9-1.254e10 cm, `T = p/Rgas/rho`
with Rgas = 4.593e7 only gives Kelvin if p, rho are cgs), and `eta_b` goes straight into
E = eta J with NO conversion. So `max_eta` is in cm^2/s.

**A better x_e would be FASTER, not slower** (benchmarked, single core, -O2, per cell):
Perna formula 82 ns (it is transcendental-bound: log10 + pow(10,x) + 2 sqrt); 6-species
Saha solved in-kernel with the T-factors hoisted **1384 ns = 17x SLOWER**; bilinear table
lookup **13.5 ns = 6x FASTER**. Tabulate log10 x_e on the same (log rho, log T) grid — the
dynamic range is 1e-25 to 1e-2. Timestep impact is also favourable: dr = 4.8e7 cm gives
dt_CFL 32-56 s and dt_diff = dr^2/2eta = 117 s at the 1e13 cap, so diffusion is not
limiting; the correction RAISES eta only in the thin upper atmosphere (eta ~1e9, dt_diff
~7e5 s, irrelevant) and LOWERS it ~1.4x in the cold deep gas where the cap actually binds.
Caveat: the dhj runs use `eos = ideal`, so there is no EOS table to piggyback on and a
standalone x_e table would be needed (still 6x under the current formula, just without the
free index reuse).

## DONE 2026-08-15: option B implemented, `ohmic_resistivity = eos`

Commits `b689dac1`, `b484df85`, `43b0861a`, docs in `dadd887a`. `perna` is untouched and
still selectable.

- `eos_metal_ionization` (default false) adds Na, K, Ca, Al, Mg, Fe as singly ionizing
  donors in `eos_composition.hpp`. They were ALREADY inside the inert metal lump and give
  one heavy particle each either way, so the only change to `n_tot` is their electrons —
  mu moves 0.004%, n_e moves orders of magnitude. Ionization energy is carried too, which
  is what keeps c_v consistent across the metal ionization zone.
- `x_e` is a FOURTH tabulated surface, `log10(n_e/n_tot)` (ITXE=12, ITNVAR 12->16, table
  14->18 MB). Deliberately NOT read by `EOSTable::Eval()` — only non-ideal MHD wants it.
  `EOS_Data::ElectronFraction(d,e,T)`; returns 0 for an ideal gas.
- `ohmic_resistivity = eos` gives `eta = 230 sqrt(T)/x_e + 5.2e11 lnL/T^1.5` (Spitzer term
  is the larger above ~7000 K). Selecting it without a general EOS is FATAL.
- **`eos_metal_mh` ([M/H], dex) DEFAULTS to `<problem>/met`** so metallicity is set in one
  place; `deep_hot_jupiter_rt` aborts if an explicit override disagrees. The user caught
  this — the first version had a linear `eos_metal_scale` defaulting to 1, which would have
  given 3.2x solar opacity with solar donors at met=0.5.
- **`eos_metal_condensation`** removes each donor below its own T_cond, in the published
  form `10^4/T_cond = a - b log10 p_bar + d (log p)^2 - c [M/H]`, log p clamped to [-8,3].
  **Coefficients are LITERATURE, not guessed** (commit `09e6fd95` replaced my anchored
  guesses, which were wrong by up to 300 K and over-suppressed n_e by 25x at 800 K):
  Fe + Mg2SiO4 from Visscher, Lodders & Fegley 2010 ApJ 7166 1060 eqs 2/18 (valid
  800-2500 K, [Fe/H]<=+0.5); Na2S + KCl from Morley+2012 ApJ 756 172; Al2O3 + CaAl4O7 from
  Wakeford+2017 MNRAS 464 4247. T_cond at 1 bar: Ca 2004, Al 1994, Fe 1838, Mg 1698,
  Na 996, K 801 K.
  **Are they current?** They are what virga/PICASO/Ackerman-Marley use. virga v1
  (Batalha+2026, arXiv 2508.15102) moved Fe/Mg/Al to Morley+2024 — worth 1-2% in T_cond —
  but KCl and Na2S are UNCHANGED, and those are the only two that matter, because K and Na
  dominate n_e wherever the refractories are condensing. Full equilibrium condensation
  (FastChem Cond, Kitzmann+2024, arXiv 2309.02337; GGchem) is the accuracy ceiling and is a
  Gibbs minimisation per node — overkill here.
  Because K and Na condense LAST and dominate n_e, the net effect is a cliff not a taper:
  at rho=1e-5, n_e unchanged above 900 K, 0.17x at 700 K, 0.013x at 600 K, ~0 by 400 K.
  It scales the ABUNDANCE not the equilibrium constant (scaling K would change the ionization
  ratio instead of removing the species), tapers as a tanh over 5% in T (a step would make
  the Hermite patch ring on the differenced node derivatives), and needs TWO passes through
  `Evaluate` because the curves want p and p is an output.
- Reference input: `inputs/mhd/deep_hot_jupiter_rt_eos.athinput`. NOTE two parameters
  change meaning under the general EOS: `Rgas` becomes unused, and `tfloor` becomes a floor
  on CODE TEMPERATURE (8.3145e7 * T_K), not p/rho — the ideal input's `tfloor = 400.0` was
  ~1e-5 K, i.e. no floor at all.
- Verified on 8 ranks against the real RT problem: met=0 -> n_e/n_tot 1.53e-4; met=0.5
  propagates automatically -> 2.65e-4 (x1.73 vs sqrt(3.16)=1.78 for trace Saha); explicit
  disagreeing `eos_metal_mh` aborts; `eos` + `eos = ideal` aborts; condensation on/off/
  [M/H]=0.5 all run. With `eos_metal_ionization` off, `solar_convection` reproduces the
  pinned binary's base state exactly.
- Built in a SEPARATE tree `/orion/u/jinma/ATHENAK/build_dhjrt` (PROBLEM=deep_hot_jupiter_rt)
  so `athenak/build` stays configured for solar_convection.

**Regression test: DONE** 2026-08-16, `tst/scripts/mhd/mhd_eos_electrons.py`
(commit `b93ed468`) — see [[eos-electron-regression-test]].

**VERIFIED END TO END 2026-08-16.** The earlier warning here — that "verified" meant only
the startup banner, and that the first real run pinned eta at `max_eta` — is resolved.
The cause was NOT the electron fraction: two out-of-bounds reads of `eta_b` wrecked the
resistive update, and one of them had also been making `ohmic_resistivity = constant` a
silent no-op. Fixed in `c0edbe49`, details in [[eos-xe-resistivity-capped]]. The tabulated
x_e now reproduces the composition model to 0.1% in a running fluid.

## The viper-GPU dhj input is converted (2026-08-16)

`run/dhj_pole_rt_viper/deep_hot_jupiter.athinput` now runs `eos = general`,
`general_eos = table`, metals + condensation on, and `ohmic_resistivity = eos`. Original
ideal-gas version backed up only in that session's scratchpad. Validated for 2 cycles
against `/orion/u/jinma/ATHENAK/build_dhjrt` (rebuilt): table 281x451 = 15 MB, banner
n_e/n_tot = 1.53e-4, dt = 13.41 s.

**Requires `-DPROBLEM=deep_hot_jupiter_rt`.** `deep_hot_jupiter.cpp` has 0 `IsGeneral`
call sites and still references the renamed `Mesh::use_grid_stretch`, so it neither
compiles nor supports the EOS.

**max_eta is an x_e FLOOR, and at their 1e12 it binds over much of a Teq = 2500 K
atmosphere.** `x_e >= 230 sqrt(T)/max_eta` ~ 1e-8 at 1e12. EOS x_e (metals + condensation):

| rho | 1000 K | 2000 K | 3160 K | 5010 K |
|---|---|---|---|---|
| 1e-6 | 9.0e-14 | 4.3e-8 | 3.0e-6 | 4.9e-5 |
| 1e-4 | 9.0e-15 | 4.6e-9 | 6.7e-7 | 1.3e-5 |
| 1e-2 | 7.6e-16 | 4.6e-10 | 8.6e-8 | 3.0e-6 |

So 1e12 caps everything below ~2000 K (upper atmosphere) to ~2800-3000 K (deep gas).

**Cost of raising the cap, MEASURED on that 64x64x128 polar grid** (constant eta, the worst
case where the cap binds in every cell): eta = 1e11 and 1e12 both give dt = 13.41 s, the
CFL value — diffusion never binds. eta = 1e13 gives 4.68 s. So the CEILING on the cost of
1e12 -> 1e13 is 2.9x, not the ~10x a naive smallest-cell estimate suggests (the implied
min dx is 3.06e7 cm, not the sub-1e7 phi cell at the pole). Kept at 1e12 for comparability
with their `perna` runs; 1e13 is the better physics if they can spend 3x.

Measuring this at all needs the `eta_b` fix in `c0edbe49` — before it, `constant` was a
silent no-op. See [[eos-xe-resistivity-capped]].
