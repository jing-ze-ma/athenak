---
name: w121-observables-lit-0929
description: WASP-121b observational comparison targets (phase curves, limb winds, energy budget) - survey in /viper/ptmp2/jinma/w121prod_0929/obs_lit/LITERATURE.md (09-29)
metadata:
  type: reference
---
Survey (31 refs, [unverified] marks snippet-only numbers): /viper/ptmp2/jinma/w121prod_0929/obs_lit/LITERATURE.md.
Key targets: JWST NIRSpec G395H phase curve (Mikal-Evans+2023: day 3924/4924 ppm, night T_b 926/1122 K, peak ~3 deg east);
NIRISS/SOSS energy budget (Splinter+2025: T_day 2717, T_night 1562 K, A_B 0.28, eps 0.25); Spitzer offsets ~0 after
Davenport+2025 reanalysis; HST (Mikal-Evans+2022: ~9 deg east); ESPRESSO limb winds (Seidel+2025: Fe morning -4.1,
evening -6.9 km/s; Na jet 13.7 -> 26.8 km/s; H-alpha +5.2 / -19.2); IGRINS phase-resolved (Wardenier+2024; H2O prefers
strong drag); Sing+2024 Kp 215.7, dayside wind ~5 km/s. Our runs: no scattering (A_B = 0) -> T_day ~8 % hot expected.

**SYNTHETIC OBSERVABLES rot 300 1x (09-29, synth_rot300/RESULTS.md, verified phase.txt; olr_dump diag branch dhj-olr-dump 2efb2851, pushed by the agent):** peak offset 17-21 deg EAST in every band (obs ~3 deg), nightside too bright (NRS1/NRS2 night 845/1291 ppm vs 136/630; eps 0.35 vs 0.25), dayside NRS2 667 ppm low, limb Fe RV +2.2/-13.2 km/s (obs -4.1/-6.9; limb mean -7.3 ~ obs -5), Na eq +7.4/-15.7 (Seidel jet 13.7/26.8). Albedo 0.277 IMPOSED by the input (not predicted). Likely missing: drag (jet 14 km/s), nightside clouds, (dayside: metallicity -> w10x). 

**10x rot 300 (09-30 ~00:20, ana_rot300_10x + synth_rot300_10x):** 10x does NOT fix the dayside NIRSpec deficit (NRS2 -641 vs -667 ppm; 3.5-4.4 um ck band stays 2729 K -> coarse ck band / missing opacity). 10x cuts nightside ~3x (NRS1/NRS2 night 278/539 vs obs 136/630), offset 6-6.5 deg E (obs ~3), eps 0.13 (obs 0.25) = overshoots; T_day/T_night 2789/1337 (obs 2717/1562). Limb RVs unchanged (evening Fe -12.2 vs -6.9). Deep converged (100 bar -0.09 K/10 rot). KE still rising +7 %/100 rot. hst: unbooked +7.6e27 (as 1x, but 10x has OUTFLOW, so not inflow-tied) + flat tot-E loss -3.8e28 erg/s (1.35 % L) every window; T isobars flat -> bookkeeping/ck gap, not diagnosed (ck_impl_tol 1e-7 in 10x vs 1e-8).

**10x LEDGER CLOSED 09-30 ~01:10 (budget_rot300_10x A/B/T, residual 2.5e20):** +7.6e27 unbooked = ck clamp reporting (CKX; here mostly the LINEARISED cadence steps, shell 72, p < 1e-7 bar); ck_impl_conserve=1 cuts it to 9.1e26 (linearised classes not covered; optional extension). The -3.8e28 flat loss is REAL deep cooling (10-430 bar, ~-0.9 K/100 rot): upward fluid enthalpy flux 1.15e28 (wall) -> 4.4e28 (10 bar) vs Lrad_bot 9.2e27 drains the initial deep adiabat (effective T_int ~1.5x input); 1x deep warms slightly instead. ck_impl_tol 1e-8 vs 1e-7 irrelevant. Open: is the deep upward enthalpy flux resolved circulation or numerical mixing (not tested).

**METALLICITY vs DRAG literature 09-30 (obs_lit/METALLICITY_VS_DRAG.md, tagged verified/abstract/unverified):** effect known (Showman+2009 HD189733b 1x->5-10x offset 50->20 deg; Lewis+2010; Kataria+2015; Drummond+2018; Zhang&Showman 2017); consensus: helps but insufficient (Parmentier+2021, Komacek+2017); UHJ magnetic drag expected (Perna+2010, Rauscher&Menou 2013, Beltz+2022, Kesseli+2024). GAP: every published WASP-121b GCM-data comparison is 1x solar varying only drag (Wardenier+2024, Frazier+2026, Davenport+2025) while retrievals give ~10x (Pelletier+2026, Evans-Soma+2025) -> drag verdict conditional on solar metallicity; a matched metallicity x drag study is new. Discriminators: Doppler winds (Seidel+2025/2026, Vaulato+2025 dKp -15+-3), offset vs wavelength (Yang+2026 4->9 deg with depth). Recommend MHD scan at ~5-10x.

USER 09-30 ~07:15: plan the 3x/10x MHD drag run -> PLAN-2026-09-30-w121-mhd-drag-metallicity.md (6 arms {3x,10x} x {0,3,10 G}, >= 10-20 rot, deps: max_eta scan, 3x state choice, machine, GO).
