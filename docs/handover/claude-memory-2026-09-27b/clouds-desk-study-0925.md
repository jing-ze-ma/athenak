---
name: clouds-desk-study-0925
description: 09-25 cloud desk study (/viper/ptmp2/jinma/clouds_0925) -- ADAM = renamed SPARC/MITgcm (11 Kataria bands, scattering Toon89 + delta-DO); our RT lacks scattering; cloud-free baseline defensible for WASP-121b; minimal test = diagnostic nightside MgSiO3 deck
metadata:
  type: project
---
ADAM ("ADvanced Atmospheric MITgcm", Mehta, Parmentier, Tan, Lee et al., arXiv:2509.23406, 2604.26911) = SPARC/MITgcm
renamed; 11 Kataria bands, premixed Lupu/Marley tables mixed by random overlap (exo_k); scattering two-stream
(delta-DO stellar, Toon+1989 thermal), explicit RT once per step. Clouds: Tan & Showman 2021 tracers (Na2S, KCl, MgSiO3),
fixed lognormal size, Batalha & Marley 2020 Mie tables; code/tables not stated public.
OUR GAP: no single-scattering albedo anywhere (Rayleigh = extinction, correlated_k.hpp ~981); beam (1-A) factor
A ~0.015 vs observed Bond 0.277 +- 0.016 (Frazier+2026, arXiv:2605.01589).
WASP-121b: JWST nightside T_b 926/1122 K; cloud-free Parmentier+2018 GCMs over-predict nightside flux; NIRISS
nightside featureless ("may be best explained by nightside clouds"); strong drag also needed (clouds and drag
degenerate on the nightside).
Effort: stellar scattering (T-independent at frozen opacity, outside Newton) 1-2 wk, low risk, +10-25 % RT; thermal
scattering LW-0 (Li 2002 scaling) 2-4 d, LW-1 VIM 2-3 wk moderate, LW-2 Toon89 4-6 wk high risk; Mie tables 3-5 d;
diagnostic cloud scheme 3-5 d; tracers 2-3 wk. Validation reference: Exo-FMS_column_ck (has a WASP-121b namelist +
DISORT/two-stream scattering solvers).

**How to apply:** cloud-free baseline is fine (every SPARC/ADAM WASP-121b model); if clouds are pursued, start with
stellar scattering + a diagnostic nightside MgSiO3 deck (r0 1 um, condensed fraction 1 and 0.1, LW-0) -- needs the user's go.
