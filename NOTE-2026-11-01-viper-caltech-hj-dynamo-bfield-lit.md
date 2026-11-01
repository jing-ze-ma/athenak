# NOTE 2026-11-01 viper -> Caltech (USER): literature on dynamo field strengths of hot / ultra-hot Jupiters, and what it means for our B seed

FYI for choosing the field strength of the MHD runs (WASP-121b and the new planets). Nothing to run.

## Bottom line
- Deep convective-dynamo predictions are mostly **10-100 G** at the top of the dynamo region (heat-flux scaling).
- Newer work (2025-2026) finds **~1 G (Jovian or less)** if the inflation heat is deposited shallow, outside the dynamo region.
- Our production MHD deck (10 G polar dipole at R_p, about 3 G at the domain bottom) sits inside that range, on the conservative side of the deep-heating picture.

## Heat-flux scaling (the standard method)
- Christensen, Holzwarth & Reiners (2009, Nature): in fast-rotating convective dynamos the field depends mainly on the convective heat flux, B ~ rho^(1/6) F^(1/3), almost independent of rotation.
- Reiners & Christensen (2010, A&A): applied to exoplanets; a few to a few tens of G. *(cited from memory, not re-checked)*
- Yadav & Thorngren (2017, ApJL): with the internal heating inferred from inflated radii, fields up to about 10x Jupiter, roughly 50-250 G for strongly inflated hot Jupiters. Sources disagree on whether they found an incident-flux dependence.
  - https://spaceref.com/?p=137005
  - https://arxiv.org/html/2411.00674v2

## Observational inference
- Cauley, Shkolnik, Llama & Lanza (2019, Nature Astronomy 3, 1128): Ca II K star-planet interaction signals in 4 hot-Jupiter systems give surface fields of **20-120 G**. That is 10-100x rotation-based scaling and consistent with heat-flux scaling.
  - Model-dependent; no direct measurement of an ultra-hot Jupiter field exists; radio confirmation is pending.
  - https://arxiv.org/abs/1907.09068

## Pushback, 2025-2026
- *Rossby number regime, convection suppression, and dynamo-generated magnetism in inflated hot Jupiters* (2025): if the inflation heat sits in the outer envelope, convection in the dynamo region can shut down, and flux and field drop to Jovian levels or below. Dynamo regime transition at Ro ~ 0.12.
  - https://arxiv.org/abs/2507.05202
- Elias-Lopez et al. (EPSC 2026 abstract): deep heating up to ~100 G; shallow heating ~1 G.
  - https://meetingorganizer.copernicus.org/EPSC2026/EPSC2026-1242.html
- Coupled dynamo + Ohmic-atmosphere model (2025): strong Ohmic dissipation shrinks the convective dynamo region and can make the field oscillate.
  - https://arxiv.org/html/2507.13991v2
- KELT-9b induction study (2022): same scaling, calibrated on Jupiter; the internal heat flux is the main uncertainty.
  - https://arxiv.org/pdf/2210.03351

## Atmospheric (shallow) dynamos, a separate mechanism
- Rogers & McElwaine (2017): a self-consistent dynamo in the hot, partly ionized atmosphere of an ultra-hot Jupiter when conductivity follows temperature.
  - https://ar5iv.labs.arxiv.org/html/1704.04197
- Boening, Dietrich & Wicht (2024): MHD of an idealized ultra-hot Jupiter atmosphere at ~2400 K; subcritical dynamo action, westward hotspot offset.
  - https://arxiv.org/html/2407.12434v1
- I found no anelastic dynamo model of an inflated hot-Jupiter interior. Closest: Duarte, Gastine & Wicht (2013) on Jupiter-like dynamos.
  - https://arxiv.org/abs/1210.3245

## Applied to WASP-121b (our 3x IC)
- T_int = 567 K gives F_int = 5.9e6 erg/cm^2/s, about 1000x Jupiter's ~5.7e3 (T_int ~ 100 K).
- F^(1/3) scaling then gives about 10x Jupiter's field (4-14 G) at the top of the dynamo region: roughly **40-140 G**, falling off as 1/r^3 to the photosphere.
  - This holds only if that flux is convective at depth; shallow heating gives ~1 G.
- Our 10 G dipole at R_p is inside the range.
- Possible later test: a 30-100 G arm. Ohmic drag scales as B^2, and the Alfven-limited dt cost rises with B; we already see it at 3-10 G.
- For the new planets, the same scaling with each planet's T_int gives a first estimate; their packages will list T_int.
