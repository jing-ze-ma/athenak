---
name: fecz-low-lum-gap
description: Literature check 2026-09-13 - NO multi-D simulation of the iron-bump convection zone (FeCZ) exists for main-sequence B stars of 8-20 Msun (L/L_Edd < 0.15). Lowest published - Debnath+2024 (2D, O8 26.9 Msun, Gamma 0.16), Schultz/Bildsten/Jiang 2022 (3D global, 35 Msun ZAMS, Gamma ~0.18), Jiang+2015 boxes at Gamma 0.29-0.57 (80/40 Msun). Pathak+ PPMstar 25 Msun moved the Fe bump inward. Cantiello+2021 ran no new sims. Observations - SLFV ubiquitous (Bowman 2019/2020); Anders+2023 core IGWs too weak (0.06 vs >10 umag); Pedersen/Bowman 2025 favour the FeCZ but SLFV persists below the ~16 Msun FeCZ threshold and in the SMC = the low-L regime is exactly the open test
metadata:
  type: reference
---
Candidate box project - low-L B-star FeCZ at 32/H_p with the general EOS (eos_radiation on), solar opacity table, two-stream RT,
WB + lhllc (needed, see [[rg-box-sweep]]); hours per model on 2 GPUs. Refs - Jiang 2015 arxiv 1509.05417; Schultz 2022
2110.13944; Debnath 2024 2401.08391; Anders 2023 2306.08023; Pedersen/Bowman 2025 2504.15861; Cantiello 2021 2102.05670.

HELIUM STARS (same search): also NO multi-D iron-zone convection sim for low-mass stripped He stars (3-10 Msun, Gamma < 0.2).
The only multi-D He-star work is Moens+2022 (arxiv 2203.01108): one 10 Msun He star at log L 5.33-5.74, Gamma 0.33-0.84,
box-in-wind FLD, wind launching, convection only qualitative; Goldberg+2025 (2508.12486) partially stripped YSGs are H/He
opacity peaks not Fe. Everything at low Gamma is 1D (Grassitelli, Ro & Matzner, Poniatowski). Needs a He-dominated EOS +
opacity table first.
