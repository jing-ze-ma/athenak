---
name: index-active-m1-he4
description: Full original index lines for the active M1 radiation and He4 presupernova work.
metadata:
  type: reference
---

## Active work (M1 radiation, He4) — full lines

- **[He4 session state 09-21 + next steps (porous envelope drains in every variant; next = rt_kappa_hsmooth test in tests_r28; he4-presn-global 12d3644f now merged into rt-integration and frozen)](he4-session-state-2026-09-21.md)**
- **[M1 RT PLAN (09-21): design note docs/dev/rad_m1_design.md; comoving-frame sources (= IDORT/Krumholz), AP-HLL, RSLA useless at depth, EOS gas-only](rad-m1-design.md)**
- **[M1 STAGE 3 SURVEY (09-21): no published implicit two-moment RT on GPUs; ARK-RT the closest; plan = matrix-free BiCGStab on Schur-reduced E + column/ADI line preconditioner](rad-m1-implicit-gpu-survey.md)**
- **[He4 RADIATION REDESIGN agreed 09-19 night: EOS weight density-only with a very low window + optical-depth-gated FORCE in the two-stream; full design + IC conversion steps](he4-rad-eos-force-split-design.md)**
- **[He4 3-D wedge cools/contracts after ~2.9 turnovers = radiative POROSITY (L_out/L 1.15 at 3.5; per-column flux reproduces F_2s to 1-5 %); resolved; steady vs runaway open](he4-porosity-luminosity-excess.md)**
- **[He4 Strang vs unsplit differ because the MLT closure is SEEDED before the two-stream runs (unsplit only); Strang is the correct one; not fixed](he4-strang-vs-unsplit-mlt-seed.md)**
- **[He4 GLOBAL (09-19 evening): expansion THERMALLY driven, user: treat as physical; fixes mlt_split_deposit + cs_max (27e0c387); NEW sp WEDGE 90x90 deg 96x256^2 periodic theta/phi with ADI transverse ported (47f80807), BUG FIXED: rt_rad_force transverse term read unfilled ghost T_g at MeshBlock edges (cs too; earlier 3-D KEh growth suspect); wedge2 clean to 2.07 turnovers (dt 3.2 s constant, no convection from the v1 seed); RUNNING tests_r11/wedge3 job 11850864 = entropy seed (vpert_sp_rand, a61b29a8, athena_v18); global 32^2 w3d job 11847855 stays 1-D-like; all unpushed; 6 retractions inside](he4-r11-floor-arms-verdict.md)**
- **[PRESUPERNOVA He STARS sized 09-17: none fits a box (FeCZ = 0.4-0.7 R, super-Eddington); needs global spherical RHD; numbers in he-presn-sizing](he-presn-sizing.md)**
- **[He4 PRESN GLOBAL MODEL PLAN (09-17): base red_giant.cpp; 2 blockers (cs ADI fatal, no r^2 in two-stream deposit); grid/cost/build order](he4-presn-global-plan.md)**
