# NOTE 2026-10-17 raven: accretor 3-D PHYSICS groundwork (problem/phys_3d), static 3-D tests

USER approved 10-17 (tests only). Mesh stays spherical-polar (USER: "still use sp not cs").
Code: branch **accretor-3d-1017 @ 5825a887** (github jing-ze-ma/athenak; PUSHED, NOT merged), off accretor-bigdom-1015
(99adf52d); only src/pgen/ry_per_accretor.cpp + design doc **section 16**. `git merge-tree` with accretor-sector-light-1017
(62119e93): clean. Raven dir /raven/ptmp/jinma/acc3d_1017 (NOTE.md, decks run/, jobs/, tools tool/).

## What is implemented (all keys read only when given; 2-D decks run the old code)
- `problem/phys_3d = true` (3-D mode, envelope, thermo = general): exact 3-D Roche potential RochePot3 (same as the 3-D lobe
  mask) in phicc0 / phi0 (all faces, etotgrav), Phi_wb (PhiWB3, spin term with r sin theta), g_r remainder, g_theta, g_phi on
  (m,k,j,i); full Coriolis -2 Omega x v (2-D band had a_r = 2 Omega v_phi, a_phi = -2 Omega v_r and the sin theta = 1
  centrifugal for every row); centrifugal inside Phi, no explicit centrifugal work; IC (WB march per (m,k,j) column), T_col,
  G(psi), a_ref, F_col, rho_amb, relaxation targets on the 3-D equipotentials; Lz bookkeeping with R = r sin theta.
- `problem/theta_bc = reflect | outflow | hse` (outflow/hse: mesh ix2/ox2_bc = user): theta ghosts walk the column's net
  potential Phi + G(psi) with the local polytrope along theta; outflow = v_theta diode, hse = closed (mirror). Radiation at the
  theta faces stays reflecting (a Marshak/dark theta face drains the deep envelope: first smoke broke in 50 cycles). Theta-face
  mass/Lz/energy booked in Mout/Jout/Ebnd. `problem/theta_sponge_deg/_tau`: wall-layer velocity sponge.
- `problem/dt_report`: global MPI MINLOC over ranks (+ [rank, theta]).
- From the WASP-121b pgens (deep_hot_jupiter_rt.cpp): the spherical-basis Coriolis, rot_potential + etotgrav (centrifugal in
  the potential, the double-count guard), the equipotential (ZEffFromPot-style) IC, the sponge form. Nothing cs-specific.

## Bitwise gate (A100, keys absent, vs v4 5b0d4bc1): IDENTICAL hst/bins/cycle-dt/dt_report lines
v4 gate A (base deck, donor_irr=false, 20 cyc) for all three code versions; 2-D cal deck (20 cyc); b3d_half 3-D plumbing deck (10 cyc).

## Static 3-D tests (stream off, donor_irr off; 420 x 58 x 512, +-20.4 deg, 2 A100 per arm; 2-D band reference on the same r/phi grid)
| arm | walls | t reached | result |
|---|---|---|---|
| s3dA | reflect | 0.0083 | 0 NC/retries; deep envelope v_theta up to 8 km/s at the walls (wall face misses the hydrostatic pressure) |
| s3dH / s3dO | hse / outflow+sponge (v2 ghosts) | 0.013 / 0.015 | outflow walls leaked deep envelope (r 6.3-7) mass, 2.6e-4 of Mdom by 0.015 |
| s3dO2 | outflow, no sponge (v2) | 0.012 | = s3dO (sponge irrelevant) |
| s3dH3 / s3dO3 | hse / outflow (v3 ghosts) | 0.037 / 0.037 | leak cut: |dMdom/M| 2e-5 / 1e-4 by 0.037; dt 1.41e-5 -> 1.29e-5 |
- Away from the walls the 3-D star = the 2-D band: equatorial rho median |3D/2D-1| 6e-4 (t 0.0075), 1.4e-3 (t 0.03); photosphere
  radius over phi identical; thin-atmosphere flows the same (max|v| 123 vs 118 km/s at t 0.03). Energy budget closes; mass to 1e-13
  (reflect).
- **OPEN: theta-wall pile-up.** In all wall types the photospheric/atmospheric gas within 4-6 rows (3-4 deg) of the walls piles up
  (x5-18 by t 0.0125, x10-80 by 0.035 at rho0 1e-6..1e-3, worst at phi ~0 where the equipotentials tilt most): the atmosphere
  drifts poleward (~1.5 km/s) and converges at the walls. Not the hydro wall type; suspect the reflecting radiation wall / lateral
  radiation force (F_theta = 0 at the face). Not production-ready until solved.
- Cost: 0.72 s/cycle per 12.5 M cells on 2 A100 (0.12 GPU-s / M cell-cycle; 049c1f39 benchmark 0.13): phys_3d is free. Half-res
  (49 M, 4 A100) stays ~15 node-h per orbit pre-impact.

## Stream-on plumbing smoke (outflow walls, stream + donor_irr s=1, 901 cycles to t 0.012, 2 A100)
0 NC, 0 retries; Mnzf 254.88 (+0.015 %); Mmsk 4e-8; dt 1.41e-5 -> 1.15-1.3e-5 (thin floor gas near the photosphere). The stream
is now vertically confined: sigma_z 1.005 (nozzle) -> 0.965 (r 14) -> 0.81 (r 13).

## Open decisions (USER)
1. Theta-wall pile-up: next step options (a) radiation outflow-type theta face (needs an implicit-BC option), (b) a wider wedge
   with walls well above the stream, (c) walls that are closed inside the star and open only for the atmosphere.
2. Wall type for production: hse (closed, smallest leak) vs outflow (lets splash/stream leave; books Mout).
3. env_lin_band_factor in the wedge: kept 0.963286 (as b3d_half; reflecting radiation walls), 1 only for a full sphere.
4. donor_irr_scale in 3-D: default 1 (physical dilution); the 2-D calibration 0.22 does not apply.
Viper: v3 also built with ROCm 7.1 (/viper/ptmp2/jinma/acc3d_1017/bin/athena_apu_ry3d_v3_r71, md5 b9229a3f), not run.
