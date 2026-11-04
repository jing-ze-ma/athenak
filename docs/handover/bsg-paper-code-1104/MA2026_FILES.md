# Ma+2026 BSG files (from L. Ma) vs our BSG setup — 10-06

Inputs: `ma2026_files/BSG/` (Linhao's IC + TOPS tables), our arm-2 run
(`/raven/ptmp/jinma/bsg_raven_1002/arm2r2b_4n/bsg3d_arm2.athinput`, IC
`/raven/u/jinma/bsg_raven_1002/bundle/ic_bsg_arm2.txt`, tables
`rosseland/planck_tops_x0.7_z0.008.txt`, `he_opac_logd_min = -21`), SK -67 133 column
`sk67/ic0/ic_sk67_re_f20.txt`. Script `ma2026_files/ana/ana_ma2026.py` (login node, CPU);
every number below is in `ma2026_files/ana/numbers.txt`. Plots: `ana/ratio_maps.png`,
`ana/ic_compare.png`, `ana/hse_ma.png`, `ana/thin_top_kPkR.png`.
Nothing in any run or existing table was touched; no git changes.

## 1. Opacity tables: same source, same mixture
- Their grid: log T 4.0646..6.5417 (34 nodes, uneven: TOPS keV list), log rho -11.523..1.0
  (100 nodes, step 0.12649). kappa_e = 0.34 except 499 nodes where kR < 0.34+1e-5
  (kappa_e = kR - 1e-5, paper eq. 7 rule; min 0.3355). kP/kR spans 4.5e-8..2350.
- All 34 of their T nodes are on our raw TOPS T list (`opac/tops_gs98_x0.7_z0.008.dat`).
  Our raw TOPS interpolated to their nodes: |log ours/theirs| median 0.0013 (kR), 0.0017 (kP);
  95 % 0.015/0.016. So it is the same TOPS query, same X=0.7 Z=0.008 GS98 mixture; the
  residual is rho-node interpolation (their rho nodes are not ours).
- Our AthenaK tables (resampled to 0.025 x 0.05 dex): median 0.0021 (kR) / 0.0078 (kP),
  95 % 0.018 / 0.060 dex; maxima in irrelevant corners (kR 0.14 dex at log T 4.07,
  log rho -0.26; kP 0.18 dex at log T 6.0, log rho -11.5). Map: `ratio_maps.png`.
- Along the columns, inside their box: kR within 0.029 dex (7 %) everywhere (max at the Fe
  bump, 40-41 Rsun); kP within 0.055 (Ma col, 44.8 Rsun), 0.064 (arm 2, 45.3), 0.097 dex
  (SK67, 41.6 Rsun). The kP differences come from our coarse resampling of the spiky kP.
- Their column's MESA (OPAL) kappa vs their TOPS table: median 0.011 dex, max 0.063 at
  26.6 Rsun — the "consistent with MESA" claim holds to ~15 %.
- Physics mapping: Athena++ uses kappa_a = kR - kappa_e (absorption), kappa_e (scattering),
  kP (thermal exchange). Total transport opacity = kR, energy exchange = kP. Our rad_m1
  table mode uses kappa_F = kR (kappa_s = 0) and kP for emission and absorption — the same
  thing. The electron table is not needed by our module.

## 2. Edge treatment and the thin top (the key point)
- Their files carry no extension; the AthenaK lookup (and, by inference, the standard
  Athena++ star pgens) CLAMPS T and rho to the table edge. Below T 11.6 kK and rho 3e-12
  both clamp to the corner node (log T 4.0646, log rho -11.523): **kR = 0.441, kP = 124,
  kP/kR = 281, constant** over the whole thin top.
- Ours (TOPS to log T 3.76 and log rho -14, then the -21 extension, slope clamped [0,1]):
  | T, rho | ours kR, kP, kP/kR | theirs clamped |
  |---|---|---|
  | 10 kK, 1e-12 | 0.41, 165, 406 | 0.44, 124, 281 |
  | 10 kK, 1e-14 | 0.31, 4.4, 14 | same |
  | 10 kK, 1e-16 | 0.30, 0.41, 1.4 | same |
  | 10 kK, 1e-18 | 0.29, 0.038, 0.13 | same |
  | 7.9 kK, 1e-18 | 0.28, 0.027, 0.095 | same |
  Our kP valley at rho 1e-13: kP 105 (7.1 kK) -> 22.9 (9.4 kK) -> 56 (11.2 kK); theirs 124
  flat. Along our arm-2 IC (floor region rho 1e-20, T 10.4-11.6 kK, 55-64 Rsun): ours
  kP/kR 6.4 -> 0.087, theirs 284 -> 281; our kP 1.9 -> 0.025 vs 125.
- Near the photosphere nothing changes: arm-2 IC tau_R = 10/3/1/0.3/0.1 (49.8-51.0 Rsun,
  T 26-13 kK, rho 3e-10..1e-11) all lie inside their box, kP/kR ours vs theirs 189/189,
  435/432, 532/523, 399/395, 386/382; tau 0.01 (51.45 Rsun, rho 9.6e-13, outside in rho):
  327 vs 345. Our arm-2 column leaves their box at r = 51.3 Rsun (rho < 3e-12, tau ~0.03)
  and 55.3 Rsun (T < 11.6 kK).
- So all our thin-top trouble (gas-T root cycle, kP valley 9-12 kK, rho < 1e-14) sits in the
  region where Ma+2026 had a constant, strongly emissive opacity. Gas-radiation thermal
  relaxation rate c kP 4 a T^3 / c_v (per mass, rho-independent; hand estimate, c_v = 1.5 k/(0.62 m_u)) at 11 kK: ~750 /s with
  kP = 124 vs ~0.15 /s with our kP = 0.025 at rho 1e-20: their top is thermally locked to
  T_rad, ours is nearly decoupled and kP/kR < 1 makes it scattering-like.

## 3. IC comparison (`ic_compare.png`, `hse_ma.png`)
- Their file: 327682 rows = two segments (16.50-19.99 and 19.99-52.49 Rsun, one duplicate
  row). F_r = L/(4 pi r^2) with L = 5.58146e38 (log L/Lsun 5.1638) to 2e-12: this is the
  paper's RADIATIVE-ONLY profile (rad. diffusion residual with their kappa column 2.5e-4).
  Spacing coarsens outward: dr > 0.01 Rsun beyond 36.3, > 0.1 Rsun beyond 45.3; the last
  step is ONE interval 48.32 -> 52.49 Rsun (T 70 -> 63 kK, rho 3.4e-9 -> 2.6e-11). The
  file is reliable only to ~48.3 Rsun; it contains no photosphere (T 63 kK at the top;
  Teff for L and 52.49 Rsun would be 15.6 kK).
- L: ours arm 2 5.5459e38 (-0.64 %), SK67 5.5587e38 (-0.41 %).
- At fixed r (Ma / arm 2 / SK67): 20 Rsun T 1.077/1.053/0.953 MK, rho 6.1e-5/5.6e-5/4.0e-5;
  40 Rsun 193/174/112 kK; 46 Rsun 98/73/15 kK. At fixed T, rho agrees within 2-8 %
  (31 % at 200 kK, Fe-bump) between Ma and both of ours; the radius of a given T is larger
  in theirs (T = 1e5 K at 45.9 / 44.2 / 40.6 Rsun): the radiative-only column is inflated
  above the Fe zone. SK67 (R_ph 45.9 Rsun) is the most compact.
- Convective zones (nabla_rad with our kR, L, GM 20 Msun > nabla_ad of gas+radiation):
  Ma 39.4-45.2 Rsun (radiative-only, superadiabatic, peak nabla_rad/nabla_ad 1.2; the He
  zone would sit in the unresolved top); arm 2 36.3-43.6 and 48.5-49.6; SK67 33.9-40.2 and
  44.1-45.3.
- HSE of their rho, T with OUR EOS (mu 0.62 ideal + a T^4/3): with their m(r) mean
  residual -0.70 % (20-30 Rsun -0.70, 36-46 -0.35, 46-50 -0.68, max 1.3 %); with the
  simulation's point mass 20 Msun -1.2 % (m(r) = 19.82 Msun at 20 Rsun, 19.998 at the top).
  With mu 0.616 (fully ionised X 0.7, Z 0.008, MESA-like) the m(r) residual drops to
  -0.15..-0.37 %. So their IC is ~1 % out of HSE in their own Athena++ setup (mu 0.62 +
  point mass), consistent with the large initial transient they report. The 48.3-52.5
  interval fails HSE by 51 % (single coarse step).

## 4. Conversion (`ma2026_files/athenak_tables/`)
- `rosseland_ma2026_x0.7_z0.008_padded.txt`, `planck_ma2026_x0.7_z0.008_padded.txt`
  (+ `electron_..._padded.txt`, not used by our module): our he_star_m1 format (comment
  lines, `# nT nD lTmin dlT lDmin dlD`, nT*nD values, T slowest). log T resampled to a
  uniform step 5.00024e-4 (6085 rows incl. pad), their 100 rho nodes kept exactly; padded
  with EDGE values down to log T 3.4996 and log rho -22.02 (183 columns) so the startup
  "column inside the tables" check passes; the pad reproduces the clamp exactly.
  Lookup error vs bilinear on their native grid: kR max 2.4e-4 dex, kP max 2.7e-3 dex
  (99.9 % 1.7e-4). Unpadded versions (`*_ma2026_x0.7_z0.008.txt`) are there too but would
  FATAL with our arm-2/SK67 columns (rho 1e-20, T 10 kK lie outside).
- Input keys (exact Ma+2026 opacity):
  ```
  <problem>
  he_opac_table   = /viper/ptmp2/jinma/bsg_1001/ma2026_files/athenak_tables/rosseland_ma2026_x0.7_z0.008_padded.txt
  he_planck_table = /viper/ptmp2/jinma/bsg_1001/ma2026_files/athenak_tables/planck_ma2026_x0.7_z0.008_padded.txt
  # REMOVE he_opac_logd_min (and he_opac_extend_hold): the pad already is the clamp
  ```
  (<rad_m1> opacity = table and the cgs units as in arm 2; EOS ideal mu 0.62 gamma 5/3 and
  GM 20 Msun already match the paper.) Not smoke-tested here (analysis-only brief).
- IC: `ic_ma2026_r20.txt`, he_ic_cols = 5 format (r rho eint F_r E=aT^4 | T fmlt=0),
  eint = rho k T/(0.62 m_u)/(2/3), 169433 rows, r 19.7125..52.4863 Rsun, no atmosphere.
  The reader needs the file to cover the mesh top (plus ghosts), so it cannot be used
  until an atmosphere is appended; and it needs mlt_flux_frozen = false (fmlt = 0).
- Their atmosphere (paper sect. 2.2.3): above r_photo keep T = T_photo and
  rho = rho_photo exp(-(r - r_photo) g / P_photo) to the density floor, v = 0, E = a T^4,
  F_r = L/(4 pi r^2). As printed the exponent is dimensionally P per unit density
  (i.e. g/(P_ph/rho_ph)). Building it needs r_photo, T_photo, which this file does not
  contain: either (a) take the file top as "photosphere" (52.49 Rsun, 63 kK, P = 4.09e4,
  almost all radiation: P_gas = 0.22) — then with total P the scale height is ~112 Rsun
  (density hardly falls to 80 Rsun), with P_gas 4e7 cm; or (b) integrate the radiative
  column outward from ~48 Rsun with the TOPS kR, mu 0.62 + radiation, L 5.581e38 to
  tau = 2/3, then attach the isothermal layer. (a) vs (b) is the open question for Linhao.

## 5. Recommendation
- Switching to their tables changes nothing in the interior (same TOPS data; <= 7 % kR,
  <= 25 % kP at isolated kP spikes, mostly our coarser resampling). For the exact Ma+2026
  setting use the padded tables above, no he_opac_logd_min. The real change is the thin
  top: kP/kR = 281 flat instead of ours 0.09-560 with the 9-12 kK valley.
- Expectation for the thin-top problems: the clamp makes the top strongly thermally coupled
  (gas T locked to T_rad, no kP valley, no kP/kR < 1 scattering-like layer), which removes
  the opacity side of the gas-T root cycle. Fix A / implicit_realisable_coupling acts on
  the unclipped flux iterate and may still be needed: test both on the gd-physfix-1005
  binary (one apudev smoke per arm). Note this is opacity regularisation of the atmosphere
  by fiat, just like he_opac_extend_hold but with the edge at 11.6 kK / 3e-12 instead of
  our -14 rho edge.
- vet_gd: its source function uses the capped kP/extinction as absorption fraction; with
  the clamp the whole top becomes S = a T^4 (fraction capped at 1) instead of partly
  scattering, so vet_gd thin-region gates must be rerun on the new tables.
- SK -67 133 IC: rebuild `make_ic_mlt_star.py` with the padded Ma tables (the interior
  barely moves; the atmosphere part above rho 3e-12 changes kP). Keep r_in = 20 Rsun
  ABSOLUTE as in the paper (their criterion: 5-8 H_p below the CZ lower edge at ~36 Rsun).
  20 Rsun is 0.40 of the 1-D R_ph ~50 Rsun; for SK67 (R_ph 45.9) it is 0.44 R_ph, base
  T 0.95 MK (vs 1.08 MK in theirs), CZ starts at 33.9 Rsun, so the margin is at least as
  large; the same 0.40 R_ph would be 18.4 Rsun.
- EOS: ours already matches (ideal mu 0.62, gamma 5/3, GM 20 Msun point mass).
