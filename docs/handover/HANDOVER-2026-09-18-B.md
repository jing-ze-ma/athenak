# HANDOVER 2026-09-18 B (viper, 10:30): productions complete, He4 global model at 2 turnovers

Read after HANDOVER-2026-09-18.md. Memory snapshot: the local memory dir is current (MEMORY.md
"CURRENT STATE 2026-09-18 10:30"); no new copy made under docs/handover (inode economy).

## Productions (both COMPLETE at 40 turnovers)
- B star prod_w7: 40 node-hours. Steady since turnover 16. Final maps: analysis_0916/patterns_0918
  (dark, absolute scales, axes in R_star): rms v_z 0.34 km/s (tau 1), 2.01 km/s (tau 1800) =
  0.82 v_MLT; Mach 0.013/0.028; rms dT/T 1.0e-3/9.5e-4 (dT/T ~ Mach^2).
- He box_w8: 62 node-hours. Driving-layer velocity came DOWN from 2.5 v_MLT (27 turnovers) to
  1.12 v_MLT (40): the overshoot was a transient. rms v_z 0.08/0.21 km/s (tau 1/30), Mach
  0.0022/0.0036, dT/T 9e-6/2.3e-5.
- Light curves / SLF spectra (analysis_0916/lightcurve_0918, RESULTS.md, RESULTS_he.md): B box
  rms 12 ppm, alpha0 3.5 ppm, nu_char 4.0 /d, gamma 3.1; disk-scaled (N_eff 280-560) alpha0 0.2 ppm,
  30-100x below the observed B0.5-B1 V stars (Bowman 2020 / Burssens 2020: HD 36960 28 ppm,
  HD 37042 53 ppm, HD 43112 5.8 ppm; nu_char 0.8-1.8 /d). He: 0.02 ppm, drift-dominated.
  Reasons (in order): ZAMS star = weakest FeCZ (Schultz+2022 ZAMS vs mid-MS differ 1000x); the
  box lacks scales > 0.1 R; dT/T ~ Mach^2 and radiative erasure; other sources in real stars.
  0.2 ppm is below every detection floor (TESS C_w 3-8 ppm).
- Pages: production https://claude.ai/artifact/SLgYFj9JDBUsNqxj5UYgPJ (v13), HRD of 3-D envelope
  sims https://claude.ai/artifact/79q5ATcM6xcDy5QLTR39Gz (v10; files analysis_0916/hrd_3d incl. dark
  build-up frames), RG prod11 dark https://claude.ai/artifact/EHrBVrgbdYKCuh4rVaSn9z (from the
  embedded rasters; run data on orion).

## Global 4 Msun presupernova He star (he4-presn-global, bench/wt_he4, unpushed)
Since HANDOVER-2026-09-18.md: tests_r10 (top-face measurement), tests_r11 (rt_top_vacuum exposed:
NULL; mlt_alpha 3: worse; 3-D one-piece run g3d to 2.2 turnovers: alive, KEh x30/turnover,
L_out 0.97, dt -> 0.07 s). dt killer = hot VOIDS: ~1 % of columns at 0.96-0.99 R above the swept
shell reach rho 1e-12, T 1e7 K, v 6e7 (not seam rows; rt_de_max caps only explicit RT updates).
RUNNING: tests_r11/d1 (dfloor 1e-11, job 11798237), d2 (dfloor 1e-10, 11798238), 2 turnovers each,
script r11_3d_arm.sh (rst every 0.2 turnover). The IC top density is 5.5e-11.
NEXT: evaluate d1/d2 (dt history, void count, L_out/L, KEh), chain the better one to 5 turnovers,
judge settling and convection onset; production only with the user's approval. Then: port
4bcdc855 + 94c7165d to rt-integration; bvals_fc seam defects (MHD); push after user OK.
Opus agent quota returns 2026-09-19 23:00 (Brussels).
