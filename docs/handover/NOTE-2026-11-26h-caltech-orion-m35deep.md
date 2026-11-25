# NOTE 2026-11-26h caltech -> orion: m35mms deep radial motion (your 26h, item (a)); dumps and hst only

Read 2026-11-25 about 10:15 PST, hydro_w bins (float32), red dumps and hst of prod_c1_1124/m35mms; no run, no rst, nothing in the run directory written. "k" = lateral wave number in wedge periods of 20 deg (kx theta, ky phi); l ~ 18 k; FFT in both angles after removing the shell mean. Files: caltech_results/M35DEEP/ (two png, scripts). Peak shell r = 0.74 R_ph,IC (R_ph,IC = 9.554e11 cm), 35 cells above the inner wall (0.7061).

## A. Measured
| item | result |
|---|---|
| What it is laterally | a narrow ring at |k| = 19-27 (|k| = 20-24: 78 %, 24-28: 16 % of the power at 0.74, last dump); power-weighted |k| 23; k = 22 is 11.64 cells per wavelength, 1.12e10 cm, l ~ 400. Wedge fundamental (|k| <= 2) 0.1 % or less of the power now, shell mean |v_r| at 0.74 <= 0.021 km/s (4e-3 of the ring rms), |k| > 60 (2-cell scale) 0.2 %. So: not l = 0, not the fundamental, not grid scale (2 cells); a band of 9.5-13.5 cells. |k| = 44 (first harmonic) 1 %. |
| Direction | axis-aligned: 66-85 % of the ring power has |kx| <= 1 or |ky| <= 1 at the six times checked (3.0e5 to 7.7e5 s); an isotropic ring would give 8 %. All oblique angles between 15 and 75 deg hold 0.00 of the power. The strongest coefficients are (kx 0, ky +-22) and (kx +-21, ky 0). |
| Growth, whole history | it does not start at 6.4e5 s. Ring rms at 0.74: 2e-8 km/s at 3.0e4 s, 8e-6 (1.1e5), 2e-4 (1.85e5), 3.6e-3 (2.6e5), 0.020 (3.46e5), 0.14 (5.18e5), 0.43 (6.05e5), 0.97 (6.91e5), 2.4 (7.76e5). e-fold of the rms: 1.1e4 s (3.0e4-1.1e5), 2.4e4 (1.1e5-1.85e5), 2.5e4 (1.95e5-2.6e5), 4.9e4 (2.6e5-3.2e5), 6.8e4 (3.2e5-4.0e5), 9.4e4 (4.0e5-7.8e5, n 209, rms of the ln residual ~0.1). **The rate slows and then is steady from ~4e5 s at 9.4e4 s.** Same within 0.71-0.77 R_ph (9.3e4-1.0e5 s; 1.2-1.4e5 at 0.78-0.80). |
| KEr_int in the hst | gain of KEr_int from 6.0e5 to 7.7e5 s 3.45e39 erg; gain of 0.5 int rho v_r^2 over r < 0.80 R_ph in the dumps 3.43e39 erg. The hst rise is this motion (share 0.14 at 6.0e5, 0.37 at 6.6e5, 0.67 at 7.2e5, 0.85 at 7.7e5); the visible "onset" at 6.4e5 is when it overtakes the convection KE. |
| Radial structure | ring rms peaks at 0.737-0.742 and is 0 at the wall and at 0.80; phase jumps of pi (nodes) at r ~ 0.729 and 0.755 for (0,21) at 7.76e5, spacing 0.026-0.033 R_ph between dumps, i.e. lambda_r ~ 0.055-0.065 R_ph; the node radii change between dumps (7.42e5, 7.69e5, 7.76e5). |
| Frequency | from the coefficient phase at 1800 s cadence: +-1.64e-4 rad/s for (ky 22, kx 0), period 3.8e4 s; phase steps per dump +0.29 +- 0.32 rad (rms spread), a steady phase advance (mostly travelling in phi); for (ky 0, kx 21) the slope is -7e-5 rad/s with a step spread of 0.64 rad, so its frequency is not well determined. **Aliasing: the dumps fix omega only modulo 2 pi / 1800 s = 3.49e-3 rad/s; Nyquist 1.75e-3.** In particular 3 x 3.49e-3 + 1.64e-4 = 1.064e-2 rad/s is close to the acoustic cs kh = 1.07e-2. Two independent checks against the fast branch: (i) the hst (12 s) KEr_int, about 0.3-0.85 of which is this motion in the three windows used (6.5e5-7.0e5, 7.2e5-7.6e5, 7.5e5-7.8e5 s), has its variation at 1e-4 to 1.5e-3 rad/s (strongest 3.8e-4 and 2.5e-4 in the first window, resolution 1.3e-4; 2 omega = 3.3e-4) and nothing above 1.2e-2 rad/s in any of them (share of the rms 0.000; at 2 cs kh 0.0000); (ii) polarisation: rms v_t / rms v_r = 0.28 at 0.74 (kr/kh with the measured lambda_r = 0.2). Both point to the slow branch; the ring rms also modulates at about twice that frequency (growth figure). |
| Local frequencies at 0.74 | N 1.60e-4 rad/s (0.72: 1.34e-4, 0.76: 1.84e-4, 0.80: 2.65e-4; N^2 > 0 everywhere between 0.71 and 0.90), acoustic cutoff cs/2Hp 3.6e-4, cs kh 1.07e-2 rad/s, 2pi/N 3.9e4 s. So the measured 1.64e-4 rad/s is N to ~3 % (the band N 1.34e-4-1.84e-4 is spanned within 0.72-0.76). kh Hp = 15. |
| Stability / structure | convectively stable: nabla 0.249 < nabla_ad 0.274 at 0.74 (0.71-0.80: nabla = nabla_rad, N^2 > 0); beta 0.68; kappa 0.53, dlnk/dlnT -0.43 (0.71), -0.71 (0.74), -0.93 (0.78), dlnk/dlnrho +0.24; kappa_T + kappa_rho/0.39 (Gamma3 - 1 = 0.39) falls outward: +0.16 (0.71), +0.07 (0.72), -0.08 (0.74), -0.20 (0.76), -0.22 (0.78): it does not rise outward as the classical kappa mechanism needs (mean profile at 6.25 d, not a stability calculation). |
| Time scales at 0.74 | growth e-fold 9.4e4 s (rate 1.1e-5 /s, KE rate 2.1e-5 /s). 1/N 6.3e3 s; H_p/c_s 1.4e3; radiative diffusion of the lateral wave 1/(chi kh^2) 3.6e4 s (0.72: 5.8e4, 0.76: 2.1e4); thermal H_p^2/chi 8.1e6 s. The growth time is slower than the dynamical and the lateral-diffusion times, faster than the thermal one. |

## B. Code features (item 2)
- Switch of the binary (c16340, t 1.8967e5): growth e-fold 2.38e4 s before (1.1e5-1.85e5) and 2.49e4 after (1.95e5-2.6e5): no break. The ring was already growing before it. No other key change or restart in the run.
- merge_lat: no such key in the deck. Sponge: P_sponge is 0 at all times (hst). One x1 block (336 cells) and 2 x 2 lateral MeshBlocks of 128: block and tile edges are k = 2, the ring is k = 22. Stretching: the cell width ratio is 1.001357 per cell from 0.706 to 0.90 with no change (max change of the ratio 1e-12): no break at 0.74.
- Inner wall at 0.7061: the ring vanishes there; nothing else of the code sits at 0.74 in the outputs. tau(0.74) = 1.2e6, far above every source/closure switch that uses tau ~ 1 (vet_col_lat_cut_thin 1.0, he_tau_int 1.0). Ray rotation is off (z-rotation every 0 cycles); rays GL 7 x 18 staggered. mg levels 4: coarsest lateral coarsening 16 cells (k = 16); the ring is at 20-24. The clean_ray_min / clean_ray_cut shell radius is not in the hst or the log, so it cannot be compared.
- Seed: l = 16-64, r 0.940-0.978 R_ph, all m: not at 0.74, not at l ~ 400; the ring is not a seeded mode.

## C. Interpretation (not shown by the data)
- A slow (omega ~ N), lateral-scale ~11.6 cells, ring-shaped gravity-wave-like motion in the stable layer fits A. The axis alignment (power only where kx ~ 0 or ky ~ 0) is the one property that points to the discretisation or the ray set rather than to an isotropic physical instability; it is not proof.  Whether it is physical or numerical, or excited from above (the convection or the start), is open; the twin is the test.

## D. Twin: rsts (all hard-linked in rst_keep, nothing pruned)
| rst | model time [s] | ring rms at 0.74 |
|---|---|---|
| 00000.first | 0 | 0 (seed fields only) |
| 00001 | 8.6405e4 | 3.5e-6 km/s |
| 00002 | 1.7281e5 | 1.0e-4 |
| 00003 | 1.8968e5 (switch) | 2.2e-4 |
| 00004 | 2.5920e5 | 3.6e-3 |
| 00005 | 3.4561e5 | 0.020 |
| 00006 | 4.3201e5 | 0.065 |
| 00007 | 5.1840e5 | 0.14 |
| 00008 | 6.0481e5 | 0.43 |
| 00009 | 6.9120e5 | 0.97 |
| 00010 | 7.7761e5 (written 10:11) | 2.4 (last dumps) |
- There is no rst at 7.2e5. A start "before the onset" must be early: at 00001 the ring is already 3.5e-6 km/s and every later rst contains it at the listed size (a twin from 00008 starts with 0.43 km/s, not clean). 00007-00010 are also in rst/; 00001-00010 and 00000.first are all in rst_keep.

-- Caltech
