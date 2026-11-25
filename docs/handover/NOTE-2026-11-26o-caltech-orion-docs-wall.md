# NOTE 2026-11-26o caltech -> orion: palette sheets and the wall-time table (your 26o items 1 and 2)

## 1. Palette sheets
PNG in caltech_results/DOCFIG/: palette_temperature.png, palette_density.png, palette_vr.png, palette_flux.png (1200 px wide, bsg final dump 00114, t = 57.0 d; README.txt there lists the candidates). The flux sheet has two maps added to the two you named (cmocean solar, cmasher amber), labelled "added". A page with the same four sheets is staged for the USER's pick. No final figure, movie or page re-render was made.

## 2. Wall time on one node of 4 H200 (read 2026-11-25 12:04:41 PST; m35mms and m35zams rows updated at 12:4x PST after both finished)
One node, 4 GPUs, 4 ranks. Wall per segment = the Slurm step time (sacct; includes the init and the final restart dump); the gap between two segments is listed apart. Each segment is counted once. Cycles and model time of the running ones are the last hst row; node-hours = hours (one node).
| star (status) | cells | segments: binary, cycles, model time reached, step wall | total step wall (restart gap not included) | gap between segments | node-hours | s/cycle early / late (first / last 1000 cycles of the segment) | model time per day of wall |
|---|---|---|---|---|---|---|---|
| bsg (**complete**, tlim 4.9248e6 s = 57.0 d) | 16,777,216 (256^3) | one: athena_sdorder3_95813cf85 (md5 30d01765), c0-29657, 4.9248e6 s, 21,699 s = 6.027 h | 6.027 h | none | 6.027 | 0.761 / 0.743 (c50-1050 / c28650-29650; mean 0.729) | 227 d/d |
| m35mms (**complete**, tlim 8.64e5 s = 10 d; ended 2026-11-25 12:13:57) | 22,020,096 (336 x 256 x 256) | (a) sdorder3_95813cf85 (30d01765), c0-16340, 1.8967e5 s, 15,477 s = 4.299 h; (b) from c16340, athena_c4speed_d979d5268 (3627c9a8), to c73280, 8.64e5 s, 56,112 s = 15.587 h | 19.886 h | 30 s (20:38:15-20:38:45) | 19.89 | (a) 0.660 / 0.974 (mean 0.942); (b) 0.950 / 0.962 (mean 0.984) | 12.1 d/d |
| m35zams (**complete**, tlim 8.64e5 s = 10 d; ended 2026-11-25 12:26:28) | 6,291,456 (384 x 128 x 128) | (a) sdorder3_95813cf85, c0-50816, 1.7687e5 s, 16,047 s = 4.457 h; (b) from c50816, c4speed_d979d5268, to c248240, 8.64e5 s, 56,225 s = 15.618 h | 20.076 h | 29 s (20:48:54-20:49:23) | 20.08 | (a) 0.289 / 0.296 (mean 0.315); (b) 0.287 / 0.296 (mean 0.285) | 12.0 d/d |
| t9l62 (**running**, 83.5 % of tlim 2.682e6 s = 31 d) | 67,108,864 (512 x 256 x 512) | (a) athena_c4_a000d803c (3ae130e9), c0-3713, 3.5958e5 s, 9,540 s = 2.650 h; (b) from c3713, c4speed_d979d5268, running: to c24587, 2.2389e6 s, 55,995 s = 15.554 h | 18.204 h | 24 s (20:31:05-20:31:29) | 18.20 | (a) 2.529 / 2.428 (mean 2.510); (b) 2.358 / 2.728 (mean 2.681) | 34.2 d/d |
| t9mlt (**running**, 66.2 % of tlim 2.682e6 s = 31 d) | 67,108,864 (512 x 256 x 512) | one so far: athena_c4speed_d979d5268 (3627c9a8), c0-18390, 1.7750e6 s, 55,043 s = 15.290 h | 15.290 h | none | 15.29 | 2.383 / 3.088 (mean 2.986) | 32.2 d/d |

- Sources: step times sacct 7333035.52 (bsg), 7333073.143 + .145 (m35mms), .144 + .146 (m35zams), 7333036.150 + .151 (t9l62), 7333072.428 (t9mlt); s/cycle from the `elapsed=` lines of prod_c1_1124/bsg/log.0, prod_c1_1124/m35mms/log.0 and log.1, prod_c1_1124/m35zams/log.0 and log.1, prod_c4_1124/t9l62/log.0 and log.1, prod_mlt_1125/t9mlt/log.0; cycles, model times and binary md5 from each run's user.hst, chunks.txt. Script caltech_results/DOCFIG/wall.py.
- Binary change: m35mms at c16340, m35zams at c50816, t9l62 at c3713 (a restart, same deck); bsg and t9mlt have one binary. The hst row count is one or two larger than the cycle number of the log (first and last rows), so I quote the log's cycle number for complete segments.

## 3. The USER's reply on the palette sheet (2026-11-25 ~12:15) and round 2
Verbatim: "i dont see much contrast in the density and temperature for the conv pattern. cmasher fusion for vr. flux lets try to use Blues since these are blue stars". So: v_r = cmasher fusion (settled); flux = Blues to try; density and temperature = the fluctuation from the shell mean is wanted.
Round 2 (bsg final dump 00114, t = 57.0 d; caltech_results/DOCFIG/r2_*.png, README.txt, r2_limits.json):
- Density and temperature: (a) delta rho / <rho> and delta T / <T> (shell mean) on the theta-r slice at phi = 30 deg (34-54 Rsun) in cmasher fusion, cmocean balance, crameri vik; (b) log10 rho and log10 T zoomed on 37-46 Rsun in cmocean thermal and cmasher ember; side by side. Limits symmetric for (a), from the 99th percentile of |fluctuation| of the convection-zone cells 38.5-44.0 Rsun only (slice: +-0.385 rho, +-0.0885 T); (b) 1-99 percentile of the same cells (log rho -7.71 to -7.09, log T 5.057 to 5.310). Floor gas (density <= 1.5 x floor) grey.
- Lateral maps on 40.70 Rsun (the v_r shell; it also has the largest rms v_r, 18.45 km/s, of the shells between 36 and 46 Rsun) and 43.52 Rsun (near the top of the convection zone), same three maps; limits +-0.282 / +-0.240 (rho), +-0.0682 / +-0.0487 (T).
- Flux (F_r / mean at tau_R = 1): Blues (low light, high dark), Blues_r (low dark, high light), cmocean dense (low light, high dark), cmocean ice (low dark, high light), all with the same limits 0.473-1.500 (2-98 percentile). Found in the earlier work: the 3-D photosphere coloured by T_eff used Blues (high = dark; vis3d/cs3d.py), the 2-D flux map of the bsg page used Blues_r (high = light).
- v_r in cmasher fusion (limits +-64.8 km/s) as the reference.
- Staged for the USER's second pick in the same page folder (Palette Sheet). No final figure, movie or page re-render.

-- Caltech
