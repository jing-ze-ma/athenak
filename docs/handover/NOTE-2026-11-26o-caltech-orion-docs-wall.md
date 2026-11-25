# NOTE 2026-11-26o caltech -> orion: palette sheets and the wall-time table (your 26o items 1 and 2)

## 1. Palette sheets
PNG in caltech_results/DOCFIG/: palette_temperature.png, palette_density.png, palette_vr.png, palette_flux.png (1200 px wide, bsg final dump 00114, t = 57.0 d; README.txt there lists the candidates). The flux sheet has two maps added to the two you named (cmocean solar, cmasher amber), labelled "added". A page with the same four sheets is staged for the USER's pick. No final figure, movie or page re-render was made.

## 2. Wall time on one node of 4 H200 (read 2026-11-25 12:04:41 PST)
One node, 4 GPUs, 4 ranks. Wall per segment = the Slurm step time (sacct; includes the init and the final restart dump); the gap between two segments is listed apart. Each segment is counted once. Cycles and model time of the running ones are the last hst row; node-hours = hours (one node).
| star (status) | cells | segments: binary, cycles, model time reached, step wall | total step wall (restart gap not included) | gap between segments | node-hours | s/cycle early / late (first / last 1000 cycles of the segment) | model time per day of wall |
|---|---|---|---|---|---|---|---|
| bsg (**complete**, tlim 4.9248e6 s = 57.0 d) | 16,777,216 (256^3) | one: athena_sdorder3_95813cf85 (md5 30d01765), c0-29657, 4.9248e6 s, 21,699 s = 6.027 h | 6.027 h | none | 6.027 | 0.761 / 0.743 (c50-1050 / c28650-29650; mean 0.729) | 227 d/d |
| m35mms (**running**, 99.3 % of tlim 8.64e5 s = 10 d) | 22,020,096 (336 x 256 x 256) | (a) sdorder3_95813cf85 (30d01765), c0-16340, 1.8967e5 s, 15,477 s = 4.299 h; (b) from c16340, athena_c4speed_d979d5268 (3627c9a8), running: to c72749, 8.5765e5 s, 55,557 s = 15.432 h | 19.732 h | 30 s (20:38:15-20:38:45) | 19.73 | (a) 0.660 / 0.974 (mean 0.942); (b) 0.950 / 1.000 (mean 0.984) | 12.1 d/d |
| m35zams (**running**, 98.2 % of tlim 8.64e5 s = 10 d) | 6,291,456 (384 x 128 x 128) | (a) sdorder3_95813cf85, c0-50816, 1.7687e5 s, 16,047 s = 4.457 h; (b) from c50816, c4speed_d979d5268, running: to c243795, 8.4853e5 s, 54,920 s = 15.256 h | 19.713 h | 29 s (20:48:54-20:49:23) | 19.71 | (a) 0.289 / 0.296 (mean 0.315); (b) 0.287 / 0.285 (mean 0.285) | 12.0 d/d |
| t9l62 (**running**, 83.5 % of tlim 2.682e6 s = 31 d) | 67,108,864 (512 x 256 x 512) | (a) athena_c4_a000d803c (3ae130e9), c0-3713, 3.5958e5 s, 9,540 s = 2.650 h; (b) from c3713, c4speed_d979d5268, running: to c24587, 2.2389e6 s, 55,995 s = 15.554 h | 18.204 h | 24 s (20:31:05-20:31:29) | 18.20 | (a) 2.529 / 2.428 (mean 2.510); (b) 2.358 / 2.728 (mean 2.681) | 34.2 d/d |
| t9mlt (**running**, 66.2 % of tlim 2.682e6 s = 31 d) | 67,108,864 (512 x 256 x 512) | one so far: athena_c4speed_d979d5268 (3627c9a8), c0-18390, 1.7750e6 s, 55,043 s = 15.290 h | 15.290 h | none | 15.29 | 2.383 / 3.088 (mean 2.986) | 32.2 d/d |

- Sources: step times sacct 7333035.52 (bsg), 7333073.143 + .145 (m35mms), .144 + .146 (m35zams), 7333036.150 + .151 (t9l62), 7333072.428 (t9mlt); s/cycle from the `elapsed=` lines of prod_c1_1124/bsg/log.0, prod_c1_1124/m35mms/log.0 and log.1, prod_c1_1124/m35zams/log.0 and log.1, prod_c4_1124/t9l62/log.0 and log.1, prod_mlt_1125/t9mlt/log.0; cycles, model times and binary md5 from each run's user.hst, chunks.txt. Script caltech_results/DOCFIG/wall.py.
- Binary change: m35mms at c16340, m35zams at c50816, t9l62 at c3713 (a restart, same deck); bsg and t9mlt have one binary. The hst row count is one or two larger than the cycle number of the log (first and last rows), so I quote the log's cycle number for complete segments.
- The m35mms and m35zams rows are updated once when they finish.

-- Caltech
