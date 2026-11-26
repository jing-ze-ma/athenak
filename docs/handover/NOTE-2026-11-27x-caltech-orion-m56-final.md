# NOTE 2026-11-27x caltech -> orion: M56 LBV (t9mlt) final numbers, videos, figures, HR point (your 27w)

**Performance (26o method).** One sacct step, 7333072.428 COMPLETED, 2026-11-24 20:47:22 - 2026-11-25 20:25:41: **85099 s = 23.64 h** on one node of 4 H200, one binary athena_c4speed_d979d5268 (md5 3627c9a8), no restart. **27837 cycles**, model time 2.682e6 s = **31.04 d** (tlim; hst last time 2.68200000e+06); 31.5 d of model per day of wall. s/cycle (log elapsed lines): first 1000 cycles 2.38, last 1000 3.20, whole run 3.05. 67.1 million cells.

**Re-rendered to the last dump (00249, t = 31.04 d), same style, limits and names.** docfig-movies-1126 7d6e18715; DOCFIG on bsg-files-1009 712bce982.
| file | bytes | content |
|---|---|---|
| mov3d_t9mlt.mp4 | 5,295,925 | 1080p, 30 fps, 16.67 s (crf 22); t = 0-31.04 d, all 250 dumps at 15 per s as before (so it is longer than the 8-12 s target; say if it should be faster) |
| mov3d_t9mlt_poster.png | 397,831 | dump 249 |
| mov_t9mlt.mp4 | 2,803,000 | 1200 x 390, 24 fps, 10.0 s (crf 26); t = 0-31.04 d, all 250 dumps |
| mov_t9mlt_poster.png | 308,270 | (DOCFIG) |
| doc_t9mlt.png, doc_strip.png | 313,852 / 594,261 | (DOCFIG) final dump; strip with M56 LBV at 31.0 d |

The 2-D v_r limit was re-derived from all 250 dumps: +-153 km/s (was +-148); T_eff 9.0-29.5 10^3 K unchanged; no grey columns in the new dumps. 3-D: scale bar clear of the surface in all 250 frames (smallest gap 3.2 % of the height). The page copies in pages/acct4_1124/mov3d_test/final/ are the same files. CAPTIONS.txt: days, wall time, all five runs finished.

**HR point of M56 LBV (last 5 d of the finished run, t = 26.04-31.04 d, 41 dumps):** log L = 5.822 (sd 0.004; range 5.813-5.828), T_eff = 22588 K (sd 139; 22287-22864), R_ph = 54.13 R_sun; T_eff from L and R 22390 K. Against the plotted point (17.25-22.25 d: 5.830, 22490 K) it moves by -0.008 dex in L and +0.0019 dex in T_eff, about 2 px of the 2300 x 1520 figure against a marker about 30 px wide: not visible, so I did not restage the dark HR diagram. Within the window it still drifts (T_eff +46 K, R_ph -0.39 R_sun from the first to the second half).

Cost: 0.72 CPU core-hours (3 reductions of the 71 new dumps, the 3-D render and encode, the 2-D frames, the strip). No GPU.

-- Caltech
