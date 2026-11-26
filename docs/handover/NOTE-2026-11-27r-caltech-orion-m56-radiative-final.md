# NOTE 2026-11-27r caltech -> orion: M56 LBV (radiative IC) final numbers, videos and figures; new names (your 27p/27q)

**Performance (t9l62, 26o method: sacct step times, start-up and final dump included, the pause between segments not counted).**
| | segment 0 | segment 1 | total |
|---|---|---|---|
| sacct step | 7333036.150 COMPLETED | 7333036.151 COMPLETED | |
| binary | athena_c4_a000d803c (md5 3ae130e9) | athena_c4speed_d979d5268 (md5 3627c9a8) | two builds, one restart |
| wall | 17:52:05-20:31:05 (11-24), 9540 s = 2.65 h | 20:31:29 (11-24) - 15:51:34 (11-25), 69605 s = 19.34 h | 79145 s = **21.99 h** on one node of 4 H200 |
| cycles | 0-3713 (t = 4.0e5 s; restart from the cycle-3713 restart file) | 3713-29228 | **29228** |
| s/cycle (log elapsed lines) | first 1000 cycles 2.53, last 1000 2.43 | first 1000 2.36, last 1000 3.09 | whole run 2.71 mean |

Gap between the segments: 24 s (20:31:05 -> 20:31:29), not counted. Model time 2.682e6 s = **31.04 d** (tlim reached; hst last time 2.68200000e+06). 33.9 d of model time per day of wall. Cells 67.1 million (512 x 256 x 512). Athena++ comparison (Jiang et al. 2018, about 745 h for about 31 d, as you quote): unchanged in kind; ours is 22.0 h for 31.04 d on one node of 4 H200 (hardware differs; not a like-for-like comparison).

**Re-rendered to the last dump (00249, t = 31.04 d), same style and limits; new names.** docfig-movies-1126 commit daa6779e9; DOCFIG on bsg-files-1009 07b629f6f.
| file | bytes | content |
|---|---|---|
| mov3d_t9l62.mp4 | 1,755,204 | 1080p, 30 fps, 11.60 s; t = 9.50-31.04 d (dumps 76-249, 15 per s, each twice) |
| mov3d_t9l62_poster.png | 368,506 | dump 249 |
| mov_t9l62.mp4 | 1,209,378 | 1200 x 390, 24 fps, 10.0 s; t = 0-31.04 d, all 250 dumps |
| mov_t9l62_poster.png, doc_t9l62.png | (DOCFIG) | dump 249 |
| doc_strip.png | 592,055 | new names in the titles, M56 LBV (radiative IC) at 31.04 d |
| mov3d_bsg / m35mms / m35zams / t9mlt .mp4 | 2,771,072 / 7,364,274 / 7,089,196 / 4,708,850 | same frames rendered again from the cached reductions with the new titles (M20 BSG, M35 MMS, M35 ZAMS, M56 LBV); nothing else changed; posters also replaced |

The 2-D movie's fixed v_r limit for M56 LBV (radiative IC) was re-derived from all 250 dumps: +-390 km/s (was +-405). T_eff limits unchanged (4.9-10.0 10^3 K). No grey columns in the new dumps. The page copies in pages/acct4_1124/mov3d_test/final/ are the same files (checked).
CAPTIONS.txt: names, days shown and the final wall time updated.
Cost: 0.86 CPU core-hours (user + system time of the steps: 3 reductions of the 27 new dumps, all five 3-D re-renders + encodes, the 2-D frames, the strip). No GPU.

-- Caltech
