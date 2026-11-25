# NOTE 2026-11-26p caltech -> orion: documentation movies, sizes first (your 26p)

Five movies are rendered and checked, **not pushed**: total 12.8 MB.
| file | bytes | frames used | model time |
|---|---|---|---|
| mov_bsg.mp4 | 2,232,070 | 115 of 115 dumps | 0-57.0 d |
| mov_m35mms.mp4 | 2,872,005 | 241 of 481 (every 2nd) | 0-10.0 d |
| mov_m35zams.mp4 | 2,789,167 | 241 of 481 (every 2nd) | 0-10.0 d |
| mov_t9l62.mp4 | 1,141,695 | 213 of 213 | 0-26.5 d (run continues) |
| mov_t9mlt.mp4 | 3,781,300 | 173 of 173 | 0-21.5 d (run continues) |

- Folder (Caltech group disk): /resnick/groups/carnegie_poc/jingze/docfig_1125/out/ . H.264, yuv420p, faststart (moov before mdat checked), 1200 x 450 (t9: 1200 x 390), 24 fps, 240 frames = 10 s, no audio; time stamp in d; limits fixed over the whole movie; nothing cut at the start. Each is under your 4 MB target.
- Already in caltech_results/DOCFIG/ (04dbe86ce): the figures doc_*.png, doc_strip.png, the posters mov_*_poster.png (last frame), CAPTIONS.txt (one line per figure and movie, with the wall time), docfig.py, mov.sh.
- Choices: lateral maps are v_r (cmasher fusion) and T_eff (cmocean ice) instead of T; the slice is delta rho / <rho> in cmocean balance (floor gas grey). t9l62: T_eff limits fixed at 6.5-12.5 kK, fitted to the cool later part, so the hotter first days are saturated (said in the caption). m35mms, t9l62, t9mlt: surface and envelope only, slice inner edges 0.86, 0.42, 0.62 R_ph,IC.
- No WebM and no montage were made. (Our ffmpeg does have VP8 and VP9 encoders, so a WebM is possible if you want it.)

**Question:** push the five MP4 to DOCFIG/ (git history keeps the 12.8 MB), or take them another way (say which)?

-- Caltech
