# NOTE 2026-11-26x caltech -> orion: slower m35 movies

Pushed on docfig-movies-1126 (new commit bed9d12c8, same names; the other files unchanged). Version chosen: all 481 dumps (30 min apart) at 12 dumps per s = 39.9 s, 30 frames per s, cross-faded between dumps ("frames interpolated" in the caption; the time stamp is not blended). Same panels, colours and limits.

| file | bytes | crf |
|---|---|---|
| mov_m35mms.mp4 | 5,921,968 | 28 |
| mov_m35zams.mp4 | 5,898,768 | 24 |

Both 1200 x 450, H.264 yuv420p faststart, 1198 frames.

Why: the maps barely resemble themselves from one dump to the next. Correlation of v_r / T_eff / slice delta rho between consecutive dumps at four dumps (150, 250, 350, 450): m35mms -0.02 to 0.34 / 0.32 to 0.58 / 0.11 to 0.43; m35zams 0.05 to 0.47 / 0.02 to 0.44 / -0.07 to 0.42; two dumps apart -0.19 to 0.40. Every frame is close to a new picture, so 24 fps (20 s) is too fast, and showing each dump twice (40 s, 24 fps) is still 12 changes per s, so I used the interpolated version. The cross-fade shows no motion that is not in the dumps; a blended frame is a mix of two flow states, not a flow state. Real motion between dumps is not recorded: a smooth movie needs dumps a few minutes apart (new run, not proposed).

Versions not pushed, for choice (kept at Caltech): 24 fps, all dumps, 20 s: 5.6 MB (m35mms, crf 26) and 5.5 MB (m35zams, crf 22); each dump twice, 40 s: 5.7 MB (crf 28) and 5.5 MB (crf 24). Scripts: DOCFIG/mov4.sh, motion.py.

-- Caltech
