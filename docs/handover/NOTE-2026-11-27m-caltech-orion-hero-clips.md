# NOTE 2026-11-27m caltech -> orion: hero clips (27f/27g/27m)

Both on docfig-movies-1126 (latest 4b5535dc8); scripts in caltech_results/DOCFIG/hero/. H.264 yuv420p, faststart, no audio, 960 x 540, both decode with 0 errors.

| file | bytes | length | content |
|---|---|---|---|
| hero_bsg.mp4 | 603,043 | 6.88 s, 24 fps (8 dumps per s, each shown 3 times) | bsg 3-D frames rendered again from the cached reductions without title and colour bars (camera as the movie, wedge centred), small time stamp, 10 R_sun bar; dumps 51-101 (t = 25.5-50.5 d) plain, then 102-105 cross-faded into 47-50, the time stamp jumps from 52.5 d back to 25.5 d at the loop |
| hero_bsg_poster.png | 151,466 | | dump 101 |
| hero_rsg.mp4 | 1,442,864 | 7.74 s, 23 fps (the source rate) | cut from viper's mov3d_rsg.mp4 (d261d7cdc), frames 115-292 of 420, cropped to the star's box over these frames (title and colour bars removed, no time stamp or scale left) and padded black to 16:9; last 6 frames cross-faded into the first, crf 28 |
| hero_rsg_poster.png | 302638 | | frame 286 |

Loop seams: chosen where the picture changes least (rms of the image difference 19.8 for bsg against 18.6 between consecutive dumps; rsg 35.9 against 24.6 between consecutive frames, so its loop point still shows a small jump softened by the cross-fade). In bsg the photosphere radius at the two ends differs a little; I looked at one cross-faded frame and saw no clear double edge.

-- Caltech
