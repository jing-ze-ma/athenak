# NOTE 2026-11-26w caltech -> orion: 3-D videos: how it was done, plan, cost, test

Test: mov3d_test_bsg.mp4 is on docfig-movies-1126 (second commit 263161e3e; the five movies are untouched). No other render is made before the USER's go.

## 1. How it was done before
Directory vis3d (Caltech, /resnick/groups/carnegie_poc/jingze/vis3d). CPU only, software rendering (Mesa llvmpipe / OSMesa), no GPU.

| | wedges: render3d.py (RENDER_VERSION v7) | global cubed-sphere runs: cs3d.py (c3), README_cs3d.txt |
|---|---|---|
| runs | bsg_li, hegiant, lbv_e1s (10-14 to 10-17) | rsgE_glob, hggl_glob (last built 11-06), rsgrc3 (rsg_rc3_page/vis3d) |
| engine | PyVista 0.49.1 + VTK 9.7.1, EGL | VTK 9.3.1 OSMesa |
| photosphere | exact tau_R = 2/3 per column (tau on the radial faces from the top, kappa_R rho constant per cell, run's own table), true scale, geometry 2 light 3x3 passes, colour T_eff = (lab F_r(r_ph)/sigma)^(1/4), F_r linear between cell centroids, no smoothing of T_eff, corner-shaded | same, true scale always (relief = 1) |
| inside | quarter-wedge notch (phi > phi_mid, theta < 90 deg) + outer side faces, v_r out to the local photosphere; inner boundary grey; no volume rendering, no density | 90 x 90 deg wedge cut from the sphere, 3 faces v_r; "layer" view: nested spheres |
| colours | T_eff Blues, v_r coolwarm | Blues (He giant) / hot (RSG), v_r coolwarm |
| limits | fixed per movie = 1-99 % of the newest dump; v_r symmetric | same (a new newest dump re-renders) |
| camera | fixed: azimuth 30 deg, elevation 30 deg from the wedge axis, view angle 30 deg; no rotation; not looping | fixed, looking into the cut |
| frames, size | 1000 x 760, one frame per dump (at most 300) | 2000 x 1520 |
| encoding | H.264 yuv420p faststart, crf 22-34, 1.3 x fps rule | same, crf raised below the cap |

Files on disk (movie.json, ffprobe): bsg_li.mp4 116 frames, 10 fps, 11.6 s, 1.5 MB; hegiant.mp4 250, 22 fps, 11.4 s, 0.75 MB; lbv_e1s.mp4 297, 26 fps, 11.4 s, 2.5 MB; hggl_3d.mp4 289, 25 fps, 11.6 s, 5.6 MB; rsgE_3d.mp4 288, 25 fps, 11.5 s, 10.3 MB (layers view 13.0 MB). So yes: what you remember (cs3d.py, Blues with high = dark) is the global He giant; the wedge version is render3d.py.
Cost from the job accounting (reduction + rendering + encoding, TotalCPU, 8 cores): bsg_li 79 frames 5.4 CPU-min (4.1 s per frame) and 117 frames 9.7 CPU-min (5.0 s); hegiant 101 frames 5.9 CPU-min (3.5 s), 189 frames 44.1 CPU-min (14 s); lbv_e1s 275 frames 51.3 CPU-min (11 s). No GPU hours. For the global runs the README and the logs give about 10 s (He giant) to 50 s per dump for the reduction and 1.5-3 s per frame at 1000 x 760; no accounting record of the total was checked.

## 2. Plan (bsg, m35mms, m35zams, t9l62, t9mlt)
Style as render3d.py v7 except the colours: T_eff cmocean ice, v_r cmasher fusion (as on the pages), black background, true scale. New: 1920 x 1080. Code: caltech_results/DOCFIG/mov3d/w3d.py (imports the geometry functions of render3d.py).
- The old notch shows two v_r faces (one vertical, one horizontal), plus the two outer side faces (seen edge-on at this camera); the test keeps it because it is the style the USER chose. A variant with one full side face and no notch is in the code (--variant face), not rendered.
- No density face is shown, so balance is not used.
- Camera: fixed, azimuth 30 deg, elevation 30 deg, as before. The t9 wedges are 90 x 180 deg and the m35zams wedge is 3.7 deg wide, so their azimuth and zoom will need to be set per star; not tried.
- Limits fixed over the movie: 1-99 % of T_eff and 99 % of |v_r| on the faces at the newest dump (bsg: T_eff 12.3-17.4 kK, v_r +-49 km/s); t9l62 T_eff 6.5-12.5 kK as in the docs movie (early frames saturate), v_r 99 % at the newest dump.
- Deep layers: cut faces only from the base of the docs slices outward: m35mms 0.86, t9l62 0.42, t9mlt 0.62 (R_ph,IC); m35mms and t9 without the deep stable layer, the inner sphere at that radius is grey. bsg and m35zams from the inner boundary of the domain.
- 8-12 s at 24 fps: bsg 115 dumps at 12 per s = 9.6 s; m35mms and m35zams every 2nd of 481 = 241 frames = 10 s; t9l62 222 and t9mlt 178 dumps = 9.3 and 7.4 s at 24 fps (t9mlt at 20 per s = 8.9 s); t9 runs are still writing, more dumps by then.
- H.264 yuv420p faststart, crf 20-26 up to 6 MB, poster = last frame.

## 3. Cost (bsg measured, the others scaled by cells) and test
Test (bsg dumps 95-114 = t 47.5-57.0 d, 20 consecutive dumps, 12 dumps per s written at 24 fps, 40 frames, 1.7 s; mov3d_test_bsg.mp4 550,461 B, 1920 x 1080, H.264 yuv420p faststart): on hpc-sm-02-09 with 4 worker processes: 23 s wall, 38.5 s user CPU; per dump 1.3-3.5 s reduction + 0.6-1.9 s rendering. The first dump (00114) was already reduced, and R_ph, T_eff and v_r numbers of that frame are the same as in the earlier single-frame test.
Estimate for the five (1 dump = 1 frame, 4-5 core-s for bsg (16.8 M cells, float64), scaled by cells): bsg 115 frames 0.15 core-h; m35mms 241 x 22 M cells 0.4; m35zams 241 x 6.3 M 0.15; t9l62 222 x 67 M 1.0; t9mlt 178 x 67 M 0.8. Sum 2.5 core-h, within a factor 2 (the memory of one t9 reduction, about 30 GB, was not tried); no GPU hours. Sizes: 27 kB per dump at crf 20 from the test: bsg 3.2 MB, m35mms about 6, m35zams about 3, t9l62 about 6, t9mlt about 5; crf is raised if one passes 6 MB. The sizes of the other four are guesses from bsg.

-- Caltech
