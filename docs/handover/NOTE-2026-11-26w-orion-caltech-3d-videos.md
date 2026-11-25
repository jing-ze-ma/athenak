# NOTE 2026-11-26w orion -> Caltech: 3-D videos of the star runs (the USER, 13:15 PST)

The USER: "can you ask caltech to create 3d videos for those runs? ask it how we did that before". So, in this order:

1. **How was it done before?** Please answer with the tool, the script and the settings (camera, rendering, colour map, what is shown: photosphere surface coloured by T_eff, volume of density fluctuation, slices), the run it was used for, the files, the length / frames / file size, and the cost (CPU / GPU hours per frame). I remember a 3-D photosphere coloured by T_eff (vis3d/cs3d.py, Blues high = dark) mentioned in your 26o; please say if that is it, and if there were others (cs3d, wedge, RSG).
2. **Plan for bsg, m35mms, m35zams, t9l62, t9mlt** (only existing dumps, CPU post-processing on free cores; nothing on the production nodes; m35zams / bsg first): the same style as before, with the USER's colours (v_r: cmasher fusion; T_eff: cmocean ice; density fluctuation: cmocean balance), limits fixed over the movie; 8-12 s at 24-30 fps, 1080p MP4 (H.264, yuv420p, faststart), target <= 6 MB each (the docs repository must stay small); a still poster. For wedge stars (bsg, m35, t9) a wedge view with the surface and one cut face; say what the camera does (slow rotation or fixed). Nothing of the deep stable layers' velocity in m35mms, t9l62, t9mlt.
3. **Tell me the estimated render cost (core-hours) and file sizes BEFORE rendering the five**; render one test (bsg, 20 frames) first and send it for the USER to look at (same orphan branch as 26v, name mov3d_test_bsg.mp4), then wait for the USER's go.

-- orion (coordinator)
