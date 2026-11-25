# NOTE 2026-11-26o orion -> Caltech: star figures for the documentation, with wall time on one GPU node (the USER's request, 11:56 PST)

(26o is not a cancel of 26m; "26n" stays reserved for that.)

The USER wants your star plots in the documentation site, emphasizing the wall time on ONE GPU node, with a prettier colormap: the plots in the pages you published use a colormap the USER finds not that pretty. Please (post-processing only; no new star runs, one CPU job per star, nothing on the production nodes; m35mms only from dumps / after it has finished):

1. **Palette sheet first** (one image per field, the same snapshot of one star, 4 candidates each, labelled): temperature / density: cmasher ember, sunburst, cosmic, or cmocean thermal; v_r and vorticity-like signed fields: cmocean balance, cmasher prinsenvlag or fusion, crameri vik; radiation flux / luminosity maps: crameri lajolla or batlowW. Perceptually uniform maps, no jet / rainbow, no plain viridis; limits symmetric about 0 for signed fields. Put the sheet in your branch as PNG under docs/handover/caltech_results/DOCFIG/ and, if you publish pages, give the USER the link. The USER picks from it; do not render the final set before.
2. **Table of wall times** for bsg, m35mms, m35zams, t9l62, t9mlt on one node of 4 H200: cells, model time reached, cycles, wall time (running time, restarts counted separately: hours, node-hours), seconds per cycle early / late, model time per day of wall, binary. Say which runs are complete and which still run. Plain numbers only, each with its source log.
3. **After the pick: final figures** per star, in the plot style of stag2/plot_style (mesa2.mplstyle + idort_style.py; no titles in the picture where the caption says it): one lateral map of v_r or T at the photosphere and one radial-lateral slice (log density or temperature), 2 to 3 panels per star, colorbar labelled in cgs or solar units, time in d, plus one panel strip for the five stars side by side for the docs landing page (wide, readable at 1200 px, dark-free background), PNG <= 600 KB each, and the plotting scripts next to them. For m35mms, t9l62, t9mlt: do NOT show the deep stable layer's velocity (the lateral ring is unexplained, 26i / 26j); show surface and envelope only, and say so in the caption text.
4. The caption text per figure in plain words with the wall time in it ("N cells, X d of model time, Y h on one node of 4 H200"), no claims on convergence or resolution.

No USER approval is needed for rendering; the docs site is pushed from Orion only after the USER has seen it.

-- orion (coordinator)
