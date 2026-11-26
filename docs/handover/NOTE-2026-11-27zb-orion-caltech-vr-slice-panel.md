# NOTE 2026-11-27zb orion -> Caltech: right panel of the star figures and 2-D movies = v_r slice (USER)

USER: in the 2-D plots and the 2-D movies, the last (right) panel should show v_r instead of delta rho/<rho>.

Please, for all five stars (bsg, m35mms, m35zams, t9l62, t9mlt):
1. Figures doc_<star>.png and doc_strip.png if it has the density row (it does not: v_r / T_eff only, then leave it): right panel = v_r [km/s] on the SAME radial-lateral slice as now (same plane, same radial range, same axes in R_sun), cell-centre values, no smoothing. Colour map: cmasher fusion as the left v_r map, symmetric limits about 0 chosen per star from the slice (e.g. 99th percentile of |v_r|), own colour bar labelled "v_r on the slice [km/s]". Floor gas (density <= 1.5 x floor) stays grey.
2. 2-D movies mov_<star>.mp4 (+ posters): same change; limits fixed over the whole movie (99th percentile of |v_r| over all dumps), everything else as now (dumps, frame rates, cross-fade).
3. CAPTIONS.txt: replace the density description of the right panel by the v_r one, with the limits used.
4. Plot style unchanged otherwise (mesa2.mplstyle + idort_style.py).
Push to caltech_results/DOCFIG (and docfig-movies-1126 for the mp4s) with a NOTE caltech -> orion; Orion puts them into the docs.
3-D movies unchanged.
