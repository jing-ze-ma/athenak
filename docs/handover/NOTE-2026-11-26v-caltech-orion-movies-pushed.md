# NOTE 2026-11-26v caltech -> orion: the five movies are pushed; the bsg T_eff ripples

**Pushed:** branch docfig-movies-1126, one orphan commit a6c654375, only the five MP4 (mov_bsg, mov_m35mms, mov_m35zams, mov_t9l62, mov_t9mlt; 12.8 MB). Nothing on bsg-files-1009. A sixth file, mov3d_test_bsg.mp4 (26w), will follow on that branch before it can be deleted.

**The ripples, in one sentence:** they are not the display and not only the tau_R = 1 reading, because stripes of the same kind are in the native-cell lab F_r at fixed radii above the photosphere (55 and 60 R_sun); whether they come from the angular (ray) set or from the radial cells the photosphere crosses I cannot tell from these data.

What was measured (bsg dump 00114, float64 bins, high-pass = value minus a 9 x 9 cell mean; region A = phi 26-36 deg, latitude 1 to 11 deg; region B = latitude 4 to 11 deg; DOCFIG/ripple/rip.png):
- Comoving F_r at tau = 2/3, at tau = 1 and in the own tau = 1 cell (the maps): the same stripes (high-pass correlation 0.99 and 0.93 with the first, region A).
- Lab F_r at 55 and 60 R_sun: stripes of the same look (contours around the same centre), 4 % high-pass rms, but their cell-by-cell correlation with the stripes at tau = 2/3 is only 0.25 and 0.23 (region A), so they are the same kind of pattern, not the same lines. At 48 R_sun (below the photosphere) the high-pass rms is 1.4 % (0.6 % in region B), correlation 0.22. The page and the movie read the same cells, so they are in both.
- The stripes are also in log rho read at tau = 2/3, but not visible in the radius r(tau = 2/3) or r(tau = 1) itself (0.1 % high-pass rms, region A; smooth in the picture); r(tau = 2/3) falls by 1.1 R_sun over region B, 7.0 radial cells.
- Weak link (region B) to the position of the tau = 2/3 surface inside its radial cell: that position (one harmonic) explains 13 % of the variance of the lab F_r at 55 and 60 R_sun and 21 % of the own-cell reading (the same fit with the position shifted by 40 cells: 0 %). In dumps 113, 110 and 104 the same test gives 0-1 %: the link does not generalise.
- Same place over dumps: no. The high-pass pattern of F_r at tau = 2/3 of dump 113 (12 h earlier), 110 and 104 correlates with dump 114's at 0.20, -0.01 and 0.18 (region A): not a pattern fixed to the grid or to fixed ray directions in (theta, phi), but the flow changes too in 12 h, so this does not exclude a ray effect.
- Not tested: a different ray set (GL 7 x 18 is the only one in these data), a static or constructed source.

-- Caltech
