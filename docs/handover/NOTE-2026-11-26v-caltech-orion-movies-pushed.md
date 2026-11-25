# NOTE 2026-11-26v caltech -> orion: the five movies are pushed; the bsg T_eff ripples

**Pushed:** branch docfig-movies-1126, one orphan commit a6c654375, only the five MP4 (mov_bsg, mov_m35mms, mov_m35zams, mov_t9l62, mov_t9mlt; 12.8 MB). Nothing on bsg-files-1009. A sixth file, mov3d_test_bsg.mp4 (26w), will follow on that branch before it can be deleted.

**The ripples, in one sentence:** they are not the display and not only the tau_R = 1 reading, because the same stripes are in the native-cell lab F_r at fixed radii well above the photosphere (55 and 60 R_sun; 4 % high-pass rms) and are absent at 48 R_sun below it (1.4 %); whether they come from the angular (ray) set or from the photosphere's surface crossing the radial cells I cannot tell from these data.

What was measured (bsg dump 00114, float64 bins, region phi 26-36 deg, latitude 4 to 11 deg unless stated; DOCFIG/ripple/rip.png):
- Comoving F_r at tau = 2/3, at tau = 1 and in the own tau = 1 cell (the maps): the same stripes (high-pass correlation 0.99 and 0.93 with the first). The page and the movie read the same cells, so they are in both.
- The stripes are also in log rho read at tau = 2/3, but not visible in the radius r(tau = 2/3) or r(tau = 1) itself (0.1 % rms, smooth); r falls by 1.1 R_sun over the region, 7.0 radial cells.
- Weak link to the position of the tau = 2/3 surface inside its radial cell: that position (one harmonic) explains 13 % of the variance of the lab F_r at 55 and 60 R_sun and 21 % of the own-cell reading (the same fit with the position shifted by 40 cells: 0 %). In dumps 113, 110 and 104 the same test gives 0-1 %: the link does not generalise.
- Same place over dumps: no. The high-pass pattern of dump 113 (12 h earlier), 110 and 104 correlates with dump 114's at 0.20, -0.01 and 0.18 in this region: not a pattern fixed to the grid or to fixed ray directions in (theta, phi), but the flow changes too in 12 h, so this does not exclude a ray effect.
- Not tested: a different ray set (GL 7 x 18 is the only one in these data), a static or constructed source.

-- Caltech
