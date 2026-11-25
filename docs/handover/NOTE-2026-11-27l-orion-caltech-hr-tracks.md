# NOTE 2026-11-27l orion -> Caltech: HR tracks of the four runs, snapshot by snapshot (15:41 PST, the USER)

The USER: "can you draw the luminosity-temperature of our 4 runs like lines connecting the value got in each snapshot? I expect a pulsation to look like a spiral or something."
Please, from the same reduction as your 27b2 table (hrd.py; no run): for bsg, m35mms, m35zams, t9mlt (not t9l62):
- one panel per star (2 x 2), log T_eff on x (reversed, hot left), log L/L_sun on y, both ZOOMED to the star's own excursion (the spreads are 0.0001-0.007 dex), points = every dump (T_eff median over columns at tau_R = 2/3 and L_top x 4 pi/Omega at the dump time, as in 27b2), connected by a line in time order, line and points coloured by time (one sequential map, colour bar in days), first dump marked as an open circle, last as a filled star; skip the first days of the ramp/transient only if you say how many and show them in light grey.
- a second version with L from the hst at every cycle is NOT possible for T_eff, so keep the dump cadence; say the cadence (bsg 12 h, m35 30 min, t9mlt 3 h) in the caption, and whether the cadence resolves the expected pulsation periods (if a period shows, give it).
- third panel set: the same four tracks on the dark HR background (viper 27i, hrd_3d_dark_1_background), full scale, small track lines, to show where they sit.
- Plot style: mesa2.mplstyle + idort_style (stag2/plot_style; ask if you do not have it) for the zoomed panels; the dark one in viper's style.
Push png + pdf as hr_tracks_{zoom,dark}.{png,pdf} to docfig-movies-1126, with a NOTE: what loops / spirals are seen, their periods if any, and the trends (e.g. t9mlt still cooling).

-- orion (coordinator)
