# NOTE 2026-11-26z orion -> Caltech: v_r at the real photosphere (the USER, 13:40 PST)

The USER asked whether the velocity plots and videos are from the photosphere (they are on a fixed sphere: bsg 40.7 R_sun in the convection zone, m35 at the initial photosphere, t9l62 at 3.8 R_ph,IC) and said yes to maps at the actual photosphere. Please, for the docs figures and the 2-D movies of ALL FIVE stars (posters too): replace the v_r panel by v_r at the photosphere of EACH COLUMN, at tau_R = 2/3 (bsg: the same tau_R as its T_eff map, or 2/3 for all: use one value for all five stars and say it), same colour map (fusion), same limits rule (fixed over the movie, symmetric), caption text "v_r at tau_R = 2/3 in each column". The T_eff and density panels stay. If the photosphere leaves the grid top or is not found in some columns, show those columns grey and say how many. Interpolation in radius linear between cell centres, no smoothing.
- Same file names (doc_*.png, mov_*.mp4, posters); the m35 ones with the slower version of your 26x (40 s); movies <= 6 MB (m35 <= 6).
- The 3-D videos (26y): the photosphere is already the surface coloured by T_eff; leave v_r on the cut faces as planned.
- Push figures to DOCFIG (new commit) and movies to docfig-movies-1126. Report the number of grey columns per star and any change of the look.
- No new run; CPU post-processing only.

-- orion (coordinator)
