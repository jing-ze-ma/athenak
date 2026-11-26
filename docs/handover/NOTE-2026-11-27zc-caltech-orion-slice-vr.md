# NOTE 2026-11-27zc caltech -> orion: right panel = v_r on the slice (your 27zb)

Done for all five stars: doc_<star>.png, mov_<star>.mp4 and posters; the strip (no density row) unchanged. Right panel: v_r on the same radial-lateral slice as before (cell values from the reductions, no smoothing), cmasher fusion, symmetric limits, colour bar "v_r on the slice [km/s]", floor gas grey (same mask as before), same radial ranges (M35 MMS and both M56 LBV without the deep stable layer). Other panels unchanged.
Limits: figure = 99th percentile of |v_r| over the slice cells drawn in that dump; movie = 99th percentile pooled over about 40 dumps across the run (and the last), fixed:
| star | movie +- [km/s] | mov_<star>.mp4 bytes | doc_<star>.png bytes |
|---|---|---|---|
| M20 BSG | 61 | 2,145,051 | 321,242 |
| M35 MMS | 137 | 5,714,563 (40 s, interpolated, crf 28) | 364,200 |
| M35 ZAMS | 71 | 5,073,806 (40 s, interpolated, crf 24) | 163,023 |
| M56 LBV (radiative IC) | 395 | 834,319 | 290,105 |
| M56 LBV | 206 | 2,480,989 | 322,139 |

Look: in M35 ZAMS and M56 LBV (radiative IC) the largest |v_r| sit at the top of the slice (the photosphere layers and the outflow above them), so the limits are set there and the interior flow shows in pale colours. That is the rule as asked (99th percentile); a limit from the convection zone alone would show the interior better; say if you want it.
Pushed: docfig-movies-1126 973d3e90c (the five mov_*.mp4); DOCFIG on bsg-files-1009 4be7848c5 (five doc_*.png, five mov_*_poster.png, docfig.py, CAPTIONS.txt). Copies in pages/acct4_1124/mov2d_26z/. Cost 0.26 CPU core-hours (02-09, 02-11, 02-17).

-- Caltech
