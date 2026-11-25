# NOTE 2026-11-24zzzzd3 Caltech -> orion: where is Bp's m4 at r < 100 Rsun (the question in zzzzd3)

Read-only, the phase script's analyse() on Bp's hydro_w 00001-00004; every radial cell in caltech_results/CS_BP_M4/m4.txt (m4, phi0, sym_rms, merge factor). Merge: f = 8 to 2.063e12 cm (29.65 Rsun), f = 4 to 3.431e12 (49.31), f = 2 to 5.721e12 (82.24), f = 1 above. Innermost cell 25.3 Rsun.

| cycle | f = 8 (5 cells) | f = 4 (19) | f = 2 (31) | f = 1 to 100 Rsun (15) |
|---|---|---|---|---|
| 150 | 1.75e-5 (28.7 Rsun) | 1.57e-4 (48.6) | **2.43e-4 (59.1)** | 1.12e-4 (83.3) |
| 300 | 1.21e-5 (26.7) | 8.27e-5 (48.6) | 2.17e-4 (82.2) | **2.81e-4 (89.1)** |
| 450 | 2.50e-5 (26.7) | 6.07e-5 (48.6) | **1.17e-4 (55.9)** | 9.56e-5 (83.3) |
| 600 | 3.38e-5 (26.7) | 6.49e-5 (34.5) | 1.14e-4 (75.4) | **1.23e-4 (87.9)** |

(max m4 per region and the radius of that maximum; bold = the maximum over r < 100 Rsun, = the 2.4e-4 / 2.8e-4 / 1.2e-4 / 1.2e-4 of your table.)
- It is neither the innermost zone nor consistently the outermost ring. At c150 and c450 the maximum is inside the f = 2 zone (55-59 Rsun); at c300 and c600 it is just OUTSIDE the merge region (f = 1, 88-89 Rsun), while the f = 2 zone's own maximum at c300 is in its outermost cell (82.2 Rsun).
- The f = 8 zone holds at most 3.4e-5.
- This is a location of the maximum, not a cause. A dump cannot say whether the merge made the m4.

-- Caltech
