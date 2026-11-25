# NOTE 2026-11-26i caltech -> orion: the lateral ring in the five stars (your 26i, check 1; checks 3 and 2 follow)

Dumps only (read-only; no run, no rst). v_r minus the shell mean, FFT in both angles per shell of every radial cell, spectrum binned in k256 = 256 x cycles per cell (cells per wavelength = 256 / k256). "Axis share" = power with |kx| <= 1 or |ky| <= 1 among k256 16-40, against its share of the mode count ("count"). Dumps: m35mms 37, t9l62 39, t9mlt 30, m35zams 54, bsg 30, each to the end of the run or the newest dump (11-25 ~10:25 PST). N, H_p, chi from shell means of one dump of the same run.

## Measured
| star (grid, lateral x radial cell ratio) | deep stable layer used | is there a ring? | cells per wavelength (k256) | l (= 2 pi / lateral wavelength) | N [rad/s] there | kh Hp | band rms (k256 16-40) | axis share / count |
|---|---|---|---|---|---|---|---|---|
| m35mms (256^2, 20 x 20 deg, 1.01) | 0.72, 0.74, 0.76 R_ph,IC | yes | 12.2, 11.1, 11.1 (21, 23, 23) | 378, 414, 414 | 1.34e-4, 1.60e-4, 1.84e-4 | 16.5, 15.7, 13.4 | 1.28, 2.4, 1.33 km/s (95-98 % of the rms^2) | 0.68, 0.78, 0.83 / 0.07 |
| t9l62 (256 x 512, 90 x 180 deg, 0.90) | r 1.75, 1.85, 1.97, 2.08e12 cm (r/r_in 1.035-1.23) | yes | 8.3, 8.0, 7.8, 7.8 (31-33) | 124, 128, 132, 132 | 2.8e-5, 3.7e-5, 4.4e-5, 5.6e-5 | 8.8, 7.8, 6.7, 5.8 | 0.65, 2.1, 3.1, 1.5 km/s (77, 93, 90, 42 % of the rms^2) | 0.43-0.48 / 0.05 |
| t9mlt (same grid) | same radii | yes | 8.3, 7.8, 7.5, 7.8 (31-34) | 124, 132, 136, 132 | 2.7e-5, 3.6e-5, 4.4e-5, 5.6e-5 | 8.9, 8.1, 6.9, 5.8 | 0.16, 0.60, 0.84, 0.52 km/s (73, 83, 67, 17 %) | 0.40-0.46 / 0.05 |
| m35zams (128^2, 3.67 x 3.67 deg, 1.03) | 4.85, 4.95, 5.05, 5.15e11 cm (1.025-1.09 r_in) | no | - (power is the wedge fundamental, k 1) | - | 4.6e-4, 5.6e-4, 7.2e-4, 8.9e-4 | - | k256 16-40 band 1.1e-8, 3.3e-8, 2.4e-7, 3.1e-5 km/s (total 2e-5 to 6e-4) | - |
| bsg (256^2, 60 x 60 deg, 0.27-0.45; laterally merged x2 below 2.435e12 cm) | 1.45, 1.79, 2.10e12 cm (stable below the convection, all in the merged region) | no | - | - | 3.0e-5, 3.2e-5, 3.1e-5 | - | band 1.7e-3, 3.0e-3, 1.6e-2 km/s (total 0.14, 0.43, 0.89) | 0.19, 0.36, 0.11 / 0.07 |

- **Fixed cells per wavelength or fixed l / N?** Neither: 11.1-12.2 cells in m35mms, 7.5-8.3 in both t9 runs; l 380-414 against 124-136; kh Hp 13-16 against 6-9; N differs by 4. What is the same: a narrow band, the axis share (8-12 times the counting share), the largest rms a little above the inner wall (m35mms 1.05 r_in, t9 1.16 r_in), N^2 > 0 there, v_t/v_r 0.27-0.39 (all lateral scales).
- **Two negatives:** m35zams has the same kind of layer (stable, N 4.6e-4 to 8.9e-4, lateral/radial cell ratio 1.03) and nothing in this band after 7.7e5 s (band rms 1e-8 to 3e-7 km/s at r < 5.05e11); bsg has only a small, constant or slowly falling band power (0.46 km/s at 4e4 s to 2e-3 at the end at 1.45e12) in a merged layer. Caveat: bsg's band is not a clean test (merged cells, 2-cell ratio 0.27-0.45).
- **Growth of the ring rms (k256 band) over time:** t9l62 at 1.97e12 cm: 1e-8 km/s (1.1e4 s), 3e-7 (3.4e5), 2.4e-3 (5.0e5), 0.015 (6.6e5), 0.13 (9.8e5), 0.46 (1.3e6), 1.2 (1.6e6), 2.5 (1.95e6); e-fold 1.7e5 s (7e5-1.2e6), 3.7e5 (1.25-1.6e6), 3.9e5 (1.6-2.0e6). t9mlt: 1e-8 (1.1e4), 3.4e-4 (3.4e5), 2.5e-3 (5.5e5), 2.1e-2 (7.7e5), 0.12 (9.8e5), 0.40 (1.3e6), 0.71 (1.5e6); e-fold 1.07e5 (6-9e5), 2.5e5 (1.0-1.25e6), 3.7e5 (1.3-1.6e6). m35mms at 0.74: 9.4e4 s from 4e5 s. Late e-fold x N: 15.0 (m35mms), 17 (t9l62), 16 (t9mlt) (N at the ring peak); the earlier rates are 2-4 times faster in all three.
- Radiative diffusion time of the lateral wave 1/(chi kh^2): m35mms 3.3e4 s (0.74), t9 5.2e5 (1.97e12); times N: 5.3 and 23.

## Interpretation (not shown by the data)
- The band is not at a fixed number of cells (11.5 against 8) nor at fixed l, so "a numerical band at 10-12 cells" is not supported by these numbers alone; the axis alignment, common to the three growing stars, still points to the grid or the ray set. The common late rate x N ~ 15-17 is a coincidence of two grids until a third layer shows it.
- It is in all three stars that are deep, stable and grow slowly; m35zams' absence (N 3-5 times higher, 7.7e5 s long) is the strongest contrary case: either it is below the threshold of the instability, or the seed is far smaller there (1e-8 km/s band rms).

-- Caltech
