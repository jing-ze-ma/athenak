---
name: he-presn-sizing
description: Presupernova He stars (Woosley 2019 Table 4) sized 09-17 in bench/hestar_presn - NONE is a box problem; FeCZ spans 0.4-0.7 R, Hp/R 0.3-0.7, super-Eddington, Prad/Pgas 27-160; needs a global spherical RHD model
metadata:
  type: project
---

User idea 09-17: simulate a He star just before core collapse. Woosley 2019 (ApJ 878, 49; "He/2"
in its Table 4 = helium HALF BURNED, not half mass loss; our production star = that point of
the 3.0 Msun model, log L 3.91, Teff 81 kK, R 0.455 Rsun). Presupernova (nominal mass loss):
3.0: M 2.45, log L 4.61, Teff 9.8 kK, R 70 Rsun (expands during C burning, <3.2 Msun only);
3.5: 2.81, 4.69, 35 kK, 5.9 Rsun; 4.0: 3.15, 4.78, 49 kK, 3.4 Rsun; 5.0: 3.82, 4.91, 72 kK, 1.8 Rsun.
Columns + sizing in bench/hestar_presn (copies of column.py/sizing.py, stars.json, analyze.out,
sizing.out). Result: FeCZ base at 0.42-0.74 R (5.6 R for 3.0!), Hp/R 0.32-0.74, kappa F/(c g)
1.1-2.4 (inflated envelope), Prad/Pgas 27-160, MLT Mach 0.35-0.55, F_conv/F 0.1-0.65. The 3.0
also has an H/He recombination CZ at 1-5e4 K with Mach 0.6. Plane-parallel box invalid for ALL;
sizing.py's numbers (1e10-1e14 cells) are meaningless there. Only route: global spherical RHD
(Jiang+2015/2018 type) with the red-giant machinery, least extreme = 5.0 or 4.0 Msun.
**Gamma_e** 0.25-0.33 vs 0.044 now. See [[he-star-strange-modes]], [[rg-box-sweep]].
