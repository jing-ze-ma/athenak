---
name: sp-mhd-energy-excess
description: sp MHD dhj gains ~1% total energy in 2 rot and 3-4x horizontal KE relative to sp hydro; cs MHD does not. Field is identical at depth, 7x stronger in the sp UPPER NIGHTSIDE, temperature floors inject the energy. Cause OPEN; every sp MHD variant shows it
metadata:
  type: project
---

Measured 2026-09-06 on sp_mhd_prod (18f5dd21, from scratch) against sp_dhj_hyd and
cs_prod_mhd_rot / cs_prod_hyd_rot, all the same input (nx1 128 refit stretch, 2.8 deg).

**Facts.** E(spMHD) - E(spHYD) = +3e35 at 0.1 rot, +3.3e36 at 1 rot, +1.4e37 at 5 rot,
+7e37 at 100 rot (1 % of E); cs MHD - cs hydro stays < 1e36. sp MHD KE2/KE3 are 3-4x sp
hydro from 0.5 rot on; cs MHD/cs hydro = 0.95-1.0. Total ME is ~1e33, so the excess is
NOT field energy. Every sp MHD run shows the SAME +2.96e35 at 0.1 rot: sp_mhd_diss (old
curl), sp_noeta (no resistivity), sp_polavg, sp_polboth, sp_diss, sp_diss3, sp_dhj_ctl (old
binary, polar_emf_diss off), sp_mhd_pole2; HLLE gives 60x more (+1.7e37 at 0.1 rot).
Today's fixes (3D curl, pole, wall) change the excess by < 0.1 %.

**Where.** Event counters at cycle ~57k: eos_tfloor 3.4e7 (sp MHD) vs 8e5 (sp hydro) vs
4e4 (cs MHD); efloor 2.3e7 vs 2.5e6. At rot 2 the floor-proxy cells sit at i = 58-98, mid
latitudes, on the NIGHTSIDE (k ~ 0 = antistellar, substellar is k = 64): 453 cells vs 37
(hydro, top rows) vs 6 (cs). Deep field IDENTICAL sp vs cs (<ME> i<40: 282.7 vs 281.1
erg/cm3, max |B| 448 vs 312 G at i ~ 22); upper atmosphere (i >= 56) 0.70 vs 0.105, 96 %
of it on the nightside, 3 % at |lat| > 60. NOT the pole (no checkerboard: row-0 alternation
0.07 of |B|), NOT the periodic phi seam (hydro is equally rough there). The sp upper
nightside is violently structured in T, rho (0.04-3x hydro) with |B| ~ 5-10 G at i = 70.

**Open.** Why the sp upper nightside takes 7x the field of cs. Differences left: the theta
stretch (f = 3, only on sp), the sp-vs-cs MHD flux/EMF paths in low beta. Discriminators
not yet run: sp with f_stretch_theta = 1 (uniform, like cs); sp with bbot = 0 (must equal
hydro; tests the MHD path's hydro); cs with a theta-like refinement. See
[[sp-hydro-vs-mhd-comparison]] (the 100-rot picture), [[measure-impact-before-claiming]].
