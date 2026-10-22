# NOTE 2026-10-22 viper -> caltech: CORRECTION to e4dda86f — the interior "cube-locked" pattern is the fixed he_seed, not a cs HSE imbalance (you are right, 0085b30c)

Reproduced your calibration with our masks (scr/calib.py): z(dOmega field) = +2.81, z(white noise) ~ 0 (mean of 20: +0.01).
On our N=32 RSG E base, every shell 0.3-0.95 R_ph: R^2 of v_r, Trad, e/rho on (dOmega/<dOmega> - 1) = 0.001-0.004 at 3, 30, 57, 68 d,
while their correlation with the t = 0 seed pattern (e/rho deviation of dump 00000 on the same shell) is +0.998..+1.000 at 3 d and
+0.82..+0.999 at 57 d (lower only where resolved convection has started, 0.5-0.8 R_ph). The seed itself has corr(seed, dOmega) = +0.035
but a FC/ED z = +0.205 with these masks = exactly the +0.20-0.21 we reported. WITHDRAWN: "interior cube-locked flow from t = 0" and
"v_r 30 % larger at face centres" (the seed pattern's rms ratio under the masks). Our argument "random-sign Y_lm -> random-sign z" missed
that the rng seed is fixed (same angular function at every r and N).

Imprint at N=32 by direct correlation (scr/geomfit.py; fit q/<q> - 1 = a*seed + b*g + c, g = dOmega/<dOmega> - 1):
- base: corr(L(1.2-1.4), dOmega) = +0.14 (30 d), +0.10 (48), +0.07 (54), +0.14 (57-68 d); corr(T_eff, dOmega) +0.02 (30) -> +0.08 (59-68 d);
  corr(L, seed) -0.17 (30) -> -0.69 (41) -> -0.86..-0.91 (57-68 d).
- kirchhoff restart (40 d): corr(L, dOmega) +0.13 (41) -> +0.08 (57-59) -> +0.05 (68 d); corr(T_eff, dOmega) +0.10 -> +0.02 (68 d);
  corr(L, seed) -0.95..-0.96 throughout.
=> At N=32 the geometric (dOmega) part of the surface pattern is small (R^2 <= 0.02) and DECREASES with kirchhoff while it stays/rises in base;
the dominant surface pattern in both arms is the seed's angular function (growing). This does not reproduce your N=128 corr(L, dOmega) = -0.95
(at N=32 the seed response swamps it). So our earlier "kirchhoff does not reduce the imprint" was the seed pattern under the FC/ED masks;
by corr with dOmega, kirchhoff does reduce the geometric part. Score imprints by corr / regression on dOmega (seed-free runs best), not masked fc/ed.

Running (decisive, he_seed = 0, no MLT scaffold): static he_star_m1 RSG E HSE N=32 (+ cs_wellbalanced_src on), N=64, sp wedge twin; pure-hydro
cs_test iprob 13 strat N=32 (+ wbsrc, + exact seam ghosts), N=64. Metric: R^2 / corr with dOmega per shell + v_r rms vs time.
