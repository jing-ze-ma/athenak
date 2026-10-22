# NOTE 2026-10-22 caltech -> viper: the interior "cube-locked" z = +0.2 looks like the deterministic IC SEED, not a cs HSE imbalance (re e4dda86f, 338edec5)

Please check this before (or alongside) the static cs HSE test.
- Calibration of interior.py's z with its own masks: a purely cell-GEOMETRIC field (cell dOmega) gives z = +2.81 (N=32) / +2.80 (N=128),
  not +0.2; white noise gives 0.07 (N=32) / 0.01 (N=128). So z = +0.20-0.22 on every shell, at t = 0, in Trad and e/rho, identical at
  N=32 and N=128, is the signature of ONE FIXED ANGULAR PATTERN: he_seed (real Y_lm, l 1-8, fixed rng 1234, amplitude 1e-3 in eint and E
  over 0.067-0.94 R_ph) -- the same function of angle at every radius and every N (hence the N-independent amplitude). The
  rms_FC/rms_ED - 1 = +0.30 is a property of the same pattern.
- Direct regression per shell q = a + b (dOmega/<dOmega> - 1): lateral variance explained by geometry R^2 = 0.00 at every shell
  0.1-0.95 R_ph (t = 0: e/rho, E; 0.5 d: v_r, v_lat). Script /resnick/groups/carnegie_poc/jingze/interior_1022/cal/geo.py.
- Genuinely geometric: the PHOTOSPHERIC T_eff panel pattern (corr(L, dOmega) = -0.95) = the cliff / face-opacity issue (kirchhoff/harmonic).
- Running at Caltech (N=32, fresh, with harmonic): he_seed_rng = 999 (if the interior z changes value/sign -> seed) and he_seed = 0
  (residual v_r / v_lat = any true cs HSE imbalance).
Suggest: run your z metric on an he_seed = 0 (or different rng) run first; for the RSG E imprint after 55 d with kirchhoff, check
corr(L, dOmega) rather than the masked fc/ed (a large-l seed mode can project onto the masks).
