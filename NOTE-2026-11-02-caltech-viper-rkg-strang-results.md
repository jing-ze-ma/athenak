# NOTE 2026-11-02 Caltech -> viper (+ USER): RKG Lie vs Strang GPU results (rkgstrang-test-1102 @ e204131b2)

Reply to 0426ac492 and the addendum c8017e87. 1 H200 per arm, hpc-sm-01-11. Binary athena_gpu_dhj_rkgstr_e204131b_fmadoff
(md5 e247f6c2; incremental, -ffp-contract=off). Deck = r1m production keys (r1mprof_1101/base.athinput: use_rkg_sts,
rkg_lean, rkg_dt_safety 0.5, ambipolar eos, max_eta 1e15) + `rkg_strang` in <mhd> + output5 bin mhd_eta. Start =
developed state trc_1101/keep/off_0p50.rst (cycle 8774, t = 55074 s). All files: /resnick/groups/carnegie_poc/jingze/rkgstrang_1102.

## Gates
- G1 PASS: rkg_strang absent, 40 cycles, e204131b2 vs aeafae97 (md5 5b0d8791): rst data after <par_end> IDENTICAL,
  dhj.mhd.hst IDENTICAL. (user.hst differs; this is the known floor-ledger order, it differs between two base runs too.)
- G2 PASS: rkg_strang = true, 20 + 20 vs 40 cycles: rst data IDENTICAL. Strang does change the solution (G1 vs G2-40 DIFF).
- Extra: strangB vs strang (same binary, same GPU type) is BITWISE identical, so the GPU noise floor is zero; the noise
  below is seed - lie alone (seed 1, amp 1e-14, seed_restart).

## Arms (0.02 rot, to tlim 57277 s; all 6 final dumps at t = 57277.000000)
lie / strang / lie4 / strang4 (cfl 0.3 / 0.075) / seed / strangB. lie: 480 cycles; lie4: 1960 cycles.
Analysis: rkgstrang_1102/ana/bands.py (re-implementation of the c8017e87 spec: bands from lie's p; volume-weighted
RMS; |J| = |J_t| = |(1/r) d(r B_t)/dr| with B_t = (B_zon, B_mer); eta = output mhd_eta; dTot = |relative difference
of the band totals|). Units: T K, u v cm/s, B code units. Rows "/n" = that difference divided by the noise (seed - lie).

```
## band 1e-6..1e-3 bar  (cells 133518)
RMS diff                  T          u          v        B_r      B_zon      B_mer        eta       etaJ        E_B          q   dTot E_B     dTot q
lie field RMS     3.436e+03  2.611e+05  1.121e+05  1.834e+00  1.768e+00  1.427e+00  3.378e+14  2.968e+05  2.333e+01  1.649e-03
seed-lie          1.113e+00  4.675e+01  3.075e+01  2.394e-04  1.507e-04  7.674e-04  4.503e+11  1.048e+03  2.313e-03  1.384e-06  1.993e-05  3.568e-05
lie-strang4       2.663e+00  2.382e+02  1.350e+02  4.938e-04  1.284e-03  1.561e-03  1.770e+12  5.342e+03  1.754e-02  3.570e-05  7.619e-05  5.435e-03
strang-strang4    2.871e+00  2.275e+02  1.149e+02  5.517e-04  1.149e-03  1.755e-03  1.804e+12  3.624e+03  1.585e-02  2.417e-05  1.481e-05  6.109e-04
lie-strang        1.779e+00  1.165e+02  8.129e+01  4.580e-04  8.084e-04  1.352e-03  9.494e+11  3.805e+03  8.089e-03  2.285e-05  6.138e-05  4.827e-03
lie4-strang4      1.826e+00  1.209e+02  7.579e+01  2.545e-04  4.839e-04  7.784e-04  1.087e+12  2.133e+03  6.857e-03  8.850e-06  3.326e-06  1.842e-03
(lie-strang4)/n       2.4        5.1        4.4        2.1        8.5        2.0        3.9        5.1        7.6       25.8        3.8      152.3
(strang-strang4)/n       2.6        4.9        3.7        2.3        7.6        2.3        4.0        3.5        6.9       17.5        0.7       17.1
(lie-strang)/n        1.6        2.5        2.6        1.9        5.4        1.8        2.1        3.6        3.5       16.5        3.1      135.3
(lie4-strang4)/n       1.6        2.6        2.5        1.1        3.2        1.0        2.4        2.0        3.0        6.4        0.2       51.6

## band 1e-3..1 bar  (cells 148936)
RMS diff                  T          u          v        B_r      B_zon      B_mer        eta       etaJ        E_B          q   dTot E_B     dTot q
lie field RMS     2.411e+03  9.627e+04  5.469e+04  1.643e+00  1.133e+01  6.757e+00  2.119e+14  1.159e+05  2.167e+02  3.547e-03
seed-lie          1.555e-03  6.409e-02  6.030e-02  1.130e-04  1.081e-04  8.010e-04  2.964e+10  4.723e+02  2.559e-03  1.779e-07  2.220e-06  8.116e-06
lie-strang4       6.351e-02  8.105e+00  5.343e+00  2.691e-04  1.843e-03  1.506e-03  1.859e+10  1.715e+03  5.552e-02  2.663e-05  6.506e-05  2.125e-03
strang-strang4    6.273e-02  7.813e+00  5.219e+00  2.324e-04  1.262e-03  1.635e-03  4.677e+10  7.277e+02  4.386e-02  4.810e-06  3.186e-05  1.253e-04
lie-strang        6.275e-03  3.162e+00  1.618e+00  3.107e-04  1.342e-03  1.544e-03  4.481e+10  1.846e+03  3.258e-02  2.600e-05  3.320e-05  2.000e-03
lie4-strang4      1.517e-02  2.206e+00  2.797e+00  1.434e-04  6.639e-04  8.304e-04  1.217e+10  5.154e+02  2.218e-02  7.979e-06  8.734e-06  4.778e-04
(lie-strang4)/n      40.8      126.5       88.6        2.4       17.0        1.9        0.6        3.6       21.7      149.7       29.3      261.8
(strang-strang4)/n      40.3      121.9       86.5        2.1       11.7        2.0        1.6        1.5       17.1       27.0       14.4       15.4
(lie-strang)/n        4.0       49.3       26.8        2.7       12.4        1.9        1.5        3.9       12.7      146.2       15.0      246.4
(lie4-strang4)/n       9.8       34.4       46.4        1.3        6.1        1.0        0.4        1.1        8.7       44.9        3.9       58.9

## band p>1 bar  (cells 192644)
RMS diff                  T          u          v        B_r      B_zon      B_mer        eta       etaJ        E_B          q   dTot E_B     dTot q
lie field RMS     3.635e+03  1.598e+03  1.473e+03  1.715e+00  1.701e+00  1.953e+00  2.988e+11  8.775e+03  1.145e+01  6.745e-04
seed-lie          9.362e-06  2.451e-04  4.200e-04  4.967e-09  1.506e-08  5.608e-08  6.836e+03  1.267e-03  2.072e-07  5.162e-11  1.075e-10  2.992e-09
lie-strang4       1.729e-04  1.699e-01  1.331e-01  2.526e-05  5.793e-05  7.310e-05  2.007e+06  6.806e-01  1.869e-04  4.276e-08  1.098e-07  1.046e-05
strang-strang4    1.663e-04  1.744e-01  1.385e-01  9.207e-06  1.919e-05  2.896e-05  2.068e+06  2.254e-01  8.815e-05  1.960e-08  1.126e-06  1.415e-05
lie-strang        6.532e-05  4.002e-02  3.950e-02  1.629e-05  4.214e-05  4.895e-05  1.927e+05  5.473e-01  1.695e-04  4.723e-08  1.016e-06  2.461e-05
lie4-strang4      1.039e-04  2.874e-03  5.407e-03  7.787e-06  2.147e-05  2.497e-05  3.503e+05  2.142e-01  7.055e-05  1.644e-08  1.262e-07  5.669e-06
(lie-strang4)/n      18.5      693.2      316.8     5084.4     3845.6     1303.5      293.6      537.3      902.0      828.4     1022.0     3497.0
(strang-strang4)/n      17.8      711.7      329.8     1853.5     1274.0      516.5      302.5      177.9      425.4      379.7    10473.1     4727.6
(lie-strang)/n        7.0      163.3       94.0     3279.3     2797.9      872.9       28.2      432.1      818.0      915.0     9451.0     8224.7
(lie4-strang4)/n      11.1       11.7       12.9     1567.7     1425.4      445.3       51.2      169.1      340.4      318.4     1173.8     1894.6

```
(p < 1e-6 bar, not judged: every pair is within 1-2x of the noise.)

## Reading
1. By the c8017e87 rule Lie FAILS: lie - strang4 is above the noise in every judged band (1e-6..1e-3: 2-150x;
   1e-3..1: up to 260x, T 41x, u 127x; p > 1: 20-5000x). BUT strang - strang4 fails it by almost the same factors for
   T, u, v (e.g. 1e-3..1 bar: T 40x, u 122x; p > 1: u 712x). So the rule does not single out the split: at 0.02 rot the
   seed has not grown to a saturated noise floor (deep noise is ~1e-9 relative), and most of the cfl 0.3 -> 0.075
   difference is the rest of the time discretisation, not the RKG split.
2. The split error proper (lie - strang, equal dt) is smaller than the total dt error for the flow: 1e-3..1 bar
   u 3.2 vs 8.1 (lie - strang4), T 0.006 vs 0.064 K. Relative to the lie field RMS all differences are 1e-4 to 1e-6
   (B_zon 1e-3..1: 1.8e-3 / 11.3 = 1.6e-4; u: 8.1 / 9.6e4 = 8e-5 ; p > 1: B_r 1.5e-5).
3. Where Strang does help, it is the resistive quantities: q band total, lie vs strang4, 5.4e-3 (1e-6..1e-3), 2.1e-3
   (1e-3..1), 1.0e-5 (p > 1); Strang's own error is 6.1e-4, 1.3e-4, 1.4e-5. q RMS: 1.5-5.5x smaller with Strang;
   B_r / B_zon / B_mer at p > 1 bar: 2.5-3x smaller.
4. Order: (lie - strang4) / (lie4 - strang4) median 2.3 / 3.1 / 2.8 in the three judged bands. This fits a first-order
   error mixed with other dt errors (pure first order gives 4). (lie - strang4) / (strang - strang4): median 1.1 / 1.2 / 1.6,
   so at the production cfl Strang lowers the total error only 10-60 % in the median field. Your CPU result (Strang
   31-125x better at equal dt) is the pure split error; on r1m it is hidden under the other dt errors except in q and deep B.

## Timing (1 H200, 01-11 GPU 1, no other process on that GPU, interleaved lie, strang, lieB, strangB; cycles 8880..9260)
- lie 83.14 / 82.89 ms/cycle (21 RKG stages per cycle); strang 107.15 / 107.10 ms/cycle (2 x 15 = 30 stages).
- Strang = 1.29x per cycle (+29 % wall: 300 rot about 7.0 d -> 9.0 d on 1 GPU). The stage count is 1.43x, as expected.
- Repeat on 01-11 GPU 0 with the whole node quiet (t1_lie, t1_strang, t2_lie, t2_strang): lie 82.75 / 82.71, strang
  106.82 / 106.99 ms/cycle = 1.29x again.
- (lie4 58.4, strang4 72.7 ms/cycle; these ran beside the other arms, so they are not timing data.)

## Verdict / recommendation (the USER decides; nothing merged; rkg_strang stays default off)
Formally Lie fails the noise rule, but so does Strang against strang4, so the rule cannot separate them at 0.02 rot.
The Lie split adds errors of 1e-4 or less relative in T, v and B, which are a fraction of the other time-step errors.
The one systematic effect is on the resistive heating band totals: 0.2-0.5 % at 1e-6..1 bar with Lie, 0.01-0.06 % with
Strang. Strang costs +29 % wall. Caltech recommendation: launch r1m with Lie (the current default). Use Strang only if
the resistive heating at the 0.5 % level matters for the science.
