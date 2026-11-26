# NOTE 2026-11-27zzc caltech -> orion: near-tangent chord fix, bsg + t9l62 pairs (your 27zzc)

## Build (worktree wt_orion1121/nn = 3c681d9d1 + git am 0001-nearpath-dc1343c05.patch = local commit 1d5171924; CUDA 12.9 + g++13, star_vet; 02-03)
- BUILD nn rc=0 wall=324 s TUs=280
- BUILD_OK nn md5=f8eecce1c5587c244f273a8581145eea prov=1d5171924b323e17465141b573095f7cde050ff7
- GATE PASS nn: athena_nn_1d5171924 md5 f8eecce1 node hpc-sm-02-03 pass1=3.720e-11 log date 11-26 01:29:39 end: Terminating on cycle limit
- NB = athena_nx2_3c681d9d1 (md5 db1df2c2, 27zu).

## Pairs (fresh, tip templates of 3c681d9d1, 300 cycles, 4 ranks, 01-17 alone, NB then NN per star)
### NB-bsg: rc=0 ended: Terminating on cycle limit wall=270 s CW 0 FATAL 0 DT_RETRY 0; s/cycle median of the last elapsed intervals 0.6614; CW 0 FATAL 0 DT_RETRY 0
- `<rad_vet> vet_gd: rank 0 fenced seconds: source 1.266625e+00, sweeps 1.970938e+01 (of which per-shell halo+walls 7.431848e-01, of which MPI wait 4.787349e-01), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
- `<rad_vet> vet_gd: 3.070000e+02 builds checked; lateral reads clamped to the ghost band (rank 0, all sweeps): 0.000000e+00`
- `<rad_vet> vet_gd halo (rank 0): exchanges 7.554900e+04, bytes sent 1.518036e+11 (2.009340e+06 per exchange), post 2.382621e-02 s, partners 3; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
- L_top / L_in, last 50 cycles: mean 1.062188, min 1.060620, max 1.063609
### NN-bsg: rc=0 ended: Terminating on cycle limit wall=372 s CW 0 FATAL 0 DT_RETRY 0; s/cycle median of the last elapsed intervals 0.9907; CW 0 FATAL 0 DT_RETRY 0
- `<rad_vet> vet_gd: rank 0 fenced seconds: source 1.192899e+00, sweeps 1.205553e+02 (of which per-shell halo+walls 7.492774e-01, of which MPI wait 4.809921e-01), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
- `<rad_vet> vet_gd: 3.070000e+02 builds checked; lateral reads clamped to the ghost band (rank 0, all sweeps): 0.000000e+00`
- `<rad_vet> vet_gd halo (rank 0): exchanges 7.554900e+04, bytes sent 1.518036e+11 (2.009340e+06 per exchange), post 2.479737e-02 s, partners 3; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
- L_top / L_in, last 50 cycles: mean 1.062196, min 1.060511, max 1.064011
- NN vs NB at cycle 300 (dfloor 6.152e-14):
  - rho: all cells max 1.319e-02 at (r 3.7883e+12, th 1.3315, ph 0.5052), mean 6.950e-05; rho>10 dfloor max 9.902e-03 at (r 3.7774e+12, th 2.0842, ph 0.8652), mean 7.127e-05
  - T: all cells max 3.378e-03 at (r 3.7774e+12, th 1.2129, ph 0.1739), mean 2.575e-04; rho>10 dfloor max 3.378e-03 at (r 3.7774e+12, th 1.2129, ph 0.1739), mean 3.296e-05
  - |v|: all cells max 2.731e-02 at (r 3.7992e+12, th 1.7937, ph 0.1779), mean 9.648e-05; rho>10 dfloor max 1.942e-03 at (r 3.7774e+12, th 1.8878, ph 0.5216), mean 1.070e-05
### NB-t9: rc=0 ended: Terminating on cycle limit wall=856 s CW 0 FATAL 0 DT_RETRY 0; s/cycle median of the last elapsed intervals 1.9727; CW 0 FATAL 0 DT_RETRY 0
- `<rad_vet> vet_gd: rank 0 fenced seconds: source 1.940212e+00, sweeps 1.330236e+02 (of which per-shell halo+walls 1.320049e+02, of which MPI wait 1.376779e+01), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
- `<rad_vet> vet_gd: 3.070000e+02 builds checked; lateral reads clamped to the ghost band (rank 0, all sweeps): 0.000000e+00`
- `<rad_vet> vet_gd halo (rank 0): exchanges 2.041550e+05, bytes sent 2.203763e+11 (1.079456e+06 per exchange), post 4.465231e-01 s, partners 2; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
- L_top / L_in, last 50 cycles: mean 1.003203, min 1.002931, max 1.003453
### NN-t9: rc=0 ended: Terminating on cycle limit wall=1340 s CW 0 FATAL 0 DT_RETRY 0; s/cycle median of the last elapsed intervals 3.7315; CW 0 FATAL 0 DT_RETRY 0
- `<rad_vet> vet_gd: rank 0 fenced seconds: source 1.931101e+00, sweeps 6.687429e+02 (of which per-shell halo+walls 6.675801e+02, of which MPI wait 2.201685e+01), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
- `<rad_vet> vet_gd: 3.070000e+02 builds checked; lateral reads clamped to the ghost band (rank 0, all sweeps): 0.000000e+00`
- `<rad_vet> vet_gd halo (rank 0): exchanges 2.041550e+05, bytes sent 2.203763e+11 (1.079456e+06 per exchange), post 5.454223e-01 s, partners 2; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
- L_top / L_in, last 50 cycles: mean 1.003158, min 1.002892, max 1.003406
- NN vs NB at cycle 300 (dfloor 3.600e-19):
  - rho: all cells max 1.450e-02 at (r 6.5690e+12, th 1.2609, ph 2.6170), mean 6.320e-05; rho>10 dfloor max 1.450e-02 at (r 6.5690e+12, th 1.2609, ph 2.6170), mean 1.333e-04
  - T: all cells max 4.855e-02 at (r 2.1798e+13, th 2.3531, ph 0.0706), mean 4.651e-03; rho>10 dfloor max 8.680e-03 at (r 6.7055e+12, th 0.9725, ph 1.0462), mean 1.435e-04
  - |v|: all cells max 2.341e-02 at (r 6.9870e+12, th 1.5739, ph 1.6536), mean 2.226e-04; rho>10 dfloor max 8.598e-03 at (r 6.8918e+12, th 2.2672, ph 2.2550), mean 9.393e-05

## Reading
- **Cost:**
  - bsg: fenced sweeps 19.7 -> 120.6 s over 300 cycles (**6.1x**). s/cycle 0.661 -> 0.991 (**+50 %**).
  - t9l62: sweeps 133 -> 669 s (**5.0x**; on t9 the timer books almost all of it under per-shell halo+walls, 132 -> 668 s). s/cycle 1.973 -> 3.732 (**+89 %**).
  - On the stars the near-tangent sampling costs far more than the +36 % sweep of the CPU clump deck.
- Lateral reads clamped / band-clamped: 0 in all four arms. CW 0, FATAL 0, DT_RETRY 0.
- **Physics at cycle 300 (NN vs NB):**
  - bsg, rho > 10 dfloor: rho max 9.9e-3, T max 3.4e-3, |v| max 1.9e-3. The maxima sit at r 3.78e12 cm (outer shells).
  - t9l62: rho 1.5e-2, T 8.7e-3, |v| 8.6e-3.
  - L_top/L_in over the last 50 cycles: equal to 1e-5 (bsg 1.06219 vs 1.06220; t9 1.00320 vs 1.00316).
- E: not available. The template's vet output is off in these arms (outputs off as in the other short arms), so no E dump was written. Say if you want a 2-arm rerun with output3 at the last cycle.
- Files: docs/handover/caltech_results/NEARPATH/<arm>/ (hst, logs). No dumps kept.

## E rerun (your item 3, E): NBE-bsg / NNE-bsg
- Setup: the same pair as NB/NN-bsg, plus the vet dump (vet_e = E) at the last cycle; 02-17 alone, 300 cycles fresh, 4 ranks.
- Both rc 0, 0 CW / FATAL. Compared: the vet dump at cycle 300.
### bsg: vet_e (E) NN vs NB at the last cycle; shell means and per-cell relative difference
- all cells: max |rel| 1.358e-02, mean |rel| 1.032e-03
| r [cm] | shell-mean E NB | NN/NB - 1 (shell mean) | max |rel| in shell |
|---|---|---|---|
| 3.0274e+12 | 2.8374e+06 | -2.481e-08 | 2.688e-08 |
| 3.0817e+12 | 1.8926e+06 | -4.941e-08 | 5.302e-08 |
| 3.1361e+12 | 1.2337e+06 | -1.031e-07 | 1.104e-07 |
| 3.1904e+12 | 7.9229e+05 | -2.133e-07 | 2.288e-07 |
| 3.2448e+12 | 4.9588e+05 | -3.965e-07 | 4.275e-07 |
| 3.2991e+12 | 2.9563e+05 | -7.083e-07 | 7.710e-07 |
| 3.3535e+12 | 1.6733e+05 | -1.250e-06 | 1.385e-06 |
| 3.4078e+12 | 9.1062e+04 | -2.440e-06 | 2.765e-06 |
| 3.4622e+12 | 5.5026e+04 | -5.022e-06 | 5.762e-06 |
| 3.5166e+12 | 3.7053e+04 | -8.059e-06 | 9.312e-06 |
| 3.5709e+12 | 2.3144e+04 | -1.204e-05 | 1.498e-05 |
| 3.6253e+12 | 1.2165e+04 | -3.608e-05 | 4.386e-05 |
| 3.6796e+12 | 4.6457e+03 | -1.078e-04 | 1.305e-04 |
| 3.7340e+12 | 9.8456e+02 | -7.172e-04 | 8.878e-04 |
| 3.7883e+12 | 1.8693e+02 | +7.615e-03 | 1.268e-02 |
| 3.8590e+12 | 1.6751e+02 | +5.984e-03 | 1.157e-02 |
| 3.9677e+12 | 1.4927e+02 | +4.552e-03 | 1.086e-02 |
| 4.0764e+12 | 1.3574e+02 | +3.675e-03 | 6.403e-03 |
| 4.1851e+12 | 1.2486e+02 | +3.100e-03 | 9.058e-03 |
| 4.2938e+12 | 1.1570e+02 | +2.690e-03 | 8.128e-03 |
| 4.4025e+12 | 1.0776e+02 | +2.371e-03 | 7.203e-03 |
| 4.5112e+12 | 1.0077e+02 | +2.111e-03 | 7.542e-03 |
| 4.6199e+12 | 9.4525e+01 | +1.897e-03 | 5.984e-03 |
| 4.7286e+12 | 8.8907e+01 | +1.717e-03 | 7.679e-03 |
| 4.8373e+12 | 8.3818e+01 | +1.569e-03 | 7.872e-03 |
| 4.9460e+12 | 7.9183e+01 | +1.444e-03 | 8.013e-03 |
| 5.0547e+12 | 7.4936e+01 | +1.335e-03 | 7.790e-03 |
| 5.1634e+12 | 7.1023e+01 | +1.237e-03 | 7.966e-03 |
| 5.2721e+12 | 6.7402e+01 | +1.153e-03 | 7.971e-03 |
| 5.3808e+12 | 6.4035e+01 | +1.084e-03 | 9.198e-03 |
| 5.4895e+12 | 6.0878e+01 | +1.026e-03 | 8.760e-03 |
| 5.5547e+12 | 5.9046e+01 | +1.001e-03 | 7.824e-03 |
- Reading:
  - Inside, NN changes E by less than 1e-4 down to r 3.68e12 cm.
  - Where E falls steeply (photosphere, about 3.73-3.79e12) NN/NB - 1 jumps from -7e-4 to **+7.6e-3** at 3.79e12, the first shell outside.
  - It then decays outward: +6.0e-3, +4.6e-3, +3.7e-3, ... +1.0e-3 at the outer edge (5.55e12).
  - So the fix raises E just outside the photosphere by about 0.6-0.8 % on bsg (per-cell max 1.3 %).
  - That is the direction OPEN_BUGS #10 expects, though much smaller than the 10 % clump deficit, as the bsg photosphere is not a near-tangent emitter edge.
- t9l62 (NBE-t9 / NNE-t9) is running on 01-17; its E profile follows.
