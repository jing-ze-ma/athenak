# NOTE 2026-11-27zzh caltech -> orion: near-tangent fix, cheap variant, GPU cost on bsg + t9l62 (your 27zzh)

## Builds (CUDA 12.9 + g++13, star_vet; 02-03). NN2 = 3c681d9d1 + git am 0001-nearpath2-e218080e1.patch (local commit 0525b402d); NN0 = the same + git apply 0002-nn0-gate-off.patch, committed locally as d1b3c6477 on its own copy of the NN2 commit (ba4bd1be4)
- BUILD nn2 rc=0 wall=549 s TUs=280
- BUILD_OK nn2 md5=2df86d90bc5cedcbb8c1729314703fa3 prov=0525b402d5d1ff16284e7fdcc98852c6dffff46c
- GATE PASS nn2: athena_nn2_0525b402d md5 2df86d90 node hpc-sm-02-03 pass1=3.720e-11 log date 11-26 03:27:14 end: Terminating on cycle limit
- BUILD nn0 rc=0 wall=543 s TUs=280
- BUILD_OK nn0 md5=ad480205df01cd5d90e684bc23c89cf2 prov=d1b3c6477b8619787eb94ac2899dc81e67ab09c1
- GATE PASS nn0: athena_nn0_d1b3c6477 md5 ad480205 node hpc-sm-02-03 pass1=3.720e-11 log date 11-26 03:27:02 end: Terminating on cycle limit

## Arms
- Fresh tip templates (3c681d9d1), 300 cycles, 4 ranks, vet dump at cycle 300. bsg pair on 02-11, t9 pair on 02-17, each node to itself (7333073; the NB arms ran on 01-17, same node type).
- NB is reused from 27zzc. NN0's user.hst is byte-identical to NB's (and NBE's) for both stars, so NN0 serves as the NB reference for the dump comparisons (NB's own dumps were already deleted).
### NB-bsg: rc=0 ended: Terminating on cycle limit wall=270 s CW 0 FATAL 0 DT_RETRY 0; s/cycle median of the last elapsed intervals 0.6614; CW 0 FATAL 0 DT_RETRY 0
- `<rad_vet> vet_gd: rank 0 fenced seconds: source 1.266625e+00, sweeps 1.970938e+01 (of which per-shell halo+walls 7.431848e-01, of which MPI wait 4.787349e-01), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
- `<rad_vet> vet_gd: 3.070000e+02 builds checked; lateral reads clamped to the ghost band (rank 0, all sweeps): 0.000000e+00`
- `<rad_vet> vet_gd halo (rank 0): exchanges 7.554900e+04, bytes sent 1.518036e+11 (2.009340e+06 per exchange), post 2.382621e-02 s, partners 3; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
- L_top / L_in, last 50 cycles: mean 1.062188, min 1.060620, max 1.063609
### NN2-bsg: rc=0 ended: Terminating on cycle limit wall=291 s CW 0 FATAL 0 DT_RETRY 0; s/cycle median of the last elapsed intervals 0.7111; CW 0 FATAL 0 DT_RETRY 0
- `<rad_vet> vet_gd: rank 0 fenced seconds: source 1.579090e+00, sweeps 3.463798e+01 (of which per-shell halo+walls 7.372445e-01, of which MPI wait 4.808497e-01), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
- `<rad_vet> vet_gd: 3.070000e+02 builds checked; lateral reads clamped to the ghost band (rank 0, all sweeps): 0.000000e+00`
- `<rad_vet> vet_gd halo (rank 0): exchanges 7.554900e+04, bytes sent 1.518036e+11 (2.009340e+06 per exchange), post 2.522638e-02 s, partners 3; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
- L_top / L_in, last 50 cycles: mean 1.062188, min 1.060620, max 1.063609
### NN0-bsg: rc=0 ended: Terminating on cycle limit wall=273 s CW 0 FATAL 0 DT_RETRY 0; s/cycle median of the last elapsed intervals 0.6650; CW 0 FATAL 0 DT_RETRY 0
- `<rad_vet> vet_gd: rank 0 fenced seconds: source 1.470326e+00, sweeps 2.008970e+01 (of which per-shell halo+walls 8.672050e-01, of which MPI wait 4.803698e-01), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
- `<rad_vet> vet_gd: 3.070000e+02 builds checked; lateral reads clamped to the ghost band (rank 0, all sweeps): 0.000000e+00`
- `<rad_vet> vet_gd halo (rank 0): exchanges 7.554900e+04, bytes sent 1.518036e+11 (2.009340e+06 per exchange), post 2.554776e-02 s, partners 3; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
- L_top / L_in, last 50 cycles: mean 1.062188, min 1.060620, max 1.063609
- NN2 vs NN0 (= NB, user.hst byte-identical) at cycle 300 (dfloor 6.152e-14):
  - rho: all cells max 0.000e+00 at (r 1.4023e+12, th 1.0492, ph 0.0020), mean 0.000e+00; rho>10 dfloor max 0.000e+00 at (r 1.4023e+12, th 1.0492, ph 0.0020), mean 0.000e+00
  - T: all cells max 0.000e+00 at (r 1.4023e+12, th 1.0492, ph 0.0020), mean 0.000e+00; rho>10 dfloor max 0.000e+00 at (r 1.4023e+12, th 1.0492, ph 0.0020), mean 0.000e+00
  - |v|: all cells max 0.000e+00 at (r 1.4023e+12, th 1.0492, ph 0.0020), mean 0.000e+00; rho>10 dfloor max 0.000e+00 at (r 1.4023e+12, th 1.0492, ph 0.0020), mean 0.000e+00
### NB-t9: rc=0 ended: Terminating on cycle limit wall=856 s CW 0 FATAL 0 DT_RETRY 0; s/cycle median of the last elapsed intervals 1.9727; CW 0 FATAL 0 DT_RETRY 0
- `<rad_vet> vet_gd: rank 0 fenced seconds: source 1.940212e+00, sweeps 1.330236e+02 (of which per-shell halo+walls 1.320049e+02, of which MPI wait 1.376779e+01), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
- `<rad_vet> vet_gd: 3.070000e+02 builds checked; lateral reads clamped to the ghost band (rank 0, all sweeps): 0.000000e+00`
- `<rad_vet> vet_gd halo (rank 0): exchanges 2.041550e+05, bytes sent 2.203763e+11 (1.079456e+06 per exchange), post 4.465231e-01 s, partners 2; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
- L_top / L_in, last 50 cycles: mean 1.003203, min 1.002931, max 1.003453
### NN2-t9: rc=0 ended: Terminating on cycle limit wall=1009 s CW 0 FATAL 0 DT_RETRY 0; s/cycle median of the last elapsed intervals 2.5655; CW 0 FATAL 0 DT_RETRY 0
- `<rad_vet> vet_gd: rank 0 fenced seconds: source 2.413292e+00, sweeps 3.105363e+02 (of which per-shell halo+walls 3.094876e+02, of which MPI wait 1.570454e+01), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
- `<rad_vet> vet_gd: 3.070000e+02 builds checked; lateral reads clamped to the ghost band (rank 0, all sweeps): 0.000000e+00`
- `<rad_vet> vet_gd halo (rank 0): exchanges 2.041550e+05, bytes sent 2.203763e+11 (1.079456e+06 per exchange), post 3.934729e-01 s, partners 2; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
- L_top / L_in, last 50 cycles: mean 1.003204, min 1.002931, max 1.003460
### NN0-t9: rc=0 ended: Terminating on cycle limit wall=819 s CW 0 FATAL 0 DT_RETRY 0; s/cycle median of the last elapsed intervals 1.9863; CW 0 FATAL 0 DT_RETRY 0
- `<rad_vet> vet_gd: rank 0 fenced seconds: source 2.398003e+00, sweeps 1.360890e+02 (of which per-shell halo+walls 1.349849e+02, of which MPI wait 1.484701e+01), tables (rotation) 0.000000e+00; in-solve rebuilds 0.000000e+00 (0.000000e+00 s)`
- `<rad_vet> vet_gd: 3.070000e+02 builds checked; lateral reads clamped to the ghost band (rank 0, all sweeps): 0.000000e+00`
- `<rad_vet> vet_gd halo (rank 0): exchanges 2.041550e+05, bytes sent 2.203763e+11 (1.079456e+06 per exchange), post 3.842018e-01 s, partners 2; band-clamped lateral reads (ALL ranks, all sweeps) 0.000000e+00`
- L_top / L_in, last 50 cycles: mean 1.003203, min 1.002931, max 1.003453
- NN2 vs NN0 (= NB, user.hst byte-identical) at cycle 300 (dfloor 3.600e-19):
  - rho: all cells max 2.334e-03 at (r 6.9870e+12, th 1.6843, ph 2.2918), mean 1.254e-07; rho>10 dfloor max 1.559e-04 at (r 6.8918e+12, th 1.3407, ph 1.2425), mean 5.619e-08
  - T: all cells max 9.520e-05 at (r 5.3135e+13, th 0.8682, ph 0.6903), mean 2.941e-06; rho>10 dfloor max 2.144e-06 at (r 6.6142e+12, th 0.7885, ph 1.5739), mean 1.343e-08
  - |v|: all cells max 1.713e-02 at (r 6.9870e+12, th 1.0278, ph 2.9667), mean 6.852e-07; rho>10 dfloor max 2.123e-04 at (r 6.8918e+12, th 1.3775, ph 2.8747), mean 4.554e-08

## Reading
- **Cost** (s/cycle, median of the last elapsed intervals; fenced sweeps per 300 cycles; source):
  - bsg: NB 0.661 / NN0 0.665 (+0.5 %) / NN2 0.711 (**+7.5 %**). Sweeps 19.7 / 20.1 / 34.6 s (NN2 1.76x). Source (holds the bounds kernel) 1.27 / 1.47 / 1.58 s.
  - t9l62: NB 1.973 / NN0 1.986 (+0.7 %) / NN2 2.566 (**+30 %**). Sweeps 133 / 136 / 311 s (NN2 2.3x; booked under per-shell halo+walls on t9). Source 1.94 / 2.40 / 2.41 s.
  - So the code's presence plus the bounds kernel (NN0 - NB) costs 0.5-0.7 %. The test and the sampling (NN2 - NN0) cost +7 % (bsg) and +29 % (t9): above the 3 % target, but far below the first version (+50 % / +89 %).
- **Physics, NN2 vs NB:**
  - bsg: **bitwise identical** at cycle 300 (rho, T, |v|, E all max |rel| 0). No chord on bsg passes the new lateral test, yet the test still costs 1.76x in the sweep.
  - t9l62: tiny changes. rho > 10 dfloor at most 1.6e-4, T 2.1e-6, |v| 2.1e-4. E just outside the photosphere +2e-6 in the shell mean (per-cell max 2.3e-2 in the outer wind).
  - L_top/L_in is unchanged (bsg identical; t9 to 1e-6).
  - The cheap variant no longer gives the +0.8 % (bsg) / +1.5 % (t9) E increase of the first version. On these stars it does nearly nothing for E outside the photosphere.
- CW 0, FATAL 0, DT_RETRY 0; clamp lines 0 in all arms.
- E profiles (NN2 vs NN0):
```
### bsg: vet_e (E) NN2 vs NN0 (= NB) at the last cycle; shell means and per-cell relative difference
- all cells: max |rel| 0.000e+00, mean |rel| 0.000e+00
| r [cm] | shell-mean E NB | NN/NB - 1 (shell mean) | max |rel| in shell |
|---|---|---|---|
| 3.0274e+12 | 2.8374e+06 | +0.000e+00 | 0.000e+00 |
| 3.0817e+12 | 1.8926e+06 | +0.000e+00 | 0.000e+00 |
| 3.1361e+12 | 1.2337e+06 | +0.000e+00 | 0.000e+00 |
| 3.1904e+12 | 7.9229e+05 | +0.000e+00 | 0.000e+00 |
| 3.2448e+12 | 4.9588e+05 | +0.000e+00 | 0.000e+00 |
| 3.2991e+12 | 2.9563e+05 | +0.000e+00 | 0.000e+00 |
| 3.3535e+12 | 1.6733e+05 | +0.000e+00 | 0.000e+00 |
| 3.4078e+12 | 9.1062e+04 | +0.000e+00 | 0.000e+00 |
| 3.4622e+12 | 5.5026e+04 | +0.000e+00 | 0.000e+00 |
| 3.5166e+12 | 3.7053e+04 | +0.000e+00 | 0.000e+00 |
| 3.5709e+12 | 2.3144e+04 | +0.000e+00 | 0.000e+00 |
| 3.6253e+12 | 1.2165e+04 | +0.000e+00 | 0.000e+00 |
| 3.6796e+12 | 4.6457e+03 | +0.000e+00 | 0.000e+00 |
| 3.7340e+12 | 9.8456e+02 | +0.000e+00 | 0.000e+00 |
| 3.7883e+12 | 1.8693e+02 | +0.000e+00 | 0.000e+00 |
| 3.8590e+12 | 1.6751e+02 | +0.000e+00 | 0.000e+00 |
| 3.9677e+12 | 1.4927e+02 | +0.000e+00 | 0.000e+00 |
| 4.0764e+12 | 1.3574e+02 | +0.000e+00 | 0.000e+00 |
| 4.1851e+12 | 1.2486e+02 | +0.000e+00 | 0.000e+00 |
| 4.2938e+12 | 1.1570e+02 | +0.000e+00 | 0.000e+00 |
| 4.4025e+12 | 1.0776e+02 | +0.000e+00 | 0.000e+00 |
| 4.5112e+12 | 1.0077e+02 | +0.000e+00 | 0.000e+00 |
| 4.6199e+12 | 9.4525e+01 | +0.000e+00 | 0.000e+00 |
| 4.7286e+12 | 8.8907e+01 | +0.000e+00 | 0.000e+00 |
| 4.8373e+12 | 8.3818e+01 | +0.000e+00 | 0.000e+00 |
| 4.9460e+12 | 7.9183e+01 | +0.000e+00 | 0.000e+00 |
| 5.0547e+12 | 7.4936e+01 | +0.000e+00 | 0.000e+00 |
| 5.1634e+12 | 7.1023e+01 | +0.000e+00 | 0.000e+00 |
| 5.2721e+12 | 6.7402e+01 | +0.000e+00 | 0.000e+00 |
| 5.3808e+12 | 6.4035e+01 | +0.000e+00 | 0.000e+00 |
| 5.4895e+12 | 6.0878e+01 | +0.000e+00 | 0.000e+00 |
| 5.5547e+12 | 5.9046e+01 | +0.000e+00 | 0.000e+00 |
### t9: vet_e (E) NN2 vs NN0 (= NB) at the last cycle; shell means and per-cell relative difference
- all cells: max |rel| 2.278e-02, mean |rel| 1.606e-05
| r [cm] | shell-mean E NB | NN/NB - 1 (shell mean) | max |rel| in shell |
|---|---|---|---|
| 4.8256e+12 | 5.1150e+06 | +1.566e-11 | 9.791e-08 |
| 5.2035e+12 | 3.6912e+06 | -4.341e-11 | 6.795e-08 |
| 5.6109e+12 | 1.2122e+06 | +3.068e-10 | 1.035e-07 |
| 6.0503e+12 | 1.3021e+05 | -8.931e-09 | 6.015e-08 |
| 6.5241e+12 | 5.5145e+02 | +1.964e-06 | 6.464e-06 |
| 7.0350e+12 | 1.9300e+02 | -7.599e-07 | 4.071e-06 |
| 7.5859e+12 | 1.5115e+02 | -4.258e-07 | 2.414e-06 |
| 8.1800e+12 | 1.2242e+02 | -4.995e-07 | 2.740e-06 |
| 8.8205e+12 | 1.0102e+02 | -8.090e-07 | 2.698e-06 |
| 9.5113e+12 | 8.4209e+01 | -1.040e-06 | 2.090e-06 |
| 1.0256e+13 | 7.0630e+01 | -1.504e-06 | 2.778e-06 |
| 1.1059e+13 | 5.9530e+01 | -1.516e-06 | 3.515e-06 |
| 1.1925e+13 | 5.0380e+01 | -3.302e-08 | 4.738e-06 |
| 1.2859e+13 | 4.2766e+01 | +5.290e-07 | 5.324e-06 |
| 1.3866e+13 | 3.6388e+01 | -2.591e-06 | 4.884e-06 |
| 1.4952e+13 | 3.1023e+01 | -5.962e-06 | 1.174e-05 |
| 1.6123e+13 | 2.6491e+01 | -3.516e-06 | 9.109e-06 |
| 1.7385e+13 | 2.2651e+01 | +2.828e-06 | 2.061e-05 |
| 1.8747e+13 | 1.9389e+01 | +3.762e-06 | 3.045e-05 |
| 2.0215e+13 | 1.6611e+01 | -4.921e-06 | 2.594e-05 |
| 2.1798e+13 | 1.4241e+01 | -1.517e-05 | 5.293e-05 |
| 2.3505e+13 | 1.2215e+01 | -1.503e-05 | 3.770e-05 |
| 2.5345e+13 | 1.0482e+01 | -4.608e-06 | 3.129e-05 |
| 2.7330e+13 | 8.9979e+00 | +3.143e-06 | 3.294e-05 |
| 2.9470e+13 | 7.7273e+00 | -2.141e-06 | 6.855e-05 |
| 3.1778e+13 | 6.6386e+00 | -1.335e-05 | 6.673e-05 |
| 3.4267e+13 | 5.7047e+00 | -1.002e-05 | 1.036e-04 |
| 3.6950e+13 | 4.9027e+00 | +2.646e-06 | 1.107e-04 |
| 3.9844e+13 | 4.2142e+00 | +2.878e-07 | 1.283e-04 |
| 4.2964e+13 | 3.6230e+00 | -2.455e-06 | 4.039e-04 |
| 4.6328e+13 | 3.1151e+00 | +5.655e-06 | 1.696e-03 |
| 4.9956e+13 | 2.6785e+00 | +6.237e-06 | 6.206e-03 |
| 5.3868e+13 | 2.3033e+00 | +3.941e-06 | 1.289e-02 |
| 5.6130e+13 | 2.1211e+00 | +5.537e-06 | 1.726e-02 |
```
- Files: docs/handover/caltech_results/NEARPATH2/<arm>/.
