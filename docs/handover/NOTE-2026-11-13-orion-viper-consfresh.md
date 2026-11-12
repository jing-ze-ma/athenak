# NOTE 2026-11-13 Orion -> viper: consfresh. (b) d6fast3 gate PASS; (a) audit = viper GPU to 3 digits; neutrality PASS

## Builds (git archive + kokkos from the Orion src_vfix copy; patches applied with git apply, both clean)
| tree | commit | patch | athena md5 |
|---|---|---|---|
| f | fbb5c663825e3861140c34430dc5f804b6d45200 | - | dbbaf51588ab44268ab0aa8e2987cd1b |
| p | d3b3f2e5531a915c69b8ede3347fd1aca8951860 | - | 72ccef10641b7ad34e1b355d246f7fb3 |
| fa | fbb5c663 | consaudit_on_fbb5c663.patch | 8ddfeb83f3f4bfe59fd1e6c309d64019 |
| pa | d3b3f2e5 | consaudit_on_d3b3f2e55.patch | cb82bd1f82190da10a8896e1ed81d781 |

Build: mpicxx gcc/13 + openmpi/4.1.8, MPI ON, PROBLEM=he_star_m1, Release, -ffp-contract=off. `strings fa | grep -c dbg_cons_audit` = 1.
Tables: hgpos-1109/files, md5 ic f9678c2a, planck 23ee543f, rosseland e2868b75 (= viper). 24 ranks per run, `mpirun --bind-to none`, binding off, OMP 1.

**Deviation (needed):** the $AM and $DG keys cannot go on the command line, because none of them is in the decks. AthenaK FATALs with
"Parameter 'am_cons_diag' in block 'hydro' on command line not found" (and the same for dbg_energy_audit). I inserted them into each run's copy of
the deck directly after the <hydro> / <rad_m1> header (addkeys.py). time/nlim stays on the command line (nlim is in the deck). Smoke (3 cycles,
4 runs) then gave rc 0 and 0 FATAL, with AMCS_CONS / CAUDITH lines printing.
G and N ran on orion060 (4 x 24 ranks), A on orion101 (2 x 24), p.exclusive.

## (b) Gate, 40 cycles from the IC (plain binaries + AM print keys)
| run | rc | FATAL | NOT ADMISSIBLE | NOT POSITIVE | stage / BE / fallbacks | Picard solves, mean / max / NC | wall s |
|---|---|---|---|---|---|---|---|
| G_D6 | 0 | 0 | 0 | 0 | 39 / 1 / 0 | 7.900000e+01, 3.278481e+00 / 5.000000e+00 / 0.000000e+00 | 288 |
| G_PROD | 0 | 0 | 0 | 0 | 39 / 1 / 0 | 8.300000e+01, 3.265060e+00 / 5.000000e+00 / 0.000000e+00 | 272 |

- **D6 gate PASS** (rc 0, FATAL 0, NOT ADMISSIBLE 0, NOT POSITIVE 0; 39 stage steps, 1 BE, 0 stage fallbacks).
- **face_fct (hesdirk2 summary, D6):** calls 78, **limited cells 50881**, max rounds 8, withheld E x vol 1.351e40, left (Q<=0) 0.
  face_mmat: checked passes 259, limited cells 0. PROD has no face_fct.

### amlz_analyze.py G_D6 G_PROD
```
G_D6 AMH stages 80 sum_dT_z 2.371e+38 max|resid_z| 2.535e+29 Lz 1.077552e+29 -> -4.619820e+38
G_D6 AMCS solves 79 max|resid_z| 1.958e+26 max|resid_xy| 2.446e+39 mean scale_z 8.771e+40
G_D6 AMCS_DIAG sums |dLr_z| |dLg_z| |dtTpl_z| |dtTrad_z| |del_z| |del_z|/(rhoVR) a_f |Lg_z|: 6.297e+40 4.859e+40 8.929e+40 6.486e+40 4.770e+39 2.791e+01 9.948e+38 5.092e+41
G_PROD AMH stages 80 sum_dT_z 1.658e+38 max|resid_z| 2.080e+29 Lz 8.254908e+28 -> -2.531348e+38
G_PROD AMCS solves 83 max|resid_z| 5.077e+25 max|resid_xy| 2.446e+39 mean scale_z 1.126e+41
G_PROD AMCS_DIAG sums |dLr_z| |dLg_z| |dtTpl_z| |dtTrad_z| |del_z| |del_z|/(rhoVR) a_f |Lg_z|: 6.091e+40 4.872e+40 9.081e+40 6.723e+40 4.789e+39 2.796e+01 -3.050e+38 5.098e+41
```

- **gas L_z after 40 cycles: D6 -4.62e38, PROD -2.53e38** (viper GPU PROD at 30 cycles: -4.5e38). **AMCS_DIAG |del_z|: D6 4.77e39, PROD 4.79e39**
  (viper GPU PROD 30 cycles: 2.2e39; 40 cycles here). D6 has the production bookkeeping: |del_z| equal to PROD within 0.4 %, |L_z| within 2x of PROD
  and of the same order as viper's production reference. AMH max |resid_z| 2.5e29 (D6) vs 2.1e29 (PROD); AMCS max |resid_z| 2.0e26 vs 5.1e25.
- **Neutrality, part 1 (G vs N, 40 rows):** the plain build with the AM print keys (G) = the diagnostic build with both keys off (N), md5-identical:
  D6 hydro 0b2bef81 = 0b2bef81, user e2be5729 = e2be5729; PROD hydro 48502d66 = 48502d66, user e458e658 = e458e658.
  The A (keys on) comparison follows with part (a).

## (a) Conservation audit, 100 cycles (diagnostic builds fa / pa, DG keys inserted into the deck): **CPU reproduces viper's GPU values to 3 digits**
A_D6 and A_PROD: rc 0, FATAL 0 (orion101, 2 x 24 ranks).

| quantity | viper GPU | Orion CPU D6 | Orion CPU PROD |
|---|---|---|---|
| RE_b cycle 2 | 2.23e38 | 2.226e38 | 2.226e38 |
| RE_b cycle 10 | 1.02e38 | 1.021e38 | 1.021e38 |
| RE_b cycles 97 / 98 / 99 | ~1e34 | 1.147e34 / 9.059e33 / 6.942e33 | 1.147e34 / 9.060e33 / 6.942e33 |
| radiation boundary BrB, cycle 2 -> 10 -> 97 | 7e37 - 4.1e39, 1.2e39 | 6.98e37 -> 4.05e39 -> 1.26e39 | 7.28e37 -> 4.15e39 -> 1.26e39 |
| accumulated L_z(gas+rad), 99 cycles | D6 -4.78e40, PROD -4.50e40 | **-4.775e40** | **-4.503e40** |
| P_r(gas+rad) change | 1.174e35 | 1.174e35 | 1.174e35 |
| fix segments Lz / Pr | 0 | **0 / 0 exactly** | **0 / 0 exactly** |
| fix segment E (sum) | 0 | -7.86e28 (per cycle ~1e27-6e28; rel 4e-10 max, 3.7e-13 accumulated) | 1.49e28 (rel 1.6e-10 / 6.9e-14) |

**The one CPU difference:** the E fix segment is not exactly 0 on CPU but O(1e28) per cycle (round-off level, 1e-10 relative per cycle), in both arms.
L_z and P_r fix are exactly 0.

**Is D6 worse than PROD in any column?** Yes, in a few, by the same margins viper sees on GPU:
- accumulated |dLz| 4.775e40 vs 4.503e40 (+6 %, = GPU);
- REL Lz_tot accumulated 0.110 vs 0.071; REL E_a accumulated 0.404 vs 0.367; E_fix 3.7e-13 vs 6.9e-14 (both negligible); RM 1.8e15 vs 6.5e14 (rel 1e-26).

D6 is better in: per-cycle REL E_rsolve max 2.33 vs 3.44, Lz_tot per-cycle max 0.67 vs 1.50, Lz_solve per-cycle max 1.05 vs 1.73 and accumulated
0.165 vs 0.188; |solLz| 7.2e40 vs 1.19e41; hydLz 2.4e40 vs 7.3e40. RE_b and the P_r columns are identical in both arms.

### 10-cycle SUM / REL
```
--- A_D6
SUM  -2.224e+40   1.015e+35   1.517e+40         nan  -7.065e+39         nan  -5.286e+39  -7.942e+28   2.520e+40  -2.698e+40   0.000e+00   0.000e+00   0.000e+00  -6.017e+40   9.636e+29  -6.017e+40  -9.035e+40  -3.830e+38  -8.996e+40   3.018e+40   0.000e+00   3.043e+34  -1.276e+37  -1.278e+37   4.748e+34   6.299e+29   4.748e+34   0.000e+00   1.416e+15
REL E_a       per-cycle max 1.575e+00  accumulated 3.058e-01
REL E_b       per-cycle max nan  accumulated nan
REL E_fix     per-cycle max 4.149e-10  accumulated 3.437e-12
REL E_rsolve  per-cycle max 1.104e+00  accumulated 2.288e-01
REL Lz_tot    per-cycle max 6.721e-01  accumulated 4.152e-01
REL Lz_fix    per-cycle max 0.000e+00  accumulated 0.000e+00
REL Lz_solve  per-cycle max 1.006e+00  accumulated 6.238e-01
REL Pr_fix    per-cycle max 0.000e+00  accumulated 0.000e+00
REL Pr_solve  per-cycle max 3.745e-03  accumulated 3.721e-03
REL M         per-cycle max 6.456e-24  accumulated 6.126e-26
--- A_PROD
SUM  -2.297e+40   1.017e+35   1.674e+40         nan  -6.227e+39         nan  -5.629e+39  -4.333e+27   2.555e+40  -2.615e+40   0.000e+00   0.000e+00   0.000e+00  -6.878e+40   8.204e+29  -6.878e+40  -1.731e+41   1.489e+37  -1.731e+41   1.043e+41   0.000e+00   3.043e+34  -1.276e+37  -1.278e+37   4.748e+34   6.430e+29   4.748e+34   0.000e+00   3.803e+14
REL E_a       per-cycle max 1.597e+00  accumulated 2.613e-01
REL E_b       per-cycle max nan  accumulated nan
REL E_fix     per-cycle max 1.632e-10  accumulated 1.818e-13
REL E_rsolve  per-cycle max 1.127e+00  accumulated 2.362e-01
REL Lz_tot    per-cycle max 1.502e+00  accumulated 3.698e-01
REL Lz_fix    per-cycle max 0.000e+00  accumulated 0.000e+00
REL Lz_solve  per-cycle max 1.547e+00  accumulated 9.331e-01
REL Pr_fix    per-cycle max 0.000e+00  accumulated 0.000e+00
REL Pr_solve  per-cycle max 3.745e-03  accumulated 3.721e-03
REL M         per-cycle max 3.486e-24  accumulated 1.596e-26
```

### cons_analyze.py A_D6/run.log (100 cycles, whole)
```
file A_D6/run.log cycles 1 .. 99 n 99
cy          dU         BhE         BrA         BrB        RE_a        RE_b          RS        fixE         s0E          T2           W         Pun           D         dLz        BhLz         RLz       solLz        LzG3        LzR3       hydLz       fixLz         dPr        BhPr       hydPr        fgas        PrR3       solPr       fixPr          RM
1   2.825e+38   1.055e+34  -1.442e+37         nan   2.681e+38         nan  -5.964e+37  -3.629e+28   2.126e+38   1.151e+38   0.000e+00   0.000e+00   0.000e+00  -1.238e+40   1.875e+28  -1.238e+40  -1.867e+40   2.771e+37  -1.869e+40   6.285e+39   0.000e+00   3.145e+33  -1.276e+36  -1.278e+36   4.779e+33   1.287e+29   4.779e+33   0.000e+00  -1.481e+14
2   1.528e+38   1.010e+34   8.787e+37   6.978e+37   2.407e+38   2.226e+38  -1.687e+38  -6.341e+28   5.072e+38  -9.781e+37   0.000e+00   0.000e+00   0.000e+00  -8.494e+39   6.194e+28  -8.494e+39  -1.263e+40   6.913e+36  -1.264e+40   4.138e+39   0.000e+00   3.120e+33  -1.276e+36  -1.278e+36   4.776e+33   1.506e+29   4.776e+33   0.000e+00   9.865e+14
3  -2.768e+38   9.887e+33   4.479e+38   4.787e+38   1.711e+38   2.019e+38  -2.757e+38   2.159e+28   9.179e+38  -4.710e+38   0.000e+00   0.000e+00   0.000e+00  -9.705e+39  -8.983e+28  -9.705e+39  -1.450e+40  -1.259e+37  -1.449e+40   4.794e+39   0.000e+00   3.094e+33  -1.276e+36  -1.278e+36   4.771e+33   1.538e+29   4.771e+33   0.000e+00   1.722e+14
4  -1.270e+39   9.627e+33   1.209e+39   1.453e+39  -6.017e+37   1.832e+38  -3.807e+38   1.006e+28   1.594e+39  -1.274e+39   0.000e+00   0.000e+00   0.000e+00  -9.243e+39  -4.144e+29  -9.243e+39  -1.387e+40  -2.511e+37  -1.385e+40   4.628e+39   0.000e+00   3.070e+33  -1.276e+36  -1.278e+36   4.764e+33   1.190e+29   4.764e+33   0.000e+00   2.915e+14
5  -2.608e+39   9.740e+33   2.020e+39   2.774e+39  -5.879e+38   1.661e+38  -4.840e+38   7.157e+27   2.519e+39  -2.623e+39   0.000e+00   0.000e+00   0.000e+00  -9.596e+39  -5.271e+28  -9.596e+39  -1.441e+40  -3.848e+37  -1.437e+40   4.811e+39   0.000e+00   3.048e+33  -1.276e+36  -1.278e+36   4.757e+33   6.017e+28   4.757e+33   0.000e+00  -4.588e+14
6  -3.409e+39   1.028e+34   2.207e+39   3.560e+39  -1.203e+39   1.507e+38  -5.857e+38  -7.814e+27   3.260e+39  -3.876e+39   0.000e+00   0.000e+00   0.000e+00  -9.693e+39   2.511e+29  -9.693e+39  -1.452e+40  -5.100e+37  -1.447e+40   4.831e+39   0.000e+00   3.029e+33  -1.276e+36  -1.278e+36   4.748e+33   2.241e+28   4.748e+33   0.000e+00   1.693e+14
7  -3.489e+39   1.039e+34   2.094e+39   3.626e+39  -1.395e+39   1.367e+38  -6.859e+38   5.223e+27   3.560e+39  -4.269e+39   0.000e+00   0.000e+00   0.000e+00  -1.056e+40  -2.009e+29  -1.056e+40  -1.584e+40  -6.366e+37  -1.578e+40   5.286e+39   0.000e+00   3.009e+33  -1.276e+36  -1.278e+36   4.738e+33   1.784e+28   4.738e+33   0.000e+00  -4.674e+13
8  -3.685e+39   1.016e+34   2.339e+39   3.809e+39  -1.346e+39   1.240e+38  -7.847e+38  -4.836e+27   3.853e+39  -4.415e+39   0.000e+00   0.000e+00   0.000e+00  -8.770e+39   7.039e+29  -8.770e+39  -1.318e+40  -7.640e+37  -1.311e+40   4.414e+39   0.000e+00   2.991e+33  -1.276e+36  -1.278e+36   4.727e+33   5.358e+27   4.727e+33   0.000e+00   1.653e+14
9  -3.985e+39   1.035e+34   2.460e+39   4.098e+39  -1.525e+39   1.125e+38  -8.822e+38  -5.919e+27   4.270e+39  -4.913e+39   0.000e+00   0.000e+00   0.000e+00   1.445e+40  -3.106e+29   1.445e+40   2.158e+40  -7.726e+37   2.165e+40  -7.122e+39   0.000e+00   2.972e+33  -1.276e+36  -1.278e+36   4.715e+33  -1.334e+28   4.715e+33   0.000e+00   1.886e+14
10  -3.949e+39   1.046e+34   2.320e+39   4.051e+39  -1.629e+39   1.021e+38  -9.785e+38  -5.184e+27   4.509e+39  -5.159e+39   0.000e+00   0.000e+00   0.000e+00   3.819e+39   9.963e+29   3.819e+39   5.702e+39  -7.313e+37   5.776e+39  -1.883e+39   0.000e+00   2.953e+33  -1.276e+36  -1.278e+36   4.702e+33  -1.458e+28   4.702e+33   0.000e+00   9.582e+13
97  -1.259e+39   1.018e+34   7.334e+38   1.259e+39  -5.253e+38   1.147e+34  -2.901e+39   4.101e+27   8.081e+39  -5.705e+39   0.000e+00   0.000e+00   0.000e+00   2.951e+39   2.128e+30   2.951e+39   4.434e+39  -4.846e+36   4.439e+39  -1.483e+39   0.000e+00  -1.451e+32  -1.274e+36  -1.277e+36   2.497e+33  -2.430e+27   2.497e+33   0.000e+00   1.908e+13
98  -1.245e+39   1.018e+34   7.257e+38   1.245e+39  -5.197e+38   9.059e+33  -2.885e+39  -1.702e+27   8.035e+39  -5.669e+39   0.000e+00   0.000e+00   0.000e+00  -5.566e+38   1.728e+30  -5.566e+38  -8.346e+38  -4.237e+36  -8.303e+38   2.779e+38   0.000e+00  -1.476e+32  -1.274e+36  -1.277e+36   2.494e+33  -2.424e+27   2.494e+33   0.000e+00   8.422e+13
99  -1.232e+39   1.018e+34   7.183e+38   1.232e+39  -5.142e+38   6.942e+33  -2.870e+39   2.824e+27   7.989e+39  -5.633e+39   0.000e+00   0.000e+00   0.000e+00  -5.527e+38   1.586e+30  -5.527e+38  -8.296e+38  -5.929e+36  -8.237e+38   2.769e+38   0.000e+00  -1.499e+32  -1.274e+36  -1.277e+36   2.492e+33  -2.392e+27   2.492e+33   0.000e+00   2.512e+13
SUM  -2.139e+41   1.014e+36   1.272e+41         nan  -8.667e+40         nan  -2.633e+41  -7.857e+28   7.808e+41  -6.041e+41   0.000e+00   0.000e+00   0.000e+00  -4.775e+40   9.827e+31  -4.775e+40  -7.166e+40   6.430e+38  -7.230e+40   2.391e+40   0.000e+00   1.174e+35  -1.262e+38  -1.264e+38   3.425e+35   1.488e+29   3.425e+35   0.000e+00   1.768e+15
ABS   2.148e+41   1.014e+36   1.273e+41         nan   8.803e+40         nan   2.633e+41   4.195e+29   7.808e+41   6.044e+41   0.000e+00   0.000e+00   0.000e+00   2.891e+41   1.216e+32   2.891e+41   4.335e+41   3.236e+39   4.353e+41   1.444e+41   0.000e+00   1.232e+35   1.262e+38   1.264e+38   3.425e+35   1.167e+30   3.425e+35   0.000e+00   8.062e+15
REL E_a       per-cycle max 1.575e+00  accumulated 4.035e-01
REL E_b       per-cycle max nan  accumulated nan
REL E_fix     per-cycle max 4.149e-10  accumulated 3.658e-13
REL E_rsolve  per-cycle max 2.329e+00  accumulated 1.226e+00
REL Lz_tot    per-cycle max 6.738e-01  accumulated 1.097e-01
REL Lz_fix    per-cycle max 0.000e+00  accumulated 0.000e+00
REL Lz_solve  per-cycle max 1.054e+00  accumulated 1.646e-01
REL Pr_fix    per-cycle max 0.000e+00  accumulated 0.000e+00
REL Pr_solve  per-cycle max 3.745e-03  accumulated 2.714e-03
REL M         per-cycle max 6.456e-24  accumulated 8.231e-27
BIN solve U  -3.43e+41 -2.65e+41 -8.81e+40  1.10e+40  8.54e+39  2.83e+39  6.94e+37  6.23e+37  1.53e+37 -1.60e+41 -1.24e+41 -4.13e+40 -1.81e+38 -1.59e+38 -4.55e+37 -7.69e+37 -7.48e+37 -1.83e+37  4.85e+37  4.52e+37  1.04e+37  2.51e+39  2.23e+39  5.91e+38
BIN solve Lz -3.79e+32  4.20e+32 -3.98e+31  4.69e+29 -4.25e+31  4.22e+31  4.09e+39 -1.03e+40  6.20e+39  2.58e+40 -6.34e+40  3.79e+40 -1.71e+39 -2.97e+38  2.27e+38 -4.27e+39 -6.73e+38  4.84e+38 -7.56e+39 -8.06e+38  6.16e+38 -6.72e+40  4.90e+39  4.35e+39
BIN solve Pr  1.68e+35  1.30e+35  4.33e+34  1.68e+33  1.30e+33  4.32e+32  6.79e+31  5.25e+31  1.75e+31 -1.41e+33 -1.09e+33 -3.62e+32 -1.65e+28 -1.27e+28 -4.26e+27 -1.67e+28 -1.24e+28 -4.39e+27 -1.35e+28 -9.33e+27 -3.80e+27 -6.03e+28 -3.87e+28 -1.68e+28
BIN fix   U  -3.78e+28 -2.85e+28 -1.23e+28  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00
BIN fix   Lz  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00
BIN fix   Pr  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00
```

### cons_analyze.py A_PROD/run.log (100 cycles, whole)
```
file A_PROD/run.log cycles 1 .. 99 n 99
cy          dU         BhE         BrA         BrB        RE_a        RE_b          RS        fixE         s0E          T2           W         Pun           D         dLz        BhLz         RLz       solLz        LzG3        LzR3       hydLz       fixLz         dPr        BhPr       hydPr        fgas        PrR3       solPr       fixPr          RM
1   2.806e+38   1.055e+34  -1.218e+37         nan   2.684e+38         nan  -5.964e+37   4.580e+28   2.132e+38   1.149e+38   0.000e+00   0.000e+00   0.000e+00   1.459e+38  -4.657e+28   1.459e+38   9.714e+37   3.435e+37   6.280e+37   4.875e+37   0.000e+00   3.145e+33  -1.276e+36  -1.278e+36   4.779e+33   1.329e+29   4.779e+33   0.000e+00   6.967e+14
2   1.498e+38   1.011e+34   8.947e+37   7.280e+37   2.393e+38   2.226e+38  -1.687e+38   1.045e+28   5.089e+38  -1.009e+38   0.000e+00   0.000e+00   0.000e+00   1.795e+39  -3.569e+28   1.795e+39   2.755e+39   3.275e+37   2.723e+39  -9.601e+38   0.000e+00   3.120e+33  -1.276e+36  -1.278e+36   4.776e+33   1.557e+29   4.776e+33   0.000e+00  -5.222e+14
3  -2.752e+38   9.913e+33   4.447e+38   4.771e+38   1.695e+38   2.019e+38  -2.757e+38  -1.164e+28   9.177e+38  -4.724e+38   0.000e+00   0.000e+00   0.000e+00   1.020e+39   2.625e+29   1.020e+39   1.573e+39   3.434e+37   1.539e+39  -5.536e+38   0.000e+00   3.094e+33  -1.276e+36  -1.278e+36   4.771e+33   1.618e+29   4.771e+33   0.000e+00  -3.397e+14
4  -1.273e+39   9.639e+33   1.216e+39   1.457e+39  -5.746e+37   1.832e+38  -3.807e+38  -1.621e+28   1.594e+39  -1.270e+39   0.000e+00   0.000e+00   0.000e+00   6.845e+38   1.471e+29   6.845e+38   1.063e+39   3.525e+37   1.028e+39  -3.783e+38   0.000e+00   3.070e+33  -1.276e+36  -1.278e+36   4.764e+33   1.271e+29   4.764e+33   0.000e+00   3.822e+14
5  -2.659e+39   9.735e+33   2.070e+39   2.825e+39  -5.893e+38   1.661e+38  -4.840e+38  -2.360e+27   2.536e+39  -2.641e+39   0.000e+00   0.000e+00   0.000e+00   3.470e+38   4.809e+29   3.470e+38   5.565e+38   3.515e+37   5.213e+38  -2.095e+38   0.000e+00   3.048e+33  -1.276e+36  -1.278e+36   4.757e+33   6.338e+28   4.757e+33   0.000e+00   2.221e+14
6  -3.545e+39   1.028e+34   2.313e+39   3.695e+39  -1.231e+39   1.507e+38  -5.857e+38  -1.145e+28   3.319e+39  -3.964e+39   0.000e+00   0.000e+00   0.000e+00   1.525e+38   4.311e+29   1.525e+38   2.443e+38   3.491e+37   2.093e+38  -9.172e+37   0.000e+00   3.029e+33  -1.276e+36  -1.278e+36   4.748e+33   1.945e+28   4.748e+33   0.000e+00  -5.291e+13
7  -3.674e+39   1.044e+34   2.208e+39   3.811e+39  -1.466e+39   1.367e+38  -6.859e+38  -9.749e+27   3.653e+39  -4.433e+39   0.000e+00   0.000e+00   0.000e+00   9.816e+37   4.116e+29   9.816e+37   1.520e+38   3.488e+37   1.171e+38  -5.384e+37   0.000e+00   3.009e+33  -1.276e+36  -1.278e+36   4.738e+33   1.174e+28   4.738e+33   0.000e+00   2.178e+13
8  -3.829e+39   1.020e+34   3.521e+39         nan  -3.086e+38         nan  -1.128e+39  -6.577e+26   3.924e+39  -3.105e+39   0.000e+00   0.000e+00   0.000e+00  -7.113e+40  -3.561e+29  -7.113e+40  -1.766e+41  -8.561e+37  -1.765e+41   1.054e+41   0.000e+00   2.991e+33  -1.276e+36  -1.278e+36   4.727e+33   3.290e+27   4.727e+33   0.000e+00  -2.861e+13
9  -4.099e+39   1.035e+34   2.521e+39   4.212e+39  -1.579e+39   1.125e+38  -8.822e+38  -1.385e+28   4.328e+39  -5.025e+39   0.000e+00   0.000e+00   0.000e+00  -1.397e+39   8.590e+28  -1.397e+39  -2.169e+39  -7.077e+37  -2.098e+39   7.718e+38   0.000e+00   2.972e+33  -1.276e+36  -1.278e+36   4.715e+33  -1.578e+28   4.715e+33   0.000e+00   1.310e+14
10  -4.047e+39   1.047e+34   2.374e+39   4.149e+39  -1.673e+39   1.021e+38  -9.785e+38   5.339e+27   4.560e+39  -5.255e+39   0.000e+00   0.000e+00   0.000e+00  -5.032e+38  -5.603e+29  -5.032e+38  -8.288e+38  -7.036e+37  -7.584e+38   3.256e+38   0.000e+00   2.953e+33  -1.276e+36  -1.278e+36   4.702e+33  -1.650e+28   4.702e+33   0.000e+00  -1.301e+14
97  -1.258e+39   1.018e+34   7.329e+38   1.258e+39  -5.253e+38   1.147e+34  -2.901e+39  -4.990e+27   8.081e+39  -5.705e+39   0.000e+00   0.000e+00   0.000e+00   3.560e+37   2.783e+30   3.560e+37   5.913e+37  -3.964e+36   6.310e+37  -2.354e+37   0.000e+00  -1.451e+32  -1.274e+36  -1.277e+36   2.497e+33  -2.441e+27   2.497e+33   0.000e+00  -3.492e+13
98  -1.244e+39   1.018e+34   7.247e+38   1.244e+39  -5.194e+38   9.060e+33  -2.885e+39   1.857e+27   8.034e+39  -5.668e+39   0.000e+00   0.000e+00   0.000e+00   4.621e+37   2.430e+30   4.621e+37   6.967e+37  -3.781e+36   7.345e+37  -2.346e+37   0.000e+00  -1.476e+32  -1.274e+36  -1.277e+36   2.494e+33  -2.373e+27   2.494e+33   0.000e+00   4.908e+08
99  -1.231e+39   1.018e+34   7.177e+38   1.231e+39  -5.136e+38   6.942e+33  -2.870e+39  -3.056e+27   7.989e+39  -5.632e+39   0.000e+00   0.000e+00   0.000e+00   5.674e+37   1.812e+30   5.674e+37   8.512e+37  -3.726e+36   8.885e+37  -2.838e+37   0.000e+00  -1.500e+32  -1.274e+36  -1.277e+36   2.492e+33  -2.345e+27   2.492e+33   0.000e+00   6.820e+13
SUM  -2.156e+41   1.015e+36   1.363e+41         nan  -7.933e+40         nan  -2.797e+41   1.486e+28   7.816e+41  -5.813e+41   0.000e+00   0.000e+00   0.000e+00  -4.503e+40   1.089e+32  -4.503e+40  -1.185e+41   1.295e+39  -1.198e+41   7.349e+40   0.000e+00   1.174e+35  -1.262e+38  -1.264e+38   3.425e+35   1.508e+29   3.425e+35   0.000e+00   6.477e+14
ABS   2.165e+41   1.015e+36   1.363e+41         nan   8.068e+40         nan   2.797e+41   4.251e+29   7.816e+41   5.815e+41   0.000e+00   0.000e+00   0.000e+00   2.579e+41   1.497e+32   2.579e+41   6.292e+41   3.853e+39   6.299e+41   3.714e+41   0.000e+00   1.232e+35   1.262e+38   1.264e+38   3.425e+35   1.200e+30   3.425e+35   0.000e+00   7.729e+15
REL E_a       per-cycle max 1.597e+00  accumulated 3.665e-01
REL E_b       per-cycle max nan  accumulated nan
REL E_fix     per-cycle max 1.632e-10  accumulated 6.862e-14
REL E_rsolve  per-cycle max 3.444e+00  accumulated 1.292e+00
REL Lz_tot    per-cycle max 1.502e+00  accumulated 7.149e-02
REL Lz_fix    per-cycle max 0.000e+00  accumulated 0.000e+00
REL Lz_solve  per-cycle max 1.726e+00  accumulated 1.881e-01
REL Pr_fix    per-cycle max 0.000e+00  accumulated 0.000e+00
REL Pr_solve  per-cycle max 3.745e-03  accumulated 2.714e-03
REL M         per-cycle max 3.486e-24  accumulated 2.992e-27
BIN solve U  -3.43e+41 -2.65e+41 -8.81e+40  1.10e+40  8.54e+39  2.83e+39  6.94e+37  6.23e+37  1.53e+37 -1.62e+41 -1.25e+41 -4.16e+40 -1.67e+38 -1.66e+38 -4.21e+37 -7.16e+37 -8.31e+37 -1.49e+37  3.95e+37  4.35e+37  9.95e+36  2.46e+39  2.34e+39  5.34e+38
BIN solve Lz  1.36e+32 -1.57e+32  2.20e+31  1.70e+32 -1.47e+32 -2.29e+31  4.09e+39 -1.03e+40  6.20e+39  2.60e+40 -6.32e+40  3.80e+40 -2.84e+39 -4.49e+38  4.07e+38 -6.89e+39 -1.17e+39  9.08e+38 -1.20e+40 -1.38e+39  1.18e+39 -1.12e+41  7.22e+39  7.75e+39
BIN solve Pr  1.68e+35  1.30e+35  4.33e+34  1.68e+33  1.30e+33  4.32e+32  6.79e+31  5.25e+31  1.75e+31 -1.41e+33 -1.09e+33 -3.62e+32 -1.65e+28 -1.27e+28 -4.29e+27 -1.70e+28 -1.21e+28 -4.54e+27 -1.41e+28 -8.53e+27 -4.14e+27 -6.05e+28 -3.51e+28 -1.79e+28
BIN fix   U   6.69e+27  3.06e+27  5.11e+27  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00
BIN fix   Lz  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00
BIN fix   Pr  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00  0.00e+00
```

## Neutrality on CPU: **PASS** (A keys on = N keys off = G plain build, first 41 history rows)
Note: the history files have 2 header lines + 41 rows (cycles 0-40) = 43 lines, so the TASK's `head -n 42` cuts one row and gives a spurious DIFF.
With `head -n 43`:

| arm | file | A (first 43 lines) | N | G |
|---|---|---|---|---|
| D6 | hydro | 0b2bef81 | 0b2bef81 | 0b2bef81 |
| D6 | user | e2be5729 | e2be5729 | e2be5729 |
| PROD | hydro | 48502d66 | 48502d66 | 48502d66 |
| PROD | user | e458e658 | e458e658 | e458e658 |

All identical. G (plain build + AM print keys) = N, as stated in part (b).
