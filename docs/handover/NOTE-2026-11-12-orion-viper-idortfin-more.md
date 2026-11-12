# NOTE 2026-11-12 Orion -> viper: idortfin-more (posted per part; this is part 3; parts 2 and 1 follow)

Binary for all parts: 421a3d024 (MPI ON, gcc/13 + openmpi/4.1, Release, -ffp-contract=off; md5 e151dc3ba4a930d2d348fee069aa50a0;
/orion/ptmp/jinma/rwfin_1112/src_421a3d024/build_none/src/athena, the idortfin-rw binary). OMP 1, OMPI_MCA_hwloc_base_binding_policy=none,
`mpirun --bind-to none`. Each part ran on its own p.exclusive node, all runs of a part at once.

## (3) pulse_cs ladders, FIN vs B: **FIN == B BITWISE**
All 14 runs: rc=0, 0 FATAL, Picard NON-CONVERGED=0. Host orion072 (T, 48 cores) and orion075 (C, 72 cores);
total wall 70 s (T) / 76 s (C).

| run | rc | wall s | last "Picard iterations" (solves, mean, max, NC) |
|---|---|---|---|
| B_C_128 | 0 | 75 | 2.850000e+02, 2.010526e+00, 3.000000e+00, 0.000000e+00 |
| B_C_32 | 0 | 10 | 7.100000e+01, 2.028169e+00, 3.000000e+00, 0.000000e+00 |
| B_C_64 | 0 | 26 | 1.430000e+02, 2.006993e+00, 3.000000e+00, 0.000000e+00 |
| B_T_1 | 0 | 17 | 3.500000e+01, 2.114286e+00, 4.000000e+00, 0.000000e+00 |
| B_T_2 | 0 | 25 | 7.100000e+01, 2.042254e+00, 3.000000e+00, 0.000000e+00 |
| B_T_4 | 0 | 39 | 1.430000e+02, 2.020979e+00, 3.000000e+00, 0.000000e+00 |
| B_T_8 | 0 | 67 | 2.850000e+02, 2.003509e+00, 3.000000e+00, 0.000000e+00 |
| FIN_C_128 | 0 | 75 | 2.850000e+02, 2.010526e+00, 3.000000e+00, 0.000000e+00 |
| FIN_C_32 | 0 | 10 | 7.100000e+01, 2.028169e+00, 3.000000e+00, 0.000000e+00 |
| FIN_C_64 | 0 | 28 | 1.430000e+02, 2.006993e+00, 3.000000e+00, 0.000000e+00 |
| FIN_T_1 | 0 | 17 | 3.500000e+01, 2.114286e+00, 4.000000e+00, 0.000000e+00 |
| FIN_T_2 | 0 | 25 | 7.100000e+01, 2.042254e+00, 3.000000e+00, 0.000000e+00 |
| FIN_T_4 | 0 | 40 | 1.430000e+02, 2.020979e+00, 3.000000e+00, 0.000000e+00 |
| FIN_T_8 | 0 | 70 | 2.850000e+02, 2.003509e+00, 3.000000e+00, 0.000000e+00 |

pulse_ord.py (T ladder):
```
B    m1_e   d 2.615e-02 8.041e-03 2.411e-03  p 1.70 1.74
B    m1_f1  d 2.026e-02 6.306e-03 1.765e-03  p 1.68 1.84
B    m1_f2  d 1.278e-01 5.316e-02 2.867e-02  p 1.27 0.89
B    m1_f3  d 1.317e-01 5.421e-02 3.011e-02  p 1.28 0.85
FIN  m1_e   d 2.615e-02 8.041e-03 2.411e-03  p 1.70 1.74
FIN  m1_f1  d 2.026e-02 6.306e-03 1.765e-03  p 1.68 1.84
FIN  m1_f2  d 1.278e-01 5.316e-02 2.867e-02  p 1.27 0.89
FIN  m1_f3  d 1.317e-01 5.421e-02 3.011e-02  p 1.28 0.85
```

- **The FIN fields equal B bitwise:** max |B - FIN| = 0 in m1_e, m1_f1..f3 for T_8, T_1 and C_128 (bin_convert data compare); the .hst files are byte-identical.
  This is expected, not a dead key:
  - this deck has kappa_p = kappa_f = 0 (scattering only) and implicit_vimp = false, so wbsolve / am_cs_fold (gas write-back) have nothing to act on;
  - the FIN log shows face_mmat limited x1 faces 0 (of 571 checked passes), vet_col_lat_posfix=lag triggers 0, and no stage fallbacks (stage_local_be idle);
  - no FIN key is in the "never read" list, and the FIN banners print.
- **FLAG rule (FIN T-order < 1.8 where B >= 1.9 on m1_e / m1_f1): not triggered,** since FIN = B.
  But note that **B itself is below 1.9 here:** m1_e p 1.70 / 1.74, m1_f1 1.68 / 1.84; the lateral m1_f2/f3 are 1.27 / 0.89 and 1.28 / 0.85.
- C ladder: rc 0 for all 6; walls C_32 10 s, C_64 26-28 s, C_128 75 s (24 ranks).
