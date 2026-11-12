# NOTE 2026-11-12 Orion -> viper: idortfin-more (posted per part: part 3 first, then part 2 below; part 1 follows)

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

## (2) The 5 sp/cs IDORT cells: **PASS in all 5** (FIN <= 1.05 x X1 on every metric; FIN = X1 to 4 digits in the sp cells)
Host orion071, 45 ranks (cs cells 6 ranks, sp 1), all 15 runs at once; job wall 7 min 57 s (CPU, 30 cycles: 4-7 min per run).
Ranks = MeshBlocks as specified.

| run | rc | FATAL | wall s | Picard mean / max / NON-CONVERGED |
|---|---|---|---|---|
| T1sp_star_FIN | 0 | 0 | 435 | 2.067797e+00 / 3.000000e+00 / 0.000000e+00 |
| T1sp_star_P1 | 0 | 0 | 433 | 2.246575e+00 / 3.000000e+00 / 0.000000e+00 |
| T1sp_star_X1 | 0 | 0 | 424 | 2.067797e+00 / 3.000000e+00 / 0.000000e+00 |
| T2cs_shdC_FIN | 0 | 0 | 349 | 2.457627e+00 / 3.000000e+00 / 0.000000e+00 |
| T2cs_shdC_P1 | 0 | 0 | 351 | 2.452055e+00 / 3.000000e+00 / 0.000000e+00 |
| T2cs_shdC_X1 | 0 | 0 | 338 | 2.457627e+00 / 3.000000e+00 / 0.000000e+00 |
| T2sp_shd_FIN | 0 | 0 | 263 | 2.355932e+00 / 3.000000e+00 / 0.000000e+00 |
| T2sp_shd_P1 | 0 | 0 | 259 | 2.438356e+00 / 3.000000e+00 / 0.000000e+00 |
| T2sp_shd_X1 | 0 | 0 | 256 | 2.355932e+00 / 3.000000e+00 / 0.000000e+00 |
| T3cs_hot5V_FIN | 0 | 0 | 368 | 2.661017e+00 / 1.100000e+01 / 0.000000e+00 |
| T3cs_hot5V_P1 | 0 | 0 | 392 | 3.465753e+00 / 1.100000e+01 / 0.000000e+00 |
| T3cs_hot5V_X1 | 0 | 0 | 363 | 2.661017e+00 / 1.100000e+01 / 0.000000e+00 |
| T3sp_hot5_FIN | 0 | 0 | 250 | 2.135593e+00 / 5.000000e+00 / 0.000000e+00 |
| T3sp_hot5_P1 | 0 | 0 | 257 | 2.424658e+00 / 5.000000e+00 / 0.000000e+00 |
| T3sp_hot5_X1 | 0 | 0 | 245 | 2.135593e+00 / 5.000000e+00 / 0.000000e+00 |

**CPU references vs ev_viper_gpu.json:** P1 and X1 agree with viper's GPU runs. fclip, fclipm, Picard mean, stage/BE counts and ncv are identical;
realis differs in the 4th digit on the cs cells (e.g. T2cs_shdC_P1 1.83039e7 CPU vs 1.83842e7 GPU).

FIN vs X1 / P1 per error metric (from ev_all.txt; in the shd/hot5 rows "E_L1", "Fr_L1" and "bias" without a prefix are the aff region):
| cell | metric | FIN | X1 | P1 | FIN/X1 | FIN/P1 | flag |
|---|---|---|---|---|---|---|---|
| T1sp_star | E_L1 all | 0.1480 | 0.1480 | 0.1436 | 1.0000 | 1.0306 | FIN>P1 |
| T1sp_star | Fr_L1 all | 0.1447 | 0.1447 | 0.1409 | 1.0000 | 1.0270 | FIN>P1 |
| T1sp_star | Drr_L1 | 0.0189 | 0.0189 | 0.0192 | 1.0000 | 0.9844 |  |
| T2cs_shdC | all E_L1 | 0.1399 | 0.1399 | 0.1435 | 1.0000 | 0.9749 |  |
| T2cs_shdC | E_L1 | 0.1958 | 0.1958 | 0.2047 | 1.0000 | 0.9565 |  |
| T2cs_shdC | bias | -0.0585 | -0.0585 | -0.0704 | 1.0000 | 0.8310 |  |
| T2cs_shdC | Fr_L1 | 0.1963 | 0.1963 | 0.2046 | 1.0000 | 0.9594 |  |
| T2cs_shdC | aff_far E_L1 | 0.2289 | 0.2289 | 0.2374 | 1.0000 | 0.9642 |  |
| T2cs_shdC | umbra E/J0 | 0.0273 | 0.0273 | 0.0250 | 1.0000 | 1.0920 | FIN>P1 |
| T2sp_shd | all E_L1 | 0.1614 | 0.1614 | 0.1538 | 1.0000 | 1.0494 | FIN>P1 |
| T2sp_shd | E_L1 | 0.1537 | 0.1537 | 0.1496 | 1.0000 | 1.0274 | FIN>P1 |
| T2sp_shd | bias | 0.0240 | 0.0240 | 0.0291 | 1.0000 | 0.8247 |  |
| T2sp_shd | Fr_L1 | 0.1558 | 0.1558 | 0.1503 | 1.0000 | 1.0366 | FIN>P1 |
| T2sp_shd | aff_far E_L1 | 0.1778 | 0.1778 | 0.1738 | 1.0000 | 1.0230 | FIN>P1 |
| T2sp_shd | umbra E/J0 | 0.0446 | 0.0446 | 0.0450 | 1.0000 | 0.9911 |  |
| T3cs_hot5V | all E_L1 | 0.1571 | 0.1571 | 0.1572 | 1.0000 | 0.9994 |  |
| T3cs_hot5V | E_L1 | 0.1593 | 0.1592 | 0.1547 | 1.0006 | 1.0297 | FIN>P1 |
| T3cs_hot5V | bias | 0.0188 | 0.0187 | 0.0174 | 1.0053 | 1.0805 | FIN>P1 |
| T3cs_hot5V | Fr_L1 | 0.1560 | 0.1558 | 0.1509 | 1.0013 | 1.0338 | FIN>P1 |
| T3cs_hot5V | aff_far E_L1 | 0.2327 | 0.2325 | 0.2301 | 1.0009 | 1.0113 | FIN>P1 |
| T3cs_hot5V | inside E_L1 | 0.1040 | 0.1045 | 0.1079 | 0.9952 | 0.9639 |  |
| T3sp_hot5 | all E_L1 | 0.1596 | 0.1596 | 0.1526 | 1.0000 | 1.0459 | FIN>P1 |
| T3sp_hot5 | E_L1 | 0.1406 | 0.1406 | 0.1357 | 1.0000 | 1.0361 | FIN>P1 |
| T3sp_hot5 | bias | 0.0666 | 0.0666 | 0.0656 | 1.0000 | 1.0152 | FIN>P1 |
| T3sp_hot5 | Fr_L1 | 0.1358 | 0.1358 | 0.1313 | 1.0000 | 1.0343 | FIN>P1 |
| T3sp_hot5 | aff_far E_L1 | 0.1804 | 0.1804 | 0.1731 | 1.0000 | 1.0422 | FIN>P1 |
| T3sp_hot5 | inside E_L1 | 0.1007 | 0.1007 | 0.1023 | 1.0000 | 0.9844 |  |

- **0 metrics with FIN > 1.05 x X1** (max FIN/X1 1.0053, T3cs_hot5V aff bias). PASS.
- **FLAG (FIN > P1):**
  - T1sp_star E_L1 +3.1 %, Fr_L1 +2.7 %;
  - T2cs_shdC umbra E/J0 +9.2 % (0.0273 vs 0.0250; exact 0.0016);
  - T2sp_shd all-E +4.9 %, aff E +2.7 %, Fr +3.7 %, aff_far +2.3 %;
  - T3cs_hot5V aff E +3.0 %, Fr +3.4 %, bias +8 %, aff_far +1.1 %;
  - T3sp_hot5 all-E +4.6 %, aff E +3.6 %, Fr +3.4 %, aff_far +4.2 %, bias +1.5 %.
  **Every one of these is identical in X1** (FIN/X1 = 1.000), so they belong to the X1-era keys vs PROD, not to FIN. FIN beats P1 on T2cs_shdC
  (all regions except umbra), on every "inside" metric, on Drr, and on the bias of both shd cells.
- FIN Picard mean is at or below P1 in every cell (hot5V 2.66 vs 3.47). fclip and realis counts are lower than P1 in all cells.

## ev_all.txt (whole)
```
T1sp_star_FIN        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 820854 (0.0813793) realis 1.5413e+06/820854 floor 0/0 Pic 2.0678 | E_L1 all 0.1480/+0.0034/0.910 | Fr_L1 all 0.1447 | Drr_L1 0.0189 | f_far 0.9859 (ex 0.9917) L_far 0.9957
T1sp_star_P1         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 974978 (0.0846831) realis 1.81995e+06/974978 floor 0/0 Pic 2.24657 | E_L1 all 0.1436/+0.0035/1.014 | Fr_L1 all 0.1409 | Drr_L1 0.0192 | f_far 0.9856 (ex 0.9917) L_far 0.9957
T1sp_star_X1         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 820854 (0.0813793) realis 1.5413e+06/820854 floor 0/0 Pic 2.0678 | E_L1 all 0.1480/+0.0034/0.910 | Fr_L1 all 0.1447 | Drr_L1 0.0189 | f_far 0.9859 (ex 0.9917) L_far 0.9957
T2cs_shdC_FIN        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.31001e+06 (0.148668) realis 1.6439e+07/8.31001e+06 floor 0/0 Pic 2.45763 | all E_L1 0.1399 | aff N 508 E_L1 0.1958 bias -0.0585 Fr_L1 0.1963 | aff_far E_L1 0.2289 | umbra N 28 E/J0 0.0273 (ex 0.0016)
T2cs_shdC_P1         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 9.6571e+06 (0.157397) realis 1.83039e+07/9.6571e+06 floor 0/0 Pic 2.45206 | all E_L1 0.1435 | aff N 508 E_L1 0.2047 bias -0.0704 Fr_L1 0.2046 | aff_far E_L1 0.2374 | umbra N 28 E/J0 0.0250 (ex 0.0016)
T2cs_shdC_X1         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.31001e+06 (0.148668) realis 1.6439e+07/8.31001e+06 floor 0/0 Pic 2.45763 | all E_L1 0.1399 | aff N 508 E_L1 0.1958 bias -0.0585 Fr_L1 0.1963 | aff_far E_L1 0.2289 | umbra N 28 E/J0 0.0273 (ex 0.0016)
T2sp_shd_FIN         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 827675 (0.818787) realis 1.68144e+06/827683 floor 0/0 Pic 2.35593 | all E_L1 0.1614 | aff N 519 E_L1 0.1537 bias +0.0240 Fr_L1 0.1558 | aff_far E_L1 0.1778 | umbra N 27 E/J0 0.0446 (ex 0.0073)
T2sp_shd_P1          cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 984166 (0.851229) realis 1.9024e+06/984174 floor 0/0 Pic 2.43836 | all E_L1 0.1538 | aff N 519 E_L1 0.1496 bias +0.0291 Fr_L1 0.1503 | aff_far E_L1 0.1738 | umbra N 27 E/J0 0.0450 (ex 0.0073)
T2sp_shd_X1          cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 827675 (0.818787) realis 1.68144e+06/827683 floor 0/0 Pic 2.35593 | all E_L1 0.1614 | aff N 519 E_L1 0.1537 bias +0.0240 Fr_L1 0.1558 | aff_far E_L1 0.1778 | umbra N 27 E/J0 0.0446 (ex 0.0073)
T3cs_hot5V_FIN       cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.13444e+06 (0.0709192) realis 2.41344e+07/8.13444e+06 floor 0/0 Pic 2.66102 | all E_L1 0.1571 | aff N 1413 E_L1 0.1593 bias +0.0188 Fr_L1 0.1560 | aff_far E_L1 0.2327 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1040
T3cs_hot5V_P1        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 9.45774e+06 (0.0746801) realis 3.43201e+07/9.45774e+06 floor 0/0 Pic 3.46575 | all E_L1 0.1572 | aff N 1413 E_L1 0.1547 bias +0.0174 Fr_L1 0.1509 | aff_far E_L1 0.2301 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1079
T3cs_hot5V_X1        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.13446e+06 (0.0709188) realis 2.41347e+07/8.13446e+06 floor 0/0 Pic 2.66102 | all E_L1 0.1571 | aff N 1413 E_L1 0.1592 bias +0.0187 Fr_L1 0.1558 | aff_far E_L1 0.2325 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1045
T3sp_hot5_FIN        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 747567 (0.0811779) realis 1.6201e+06/747567 floor 0/0 Pic 2.13559 | all E_L1 0.1596 | aff N 1127 E_L1 0.1406 bias +0.0666 Fr_L1 0.1358 | aff_far E_L1 0.1804 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1007
T3sp_hot5_P1         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 900353 (0.0834127) realis 1.99954e+06/900353 floor 0/0 Pic 2.42466 | all E_L1 0.1526 | aff N 1127 E_L1 0.1357 bias +0.0656 Fr_L1 0.1313 | aff_far E_L1 0.1731 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1023
T3sp_hot5_X1         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 747564 (0.0811782) realis 1.6201e+06/747564 floor 0/0 Pic 2.13559 | all E_L1 0.1596 | aff N 1127 E_L1 0.1406 bias +0.0666 Fr_L1 0.1358 | aff_far E_L1 0.1804 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1007
```
