# NOTE 2026-11-12 Orion -> viper: idortd6 on the d6fast build (D6f), posted per part (this: part 2; part 1 follows)

**Binary (d6fast):** 5f06e171b236a2215b24fbcdfacab21b927039a5 (fork/d6fast-1112), built from a git archive with exactly the CMake options of the
FIN build (mpicxx gcc/13 openmpi/4.1.8, Athena_ENABLE_MPI=ON, Release, CMAKE_CXX_FLAGS=-ffp-contract=off, PROBLEM built_in_pgens);
**md5 cbb43fe701b280265096aba2ad708c40**. FIN references are my 421a3d024 runs (md5 e151dc3b...).
USER 11-12: run d6fast only. D6 on 421a3d024 was stopped. Only its 5 curv cells (and the N16 smoke) had already finished; they appear below as
"D6 (421a3d024)", and the D6 wave grid was not run. All D6f runs ran in ONE p.shared job (orion001, 28 cores; the exclusive queue was ~5 h):
the curv cells (15 ranks, --oversubscribe off-node not needed), vetf with --np 13 heaviest-first, then edd with --np 15. Binding off; run dirs tagged D6f.

## (2) The 5 sp/cs cells, D6f: **PASS in all 5** (D6f <= 1.05 x X1 on every metric; D6f = FIN = D6(421a3d024) to 4 digits)
| run | rc | FATAL | wall s | Picard mean / max / NC |
|---|---|---|---|---|
| T1sp_star_D6f | 0 | 0 | 588 | 2.067797e+00 / 3.000000e+00 / 0.000000e+00 |
| T2cs_shdC_D6f | 0 | 0 | 534 | 2.457627e+00 / 3.000000e+00 / 0.000000e+00 |
| T2sp_shd_D6f | 0 | 0 | 367 | 2.355932e+00 / 3.000000e+00 / 0.000000e+00 |
| T3cs_hot5V_D6f | 0 | 0 | 534 | 2.661017e+00 / 1.100000e+01 / 0.000000e+00 |
| T3sp_hot5_D6f | 0 | 0 | 363 | 2.135593e+00 / 5.000000e+00 / 0.000000e+00 |
(Wall is higher than the FIN runs because they shared 28 cores with the vetf list; the Picard counts equal FIN.)

| cell | metric | D6f | FIN | X1 | P1 | D6 (421a3d024) | D6f/X1 | D6f/FIN | D6f/P1 | flag |
|---|---|---|---|---|---|---|---|---|---|---|
| T1sp_star | E_L1 all | 0.1480 | 0.1480 | 0.1480 | 0.1436 | 0.1480 | 1.0000 | 1.0000 | 1.0306 | D6f>P1 |
| T1sp_star | Fr_L1 all | 0.1447 | 0.1447 | 0.1447 | 0.1409 | 0.1447 | 1.0000 | 1.0000 | 1.0270 | D6f>P1 |
| T1sp_star | Drr_L1 | 0.0189 | 0.0189 | 0.0189 | 0.0192 | 0.0189 | 1.0000 | 1.0000 | 0.9844 |  |
| T2cs_shdC | all E_L1 | 0.1399 | 0.1399 | 0.1399 | 0.1435 | 0.1399 | 1.0000 | 1.0000 | 0.9749 |  |
| T2cs_shdC | E_L1 | 0.1958 | 0.1958 | 0.1958 | 0.2047 | 0.1958 | 1.0000 | 1.0000 | 0.9565 |  |
| T2cs_shdC | bias | -0.0585 | -0.0585 | -0.0585 | -0.0704 | -0.0585 | 1.0000 | 1.0000 | 0.8310 |  |
| T2cs_shdC | Fr_L1 | 0.1963 | 0.1963 | 0.1963 | 0.2046 | 0.1963 | 1.0000 | 1.0000 | 0.9594 |  |
| T2cs_shdC | aff_far E_L1 | 0.2289 | 0.2289 | 0.2289 | 0.2374 | 0.2289 | 1.0000 | 1.0000 | 0.9642 |  |
| T2cs_shdC | umbra E/J0 | 0.0273 | 0.0273 | 0.0273 | 0.0250 | 0.0273 | 1.0000 | 1.0000 | 1.0920 | D6f>P1 |
| T2sp_shd | all E_L1 | 0.1614 | 0.1614 | 0.1614 | 0.1538 | 0.1614 | 1.0000 | 1.0000 | 1.0494 | D6f>P1 |
| T2sp_shd | E_L1 | 0.1537 | 0.1537 | 0.1537 | 0.1496 | 0.1537 | 1.0000 | 1.0000 | 1.0274 | D6f>P1 |
| T2sp_shd | bias | 0.0240 | 0.0240 | 0.0240 | 0.0291 | 0.0240 | 1.0000 | 1.0000 | 0.8247 |  |
| T2sp_shd | Fr_L1 | 0.1558 | 0.1558 | 0.1558 | 0.1503 | 0.1558 | 1.0000 | 1.0000 | 1.0366 | D6f>P1 |
| T2sp_shd | aff_far E_L1 | 0.1778 | 0.1778 | 0.1778 | 0.1738 | 0.1778 | 1.0000 | 1.0000 | 1.0230 | D6f>P1 |
| T2sp_shd | umbra E/J0 | 0.0446 | 0.0446 | 0.0446 | 0.0450 | 0.0446 | 1.0000 | 1.0000 | 0.9911 |  |
| T3cs_hot5V | all E_L1 | 0.1571 | 0.1571 | 0.1571 | 0.1572 | 0.1571 | 1.0000 | 1.0000 | 0.9994 |  |
| T3cs_hot5V | E_L1 | 0.1593 | 0.1593 | 0.1592 | 0.1547 | 0.1593 | 1.0006 | 1.0000 | 1.0297 | D6f>P1 |
| T3cs_hot5V | bias | 0.0188 | 0.0188 | 0.0187 | 0.0174 | 0.0188 | 1.0053 | 1.0000 | 1.0805 | D6f>P1 |
| T3cs_hot5V | Fr_L1 | 0.1560 | 0.1560 | 0.1558 | 0.1509 | 0.1560 | 1.0013 | 1.0000 | 1.0338 | D6f>P1 |
| T3cs_hot5V | aff_far E_L1 | 0.2327 | 0.2327 | 0.2325 | 0.2301 | 0.2327 | 1.0009 | 1.0000 | 1.0113 | D6f>P1 |
| T3cs_hot5V | inside E_L1 | 0.1040 | 0.1040 | 0.1045 | 0.1079 | 0.1040 | 0.9952 | 1.0000 | 0.9639 |  |
| T3sp_hot5 | all E_L1 | 0.1596 | 0.1596 | 0.1596 | 0.1526 | 0.1596 | 1.0000 | 1.0000 | 1.0459 | D6f>P1 |
| T3sp_hot5 | E_L1 | 0.1406 | 0.1406 | 0.1406 | 0.1357 | 0.1406 | 1.0000 | 1.0000 | 1.0361 | D6f>P1 |
| T3sp_hot5 | bias | 0.0666 | 0.0666 | 0.0666 | 0.0656 | 0.0666 | 1.0000 | 1.0000 | 1.0152 | D6f>P1 |
| T3sp_hot5 | Fr_L1 | 0.1358 | 0.1358 | 0.1358 | 0.1313 | 0.1358 | 1.0000 | 1.0000 | 1.0343 | D6f>P1 |
| T3sp_hot5 | aff_far E_L1 | 0.1804 | 0.1804 | 0.1804 | 0.1731 | 0.1804 | 1.0000 | 1.0000 | 1.0422 | D6f>P1 |
| T3sp_hot5 | inside E_L1 | 0.1008 | 0.1007 | 0.1007 | 0.1023 | 0.1008 | 1.0010 | 1.0010 | 0.9853 |  |

- **0 metrics with D6f > 1.05 x X1.** Max D6f/X1 is 1.0053 (T3cs_hot5V aff bias, the same as FIN). PASS.
- **D6f / FIN = 1.0000 on every metric** except T3sp_hot5 inside E_L1 (0.1008 vs 0.1007, +0.1 %). D6f = D6 (421a3d024) on every printed digit, so the
  d6fast shell kernel and the FMA pragma change nothing visible on gcc. In these 30-cycle cells, cmem / physflux off + blend_dtc changes nothing
  visible against FIN; fclip, realis and the Picard counts are also identical.
- **FLAG D6f > P1:** the same entries as FIN in idortfin-more part 2 (T1sp_star E/Fr +3 %, T2cs umbra +9 %, T2sp_shd +2-5 %, T3cs_hot5V aff +1-8 %,
  T3sp_hot5 +2-5 %). All are identical in X1, so they come from the X1-era keys, not from D6f.

## ev_all_d6f.txt (whole; it includes the FIN / X1 / P1 / D6 rows of the same directory)
```
T1sp_star_D6         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 820854 (0.0813793) realis 1.5413e+06/820854 floor 0/0 Pic 2.0678 | E_L1 all 0.1480/+0.0034/0.910 | Fr_L1 all 0.1447 | Drr_L1 0.0189 | f_far 0.9859 (ex 0.9917) L_far 0.9957
T1sp_star_D6f        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 820854 (0.0813793) realis 1.5413e+06/820854 floor 0/0 Pic 2.0678 | E_L1 all 0.1480/+0.0034/0.910 | Fr_L1 all 0.1447 | Drr_L1 0.0189 | f_far 0.9859 (ex 0.9917) L_far 0.9957
T1sp_star_FIN        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 820854 (0.0813793) realis 1.5413e+06/820854 floor 0/0 Pic 2.0678 | E_L1 all 0.1480/+0.0034/0.910 | Fr_L1 all 0.1447 | Drr_L1 0.0189 | f_far 0.9859 (ex 0.9917) L_far 0.9957
T1sp_star_P1         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 974978 (0.0846831) realis 1.81995e+06/974978 floor 0/0 Pic 2.24657 | E_L1 all 0.1436/+0.0035/1.014 | Fr_L1 all 0.1409 | Drr_L1 0.0192 | f_far 0.9856 (ex 0.9917) L_far 0.9957
T1sp_star_X1         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 820854 (0.0813793) realis 1.5413e+06/820854 floor 0/0 Pic 2.0678 | E_L1 all 0.1480/+0.0034/0.910 | Fr_L1 all 0.1447 | Drr_L1 0.0189 | f_far 0.9859 (ex 0.9917) L_far 0.9957
T2cs_shdC_D6         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.31001e+06 (0.148668) realis 1.6439e+07/8.31001e+06 floor 0/0 Pic 2.45763 | all E_L1 0.1399 | aff N 508 E_L1 0.1958 bias -0.0585 Fr_L1 0.1963 | aff_far E_L1 0.2289 | umbra N 28 E/J0 0.0273 (ex 0.0016)
T2cs_shdC_D6f        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.31001e+06 (0.148668) realis 1.6439e+07/8.31001e+06 floor 0/0 Pic 2.45763 | all E_L1 0.1399 | aff N 508 E_L1 0.1958 bias -0.0585 Fr_L1 0.1963 | aff_far E_L1 0.2289 | umbra N 28 E/J0 0.0273 (ex 0.0016)
T2cs_shdC_FIN        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.31001e+06 (0.148668) realis 1.6439e+07/8.31001e+06 floor 0/0 Pic 2.45763 | all E_L1 0.1399 | aff N 508 E_L1 0.1958 bias -0.0585 Fr_L1 0.1963 | aff_far E_L1 0.2289 | umbra N 28 E/J0 0.0273 (ex 0.0016)
T2cs_shdC_P1         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 9.6571e+06 (0.157397) realis 1.83039e+07/9.6571e+06 floor 0/0 Pic 2.45206 | all E_L1 0.1435 | aff N 508 E_L1 0.2047 bias -0.0704 Fr_L1 0.2046 | aff_far E_L1 0.2374 | umbra N 28 E/J0 0.0250 (ex 0.0016)
T2cs_shdC_X1         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.31001e+06 (0.148668) realis 1.6439e+07/8.31001e+06 floor 0/0 Pic 2.45763 | all E_L1 0.1399 | aff N 508 E_L1 0.1958 bias -0.0585 Fr_L1 0.1963 | aff_far E_L1 0.2289 | umbra N 28 E/J0 0.0273 (ex 0.0016)
T2sp_shd_D6          cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 827675 (0.818787) realis 1.68144e+06/827683 floor 0/0 Pic 2.35593 | all E_L1 0.1614 | aff N 519 E_L1 0.1537 bias +0.0240 Fr_L1 0.1558 | aff_far E_L1 0.1778 | umbra N 27 E/J0 0.0446 (ex 0.0073)
T2sp_shd_D6f         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 827675 (0.818787) realis 1.68144e+06/827683 floor 0/0 Pic 2.35593 | all E_L1 0.1614 | aff N 519 E_L1 0.1537 bias +0.0240 Fr_L1 0.1558 | aff_far E_L1 0.1778 | umbra N 27 E/J0 0.0446 (ex 0.0073)
T2sp_shd_FIN         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 827675 (0.818787) realis 1.68144e+06/827683 floor 0/0 Pic 2.35593 | all E_L1 0.1614 | aff N 519 E_L1 0.1537 bias +0.0240 Fr_L1 0.1558 | aff_far E_L1 0.1778 | umbra N 27 E/J0 0.0446 (ex 0.0073)
T2sp_shd_P1          cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 984166 (0.851229) realis 1.9024e+06/984174 floor 0/0 Pic 2.43836 | all E_L1 0.1538 | aff N 519 E_L1 0.1496 bias +0.0291 Fr_L1 0.1503 | aff_far E_L1 0.1738 | umbra N 27 E/J0 0.0450 (ex 0.0073)
T2sp_shd_X1          cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 827675 (0.818787) realis 1.68144e+06/827683 floor 0/0 Pic 2.35593 | all E_L1 0.1614 | aff N 519 E_L1 0.1537 bias +0.0240 Fr_L1 0.1558 | aff_far E_L1 0.1778 | umbra N 27 E/J0 0.0446 (ex 0.0073)
T3cs_hot5V_D6        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.13444e+06 (0.0709192) realis 2.41346e+07/8.13444e+06 floor 0/0 Pic 2.66102 | all E_L1 0.1571 | aff N 1413 E_L1 0.1593 bias +0.0188 Fr_L1 0.1560 | aff_far E_L1 0.2327 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1040
T3cs_hot5V_D6f       cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.13444e+06 (0.0709192) realis 2.41346e+07/8.13444e+06 floor 0/0 Pic 2.66102 | all E_L1 0.1571 | aff N 1413 E_L1 0.1593 bias +0.0188 Fr_L1 0.1560 | aff_far E_L1 0.2327 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1040
T3cs_hot5V_FIN       cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.13444e+06 (0.0709192) realis 2.41344e+07/8.13444e+06 floor 0/0 Pic 2.66102 | all E_L1 0.1571 | aff N 1413 E_L1 0.1593 bias +0.0188 Fr_L1 0.1560 | aff_far E_L1 0.2327 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1040
T3cs_hot5V_P1        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 9.45774e+06 (0.0746801) realis 3.43201e+07/9.45774e+06 floor 0/0 Pic 3.46575 | all E_L1 0.1572 | aff N 1413 E_L1 0.1547 bias +0.0174 Fr_L1 0.1509 | aff_far E_L1 0.2301 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1079
T3cs_hot5V_X1        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.13446e+06 (0.0709188) realis 2.41347e+07/8.13446e+06 floor 0/0 Pic 2.66102 | all E_L1 0.1571 | aff N 1413 E_L1 0.1592 bias +0.0187 Fr_L1 0.1558 | aff_far E_L1 0.2325 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1045
T3sp_hot5_D6         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 747573 (0.0811779) realis 1.6201e+06/747573 floor 0/0 Pic 2.13559 | all E_L1 0.1596 | aff N 1127 E_L1 0.1406 bias +0.0666 Fr_L1 0.1358 | aff_far E_L1 0.1804 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1008
T3sp_hot5_D6f        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 747573 (0.0811779) realis 1.6201e+06/747573 floor 0/0 Pic 2.13559 | all E_L1 0.1596 | aff N 1127 E_L1 0.1406 bias +0.0666 Fr_L1 0.1358 | aff_far E_L1 0.1804 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1008
T3sp_hot5_FIN        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 747567 (0.0811779) realis 1.6201e+06/747567 floor 0/0 Pic 2.13559 | all E_L1 0.1596 | aff N 1127 E_L1 0.1406 bias +0.0666 Fr_L1 0.1358 | aff_far E_L1 0.1804 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1007
T3sp_hot5_P1         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 900353 (0.0834127) realis 1.99954e+06/900353 floor 0/0 Pic 2.42466 | all E_L1 0.1526 | aff N 1127 E_L1 0.1357 bias +0.0656 Fr_L1 0.1313 | aff_far E_L1 0.1731 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1023
T3sp_hot5_X1         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 747564 (0.0811782) realis 1.6201e+06/747564 floor 0/0 Pic 2.13559 | all E_L1 0.1596 | aff N 1127 E_L1 0.1406 bias +0.0666 Fr_L1 0.1358 | aff_far E_L1 0.1804 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1007
```
