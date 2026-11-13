# NOTE 2026-11-13 Orion -> viper: idortd6int: WAITING for the push of fork/d6int-1113 >= 62abe7ec1

GitHub `jing-ze-ma/athenak` (Orion's origin) has `refs/heads/d6int-1113` = **bdc90697bb26** (dhjcuda2-1113 hydro flux kernels, ~2 h old).
It does NOT contain 35af2e3bc, and **62abe7ec1 does not exist on GitHub** (`git cat-file -t 62abe7ec1`: not a valid object, after a fresh fetch).
The TASK requires a d6int build that contains 62abe7ec1, so please push fork/d6int-1113 (with 62abe7ec1) to GitHub.

Orion watches the branch and starts (a), (b) and (c) as soon as the tip contains 62abe7ec1. Meanwhile I am preparing: refs unpacked and
md5-checked, the FIN OpenMP build on 421a3d024 (wt_fin build_none has Athena_ENABLE_OPENMP=OFF), the part (c) scripts, and the wave lists.

---
# idortd6int RESULTS (d6int build), posted per part: (b) below; (a) and (c) follow

## Header
- **d6int = fork/d6int-1113 at build time: 62abe7ec1099d0d30780aaba332db2d1fdb6514a** (pushed ~05:15 PST after the waiting NOTE above).
  It contains 62abe7ec1, **but it does NOT contain 35af2e3bc** (`git merge-base --is-ancestor 35af2e3bc 62abe7ec1` is false). That does not matter here:
  all my builds are MPI ON, so the MPI-off fix is not needed.
- Builds (git archive + kokkos from the Orion src_vfix copy; CMake = the wt_fin build_none options: mpicxx gcc/13 openmpi/4.1.8, Athena_ENABLE_MPI=ON,
  CMAKE_BUILD_TYPE=Release, CMAKE_CXX_FLAGS=-ffp-contract=off, no PROBLEM; wt_fin build_none has Athena_ENABLE_OPENMP=OFF, so OMP builds were made):

| binary | md5 | athena -c openmp lines |
|---|---|---|
| d6int build_none (EXE, EXE_P) | 1d15e7add5e31bfd0491c695b092b1c7 | - |
| d6int build_omp (+Athena_ENABLE_OPENMP=ON) | 1da22bebf0a4c5d00482fc960bfb2b8b | 1 |
| wt_fin 421a3d024 build_none (EXE_F, FIN references) | e151dc3ba4a930d2d348fee069aa50a0 | - |
| 421a3d024 build_omp (EXE_OMP_F) | 422b8f143b3ee404f19af7cbd3458e88 | 1 |

- Runs: 5 p.exclusive nodes started 05:53 / 05:58 PST: (a)+(b) orion014; (c) mainrk orion015, rsx orion022, spxx orion023, omp orion025.
  Singletons with OMPI_MCA_hwloc_base_binding_policy=none (checked: 1 process per CPU).
- Smokes: (a) N16 vetf 18/18, edd 18/18; (b) T2cs_shdC + T1sp_star nlim 3, rc 0 0, FATAL 0; (c) list_smoke 5/5 rows "p rc=0 fatal=0 ... f rc=0 fatal=0",
  OMP 14/14 rc=0 fatal=0, `grep -L implicit_blend_dtc=1e-6` on the p timekeys prints nothing.

## (b) The 5 sp/cs cells, D6 on d6int: **PASS in all 5** (D6i <= 1.05 x X1; D6i = FIN = D6f to 4 digits)
| run | rc | FATAL | wall s | Picard mean / max / NC |
|---|---|---|---|---|
| T1sp_star_D6i | 0 | 0 | 420 | 2.067797e+00 / 3.000000e+00 / 0.000000e+00 |
| T2cs_shdC_D6i | 0 | 0 | 335 | 2.457627e+00 / 3.000000e+00 / 0.000000e+00 |
| T2sp_shd_D6i | 0 | 0 | 251 | 2.355932e+00 / 3.000000e+00 / 0.000000e+00 |
| T3cs_hot5V_D6i | 0 | 0 | 359 | 2.661017e+00 / 1.100000e+01 / 0.000000e+00 |
| T3sp_hot5_D6i | 0 | 0 | 246 | 2.135593e+00 / 5.000000e+00 / 0.000000e+00 |

| cell | metric | D6i | FIN | X1 | P1 | D6f | D6i/X1 | D6i/FIN | D6i/P1 | flag |
|---|---|---|---|---|---|---|---|---|---|---|
| T1sp_star | E_L1 all | 0.1480 | 0.1480 | 0.1480 | 0.1436 | 0.1480 | 1.0000 | 1.0000 | 1.0306 | D6i>P1 |
| T1sp_star | Fr_L1 all | 0.1447 | 0.1447 | 0.1447 | 0.1409 | 0.1447 | 1.0000 | 1.0000 | 1.0270 | D6i>P1 |
| T1sp_star | Drr_L1 | 0.0189 | 0.0189 | 0.0189 | 0.0192 | 0.0189 | 1.0000 | 1.0000 | 0.9844 |  |
| T2cs_shdC | all E_L1 | 0.1399 | 0.1399 | 0.1399 | 0.1435 | 0.1399 | 1.0000 | 1.0000 | 0.9749 |  |
| T2cs_shdC | E_L1 | 0.1958 | 0.1958 | 0.1958 | 0.2047 | 0.1958 | 1.0000 | 1.0000 | 0.9565 |  |
| T2cs_shdC | bias | -0.0585 | -0.0585 | -0.0585 | -0.0704 | -0.0585 | 1.0000 | 1.0000 | 0.8310 |  |
| T2cs_shdC | Fr_L1 | 0.1963 | 0.1963 | 0.1963 | 0.2046 | 0.1963 | 1.0000 | 1.0000 | 0.9594 |  |
| T2cs_shdC | aff_far E_L1 | 0.2289 | 0.2289 | 0.2289 | 0.2374 | 0.2289 | 1.0000 | 1.0000 | 0.9642 |  |
| T2cs_shdC | umbra E/J0 | 0.0273 | 0.0273 | 0.0273 | 0.0250 | 0.0273 | 1.0000 | 1.0000 | 1.0920 | D6i>P1 |
| T2sp_shd | all E_L1 | 0.1614 | 0.1614 | 0.1614 | 0.1538 | 0.1614 | 1.0000 | 1.0000 | 1.0494 | D6i>P1 |
| T2sp_shd | E_L1 | 0.1537 | 0.1537 | 0.1537 | 0.1496 | 0.1537 | 1.0000 | 1.0000 | 1.0274 | D6i>P1 |
| T2sp_shd | bias | 0.0240 | 0.0240 | 0.0240 | 0.0291 | 0.0240 | 1.0000 | 1.0000 | 0.8247 |  |
| T2sp_shd | Fr_L1 | 0.1558 | 0.1558 | 0.1558 | 0.1503 | 0.1558 | 1.0000 | 1.0000 | 1.0366 | D6i>P1 |
| T2sp_shd | aff_far E_L1 | 0.1778 | 0.1778 | 0.1778 | 0.1738 | 0.1778 | 1.0000 | 1.0000 | 1.0230 | D6i>P1 |
| T2sp_shd | umbra E/J0 | 0.0446 | 0.0446 | 0.0446 | 0.0450 | 0.0446 | 1.0000 | 1.0000 | 0.9911 |  |
| T3cs_hot5V | all E_L1 | 0.1571 | 0.1571 | 0.1571 | 0.1572 | 0.1571 | 1.0000 | 1.0000 | 0.9994 |  |
| T3cs_hot5V | E_L1 | 0.1593 | 0.1593 | 0.1592 | 0.1547 | 0.1593 | 1.0006 | 1.0000 | 1.0297 | D6i>P1 |
| T3cs_hot5V | bias | 0.0188 | 0.0188 | 0.0187 | 0.0174 | 0.0188 | 1.0053 | 1.0000 | 1.0805 | D6i>P1 |
| T3cs_hot5V | Fr_L1 | 0.1560 | 0.1560 | 0.1558 | 0.1509 | 0.1560 | 1.0013 | 1.0000 | 1.0338 | D6i>P1 |
| T3cs_hot5V | aff_far E_L1 | 0.2327 | 0.2327 | 0.2325 | 0.2301 | 0.2327 | 1.0009 | 1.0000 | 1.0113 | D6i>P1 |
| T3cs_hot5V | inside E_L1 | 0.1040 | 0.1040 | 0.1045 | 0.1079 | 0.1040 | 0.9952 | 1.0000 | 0.9639 |  |
| T3sp_hot5 | all E_L1 | 0.1596 | 0.1596 | 0.1596 | 0.1526 | 0.1596 | 1.0000 | 1.0000 | 1.0459 | D6i>P1 |
| T3sp_hot5 | E_L1 | 0.1406 | 0.1406 | 0.1406 | 0.1357 | 0.1406 | 1.0000 | 1.0000 | 1.0361 | D6i>P1 |
| T3sp_hot5 | bias | 0.0666 | 0.0666 | 0.0666 | 0.0656 | 0.0666 | 1.0000 | 1.0000 | 1.0152 | D6i>P1 |
| T3sp_hot5 | Fr_L1 | 0.1358 | 0.1358 | 0.1358 | 0.1313 | 0.1358 | 1.0000 | 1.0000 | 1.0343 | D6i>P1 |
| T3sp_hot5 | aff_far E_L1 | 0.1804 | 0.1804 | 0.1804 | 0.1731 | 0.1804 | 1.0000 | 1.0000 | 1.0422 | D6i>P1 |
| T3sp_hot5 | inside E_L1 | 0.1008 | 0.1007 | 0.1007 | 0.1023 | 0.1008 | 1.0010 | 1.0010 | 0.9853 |  |

- 0 metrics with D6i > 1.05 x X1 (max 1.0053, T3cs_hot5V aff bias, the same as FIN and D6f). The D6i/FIN and D6i/D6f columns are 1.0000 except
  T3sp_hot5 inside (+0.1 %, the same as D6f). The D6i > P1 flags are the X1-era entries already reported for FIN / D6f.

### ev_all_d6i.txt (whole)
```
T1sp_star_D6         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 820854 (0.0813793) realis 1.5413e+06/820854 floor 0/0 Pic 2.0678 | E_L1 all 0.1480/+0.0034/0.910 | Fr_L1 all 0.1447 | Drr_L1 0.0189 | f_far 0.9859 (ex 0.9917) L_far 0.9957
T1sp_star_D6f        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 820854 (0.0813793) realis 1.5413e+06/820854 floor 0/0 Pic 2.0678 | E_L1 all 0.1480/+0.0034/0.910 | Fr_L1 all 0.1447 | Drr_L1 0.0189 | f_far 0.9859 (ex 0.9917) L_far 0.9957
T1sp_star_D6i        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 820854 (0.0813793) realis 1.5413e+06/820854 floor 0/0 Pic 2.0678 | E_L1 all 0.1480/+0.0034/0.910 | Fr_L1 all 0.1447 | Drr_L1 0.0189 | f_far 0.9859 (ex 0.9917) L_far 0.9957
T1sp_star_FIN        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 820854 (0.0813793) realis 1.5413e+06/820854 floor 0/0 Pic 2.0678 | E_L1 all 0.1480/+0.0034/0.910 | Fr_L1 all 0.1447 | Drr_L1 0.0189 | f_far 0.9859 (ex 0.9917) L_far 0.9957
T1sp_star_P1         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 974978 (0.0846831) realis 1.81995e+06/974978 floor 0/0 Pic 2.24657 | E_L1 all 0.1436/+0.0035/1.014 | Fr_L1 all 0.1409 | Drr_L1 0.0192 | f_far 0.9856 (ex 0.9917) L_far 0.9957
T1sp_star_X1         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 820854 (0.0813793) realis 1.5413e+06/820854 floor 0/0 Pic 2.0678 | E_L1 all 0.1480/+0.0034/0.910 | Fr_L1 all 0.1447 | Drr_L1 0.0189 | f_far 0.9859 (ex 0.9917) L_far 0.9957
T2cs_shdC_D6         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.31001e+06 (0.148668) realis 1.6439e+07/8.31001e+06 floor 0/0 Pic 2.45763 | all E_L1 0.1399 | aff N 508 E_L1 0.1958 bias -0.0585 Fr_L1 0.1963 | aff_far E_L1 0.2289 | umbra N 28 E/J0 0.0273 (ex 0.0016)
T2cs_shdC_D6f        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.31001e+06 (0.148668) realis 1.6439e+07/8.31001e+06 floor 0/0 Pic 2.45763 | all E_L1 0.1399 | aff N 508 E_L1 0.1958 bias -0.0585 Fr_L1 0.1963 | aff_far E_L1 0.2289 | umbra N 28 E/J0 0.0273 (ex 0.0016)
T2cs_shdC_D6i        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.31001e+06 (0.148668) realis 1.6439e+07/8.31001e+06 floor 0/0 Pic 2.45763 | all E_L1 0.1399 | aff N 508 E_L1 0.1958 bias -0.0585 Fr_L1 0.1963 | aff_far E_L1 0.2289 | umbra N 28 E/J0 0.0273 (ex 0.0016)
T2cs_shdC_FIN        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.31001e+06 (0.148668) realis 1.6439e+07/8.31001e+06 floor 0/0 Pic 2.45763 | all E_L1 0.1399 | aff N 508 E_L1 0.1958 bias -0.0585 Fr_L1 0.1963 | aff_far E_L1 0.2289 | umbra N 28 E/J0 0.0273 (ex 0.0016)
T2cs_shdC_P1         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 9.6571e+06 (0.157397) realis 1.83039e+07/9.6571e+06 floor 0/0 Pic 2.45206 | all E_L1 0.1435 | aff N 508 E_L1 0.2047 bias -0.0704 Fr_L1 0.2046 | aff_far E_L1 0.2374 | umbra N 28 E/J0 0.0250 (ex 0.0016)
T2cs_shdC_X1         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.31001e+06 (0.148668) realis 1.6439e+07/8.31001e+06 floor 0/0 Pic 2.45763 | all E_L1 0.1399 | aff N 508 E_L1 0.1958 bias -0.0585 Fr_L1 0.1963 | aff_far E_L1 0.2289 | umbra N 28 E/J0 0.0273 (ex 0.0016)
T2sp_shd_D6          cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 827675 (0.818787) realis 1.68144e+06/827683 floor 0/0 Pic 2.35593 | all E_L1 0.1614 | aff N 519 E_L1 0.1537 bias +0.0240 Fr_L1 0.1558 | aff_far E_L1 0.1778 | umbra N 27 E/J0 0.0446 (ex 0.0073)
T2sp_shd_D6f         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 827675 (0.818787) realis 1.68144e+06/827683 floor 0/0 Pic 2.35593 | all E_L1 0.1614 | aff N 519 E_L1 0.1537 bias +0.0240 Fr_L1 0.1558 | aff_far E_L1 0.1778 | umbra N 27 E/J0 0.0446 (ex 0.0073)
T2sp_shd_D6i         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 827675 (0.818787) realis 1.68144e+06/827683 floor 0/0 Pic 2.35593 | all E_L1 0.1614 | aff N 519 E_L1 0.1537 bias +0.0240 Fr_L1 0.1558 | aff_far E_L1 0.1778 | umbra N 27 E/J0 0.0446 (ex 0.0073)
T2sp_shd_FIN         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 827675 (0.818787) realis 1.68144e+06/827683 floor 0/0 Pic 2.35593 | all E_L1 0.1614 | aff N 519 E_L1 0.1537 bias +0.0240 Fr_L1 0.1558 | aff_far E_L1 0.1778 | umbra N 27 E/J0 0.0446 (ex 0.0073)
T2sp_shd_P1          cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 984166 (0.851229) realis 1.9024e+06/984174 floor 0/0 Pic 2.43836 | all E_L1 0.1538 | aff N 519 E_L1 0.1496 bias +0.0291 Fr_L1 0.1503 | aff_far E_L1 0.1738 | umbra N 27 E/J0 0.0450 (ex 0.0073)
T2sp_shd_X1          cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 827675 (0.818787) realis 1.68144e+06/827683 floor 0/0 Pic 2.35593 | all E_L1 0.1614 | aff N 519 E_L1 0.1537 bias +0.0240 Fr_L1 0.1558 | aff_far E_L1 0.1778 | umbra N 27 E/J0 0.0446 (ex 0.0073)
T3cs_hot5V_D6        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.13444e+06 (0.0709192) realis 2.41346e+07/8.13444e+06 floor 0/0 Pic 2.66102 | all E_L1 0.1571 | aff N 1413 E_L1 0.1593 bias +0.0188 Fr_L1 0.1560 | aff_far E_L1 0.2327 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1040
T3cs_hot5V_D6f       cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.13444e+06 (0.0709192) realis 2.41346e+07/8.13444e+06 floor 0/0 Pic 2.66102 | all E_L1 0.1571 | aff N 1413 E_L1 0.1593 bias +0.0188 Fr_L1 0.1560 | aff_far E_L1 0.2327 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1040
T3cs_hot5V_D6i       cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.13444e+06 (0.0709192) realis 2.41346e+07/8.13444e+06 floor 0/0 Pic 2.66102 | all E_L1 0.1571 | aff N 1413 E_L1 0.1593 bias +0.0188 Fr_L1 0.1560 | aff_far E_L1 0.2327 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1040
T3cs_hot5V_FIN       cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.13444e+06 (0.0709192) realis 2.41344e+07/8.13444e+06 floor 0/0 Pic 2.66102 | all E_L1 0.1571 | aff N 1413 E_L1 0.1593 bias +0.0188 Fr_L1 0.1560 | aff_far E_L1 0.2327 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1040
T3cs_hot5V_P1        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 9.45774e+06 (0.0746801) realis 3.43201e+07/9.45774e+06 floor 0/0 Pic 3.46575 | all E_L1 0.1572 | aff N 1413 E_L1 0.1547 bias +0.0174 Fr_L1 0.1509 | aff_far E_L1 0.2301 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1079
T3cs_hot5V_X1        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.13446e+06 (0.0709188) realis 2.41347e+07/8.13446e+06 floor 0/0 Pic 2.66102 | all E_L1 0.1571 | aff N 1413 E_L1 0.1592 bias +0.0187 Fr_L1 0.1558 | aff_far E_L1 0.2325 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1045
T3sp_hot5_D6         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 747573 (0.0811779) realis 1.6201e+06/747573 floor 0/0 Pic 2.13559 | all E_L1 0.1596 | aff N 1127 E_L1 0.1406 bias +0.0666 Fr_L1 0.1358 | aff_far E_L1 0.1804 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1008
T3sp_hot5_D6f        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 747573 (0.0811779) realis 1.6201e+06/747573 floor 0/0 Pic 2.13559 | all E_L1 0.1596 | aff N 1127 E_L1 0.1406 bias +0.0666 Fr_L1 0.1358 | aff_far E_L1 0.1804 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1008
T3sp_hot5_D6i        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 747573 (0.0811779) realis 1.6201e+06/747573 floor 0/0 Pic 2.13559 | all E_L1 0.1596 | aff N 1127 E_L1 0.1406 bias +0.0666 Fr_L1 0.1358 | aff_far E_L1 0.1804 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1008
T3sp_hot5_FIN        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 747567 (0.0811779) realis 1.6201e+06/747567 floor 0/0 Pic 2.13559 | all E_L1 0.1596 | aff N 1127 E_L1 0.1406 bias +0.0666 Fr_L1 0.1358 | aff_far E_L1 0.1804 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1007
T3sp_hot5_P1         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 900353 (0.0834127) realis 1.99954e+06/900353 floor 0/0 Pic 2.42466 | all E_L1 0.1526 | aff N 1127 E_L1 0.1357 bias +0.0656 Fr_L1 0.1313 | aff_far E_L1 0.1731 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1023
T3sp_hot5_X1         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 747564 (0.0811782) realis 1.6201e+06/747564 floor 0/0 Pic 2.13559 | all E_L1 0.1596 | aff N 1127 E_L1 0.1406 bias +0.0666 Fr_L1 0.1358 | aff_far E_L1 0.1804 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1007
```
