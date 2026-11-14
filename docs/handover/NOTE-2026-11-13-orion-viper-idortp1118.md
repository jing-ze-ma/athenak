# NOTE 2026-11-13 Orion -> viper: idortp1118 (prod-1118 + production keys), posted per part: (b) first; (a), (c) follow

## Header
- **prod-1118 = fork/prod-1118 = 1d1972805568f2e551b01914f6f01ba80d0ed500** (git archive; kokkos from the Orion src_vfix copy). CMake = the d6int
  build_none options (mpicxx gcc/13 openmpi/4.1.8, MPI ON, Release, -ffp-contract=off); build_omp = + Athena_ENABLE_OPENMP=ON.
  md5: build_none **55d28640edb21a9971e574fae635f11c**, build_omp **579d58f90d216e2f0c8f9e7cc8ec977a** (athena -c: OpenMP line present).
  References reused, not rerun: d6int 62abe7ec1 (1d15e7ad...) and FIN 421a3d024 (e151dc3b...).
- **UPDATE applied (TASK ...-UPDATE): arm q in (a) and (c) = EXTRA_HG.** Part (b) runs both P18 (EXTRA) and P18H (EXTRA_HG). No q0 runs exist:
  the update arrived before any q run started. Only the (c) smoke ran with plain EXTRA (5/5 rc=0).
- **Key injection (needed):** zalesak_cap_rounds, implicit_lzx and the three HG keys are in NONE of the 153 decks_d6 (time2_stage_local_be is in all of them).
  So one_q.sh as posted (EXTRA on the command line) would FATAL ("Parameter ... on command line not found"). I used an Orion copy
  (one_q_orion.sh + splitkeys.py): keys absent from the deck are inserted right after <rad_m1> in the run's copy of in.athinput, and keys present
  (time2_stage_local_be=false) stay on the command line, where they override the deck's true. The same splitkeys was used for the (b) decks. For (a), the
  tokens rm_zalesak_cap_rounds=4 rm_implicit_lzx=true rm_time2_stage_local_be=false (+ rm_implicit_vimp_lag_local=true rm_time2_adm_poslim=true
  rm_time2_dtr_poslim=true) were appended to every case line (drive.py writes them after the TIME/D6 keys; checked: the deck shows stage_local_be = false).
- Runs: 5 p.exclusive nodes, 21:30 CET 11-13 = 12:30 PST (queue wait ~2.5 h on priority): ab orion051, c main 052, rsx 053, spxx 054, omp 055.

## (b) 5 sp/cs cells, P18 (EXTRA) and P18H (EXTRA_HG): **PASS in all 5** (every metric <= 1.05 x X1; ev metrics = D6i to 4 digits)
| run | rc | FATAL | wall s | Picard mean / max / NC | cmdline key |
|---|---|---|---|---|---|
| T1sp_star_P18 | 0 | 0 | 428 | 2.067797e+00 / 3.000000e+00 / 0.000000e+00 | rad_m1/time2_stage_local_be=false |
| T2cs_shdC_P18 | 0 | 0 | 360 | 2.754098e+00 / 7.000000e+00 / 0.000000e+00 | rad_m1/time2_stage_local_be=false |
| T2sp_shd_P18 | 0 | 0 | 276 | 2.770492e+00 / 8.000000e+00 / 0.000000e+00 | rad_m1/time2_stage_local_be=false |
| T3cs_hot5V_P18 | 0 | 0 | 364 | 2.661017e+00 / 1.100000e+01 / 0.000000e+00 | rad_m1/time2_stage_local_be=false |
| T3sp_hot5_P18 | 0 | 0 | 251 | 2.135593e+00 / 5.000000e+00 / 0.000000e+00 | rad_m1/time2_stage_local_be=false |
| T1sp_star_P18H | 0 | 0 | 429 | 2.067797e+00 / 3.000000e+00 / 0.000000e+00 | rad_m1/time2_stage_local_be=false |
| T2cs_shdC_P18H | 0 | 0 | 365 | 2.754098e+00 / 7.000000e+00 / 0.000000e+00 | rad_m1/time2_stage_local_be=false |
| T2sp_shd_P18H | 0 | 0 | 280 | 2.770492e+00 / 8.000000e+00 / 0.000000e+00 | rad_m1/time2_stage_local_be=false |
| T3cs_hot5V_P18H | 0 | 0 | 371 | 2.661017e+00 / 1.100000e+01 / 0.000000e+00 | rad_m1/time2_stage_local_be=false |
| T3sp_hot5_P18H | 0 | 0 | 253 | 2.135593e+00 / 5.000000e+00 / 0.000000e+00 | rad_m1/time2_stage_local_be=false |

**q vs d6int, bitwise (all 4 .bin outputs per cell, bin_convert data compare):**
- T1sp_star, T3cs_hot5V, T3sp_hot5: **P18 = P18H = D6i bitwise.**
- T2cs_shdC: P18 = P18H; **both != D6i** in ta.m1.00001 (m1 vars): max|d|/max|x| 4e-9 to 9e-9, L1 rel 2.5e-12 to 6e-12, 68-170 of 442368 cells.
- T2sp_shd: P18 = P18H; **both != D6i** in ta.m1.00001: max|d|/max|x| 2.2e-7 to 4.7e-7, L1 rel 1.8e-8 to 4.8e-8, ~20 % of cells (10750-13135 / 55296).
- These two cells also run differently: face_mmat limited x1 faces 32 -> 475 (shdC) and 4 -> 399 (shd); Picard mean 2.46 -> 2.75 (+12 %) and
  2.36 -> 2.77 (+17 %), below the 1.25x FLAG_PICARD line. The cause is not isolated. Candidates: time2_stage_local_be=false (D6 decks have true)
  and the Zalesak cap. hot5V's 2001 limited faces are identical in D6i and P18.
- The hydro .hst files are md5-identical in all 15 runs (D6i, P18, P18H): these cells keep the gas fixed. L_z / energy columns: the ev.py rows carry
  no L_z column; the energy / flux metrics are below.

**Positivity / lzx counters (P18H logs, all 5 cells):** implicit_lzx "inert (needs the cubed sphere or spherical polar and coupling with
gas_feedback)" (these cells have no gas feedback); implicit_vimp_lag_local lag-path triggers 0; time2_dtr_poslim BE solves with eint <= 0 = 0;
time2_adm_poslim admitted stages 0. **The HG keys never fire on these cells**, hence P18 = P18H.

| cell | metric | P18 | P18H | D6i | FIN | X1 | P1 | P18/X1 | P18H/X1 | P18/D6i | P18H/D6i | flag |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T1sp_star | E_L1 all | 0.1480 | 0.1480 | 0.1480 | 0.1480 | 0.1480 | 0.1436 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |  |
| T1sp_star | Fr_L1 all | 0.1447 | 0.1447 | 0.1447 | 0.1447 | 0.1447 | 0.1409 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |  |
| T1sp_star | Drr_L1 | 0.0189 | 0.0189 | 0.0189 | 0.0189 | 0.0189 | 0.0192 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |  |
| T2cs_shdC | all E_L1 | 0.1399 | 0.1399 | 0.1399 | 0.1399 | 0.1399 | 0.1435 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |  |
| T2cs_shdC | E_L1 | 0.1958 | 0.1958 | 0.1958 | 0.1958 | 0.1958 | 0.2047 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |  |
| T2cs_shdC | bias | -0.0585 | -0.0585 | -0.0585 | -0.0585 | -0.0585 | -0.0704 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |  |
| T2cs_shdC | Fr_L1 | 0.1963 | 0.1963 | 0.1963 | 0.1963 | 0.1963 | 0.2046 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |  |
| T2cs_shdC | aff_far E_L1 | 0.2289 | 0.2289 | 0.2289 | 0.2289 | 0.2289 | 0.2374 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |  |
| T2cs_shdC | umbra E/J0 | 0.0273 | 0.0273 | 0.0273 | 0.0273 | 0.0273 | 0.0250 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |  |
| T2sp_shd | all E_L1 | 0.1614 | 0.1614 | 0.1614 | 0.1614 | 0.1614 | 0.1538 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |  |
| T2sp_shd | E_L1 | 0.1537 | 0.1537 | 0.1537 | 0.1537 | 0.1537 | 0.1496 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |  |
| T2sp_shd | bias | 0.0240 | 0.0240 | 0.0240 | 0.0240 | 0.0240 | 0.0291 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |  |
| T2sp_shd | Fr_L1 | 0.1558 | 0.1558 | 0.1558 | 0.1558 | 0.1558 | 0.1503 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |  |
| T2sp_shd | aff_far E_L1 | 0.1778 | 0.1778 | 0.1778 | 0.1778 | 0.1778 | 0.1738 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |  |
| T2sp_shd | umbra E/J0 | 0.0446 | 0.0446 | 0.0446 | 0.0446 | 0.0446 | 0.0450 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |  |
| T3cs_hot5V | all E_L1 | 0.1571 | 0.1571 | 0.1571 | 0.1571 | 0.1571 | 0.1572 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |  |
| T3cs_hot5V | E_L1 | 0.1593 | 0.1593 | 0.1593 | 0.1593 | 0.1592 | 0.1547 | 1.0006 | 1.0006 | 1.0000 | 1.0000 |  |
| T3cs_hot5V | bias | 0.0188 | 0.0188 | 0.0188 | 0.0188 | 0.0187 | 0.0174 | 1.0053 | 1.0053 | 1.0000 | 1.0000 |  |
| T3cs_hot5V | Fr_L1 | 0.1560 | 0.1560 | 0.1560 | 0.1560 | 0.1558 | 0.1509 | 1.0013 | 1.0013 | 1.0000 | 1.0000 |  |
| T3cs_hot5V | aff_far E_L1 | 0.2327 | 0.2327 | 0.2327 | 0.2327 | 0.2325 | 0.2301 | 1.0009 | 1.0009 | 1.0000 | 1.0000 |  |
| T3cs_hot5V | inside E_L1 | 0.1040 | 0.1040 | 0.1040 | 0.1040 | 0.1045 | 0.1079 | 0.9952 | 0.9952 | 1.0000 | 1.0000 |  |
| T3sp_hot5 | all E_L1 | 0.1596 | 0.1596 | 0.1596 | 0.1596 | 0.1596 | 0.1526 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |  |
| T3sp_hot5 | E_L1 | 0.1406 | 0.1406 | 0.1406 | 0.1406 | 0.1406 | 0.1357 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |  |
| T3sp_hot5 | bias | 0.0666 | 0.0666 | 0.0666 | 0.0666 | 0.0666 | 0.0656 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |  |
| T3sp_hot5 | Fr_L1 | 0.1358 | 0.1358 | 0.1358 | 0.1358 | 0.1358 | 0.1313 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |  |
| T3sp_hot5 | aff_far E_L1 | 0.1804 | 0.1804 | 0.1804 | 0.1804 | 0.1804 | 0.1731 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |  |
| T3sp_hot5 | inside E_L1 | 0.1008 | 0.1008 | 0.1008 | 0.1007 | 0.1007 | 0.1023 | 1.0010 | 1.0010 | 1.0000 | 1.0000 |  |

### ev_all_p18.txt (whole)
```
T1sp_star_D6         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 820854 (0.0813793) realis 1.5413e+06/820854 floor 0/0 Pic 2.0678 | E_L1 all 0.1480/+0.0034/0.910 | Fr_L1 all 0.1447 | Drr_L1 0.0189 | f_far 0.9859 (ex 0.9917) L_far 0.9957
T1sp_star_D6f        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 820854 (0.0813793) realis 1.5413e+06/820854 floor 0/0 Pic 2.0678 | E_L1 all 0.1480/+0.0034/0.910 | Fr_L1 all 0.1447 | Drr_L1 0.0189 | f_far 0.9859 (ex 0.9917) L_far 0.9957
T1sp_star_D6i        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 820854 (0.0813793) realis 1.5413e+06/820854 floor 0/0 Pic 2.0678 | E_L1 all 0.1480/+0.0034/0.910 | Fr_L1 all 0.1447 | Drr_L1 0.0189 | f_far 0.9859 (ex 0.9917) L_far 0.9957
T1sp_star_FIN        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 820854 (0.0813793) realis 1.5413e+06/820854 floor 0/0 Pic 2.0678 | E_L1 all 0.1480/+0.0034/0.910 | Fr_L1 all 0.1447 | Drr_L1 0.0189 | f_far 0.9859 (ex 0.9917) L_far 0.9957
T1sp_star_P1         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 974978 (0.0846831) realis 1.81995e+06/974978 floor 0/0 Pic 2.24657 | E_L1 all 0.1436/+0.0035/1.014 | Fr_L1 all 0.1409 | Drr_L1 0.0192 | f_far 0.9856 (ex 0.9917) L_far 0.9957
T1sp_star_P18        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 820854 (0.0813793) realis 1.5413e+06/820854 floor 0/0 Pic 2.0678 | E_L1 all 0.1480/+0.0034/0.910 | Fr_L1 all 0.1447 | Drr_L1 0.0189 | f_far 0.9859 (ex 0.9917) L_far 0.9957
T1sp_star_P18H       cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 820854 (0.0813793) realis 1.5413e+06/820854 floor 0/0 Pic 2.0678 | E_L1 all 0.1480/+0.0034/0.910 | Fr_L1 all 0.1447 | Drr_L1 0.0189 | f_far 0.9859 (ex 0.9917) L_far 0.9957
T1sp_star_X1         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 820854 (0.0813793) realis 1.5413e+06/820854 floor 0/0 Pic 2.0678 | E_L1 all 0.1480/+0.0034/0.910 | Fr_L1 all 0.1447 | Drr_L1 0.0189 | f_far 0.9859 (ex 0.9917) L_far 0.9957
T2cs_shdC_D6         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.31001e+06 (0.148668) realis 1.6439e+07/8.31001e+06 floor 0/0 Pic 2.45763 | all E_L1 0.1399 | aff N 508 E_L1 0.1958 bias -0.0585 Fr_L1 0.1963 | aff_far E_L1 0.2289 | umbra N 28 E/J0 0.0273 (ex 0.0016)
T2cs_shdC_D6f        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.31001e+06 (0.148668) realis 1.6439e+07/8.31001e+06 floor 0/0 Pic 2.45763 | all E_L1 0.1399 | aff N 508 E_L1 0.1958 bias -0.0585 Fr_L1 0.1963 | aff_far E_L1 0.2289 | umbra N 28 E/J0 0.0273 (ex 0.0016)
T2cs_shdC_D6i        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.31001e+06 (0.148668) realis 1.6439e+07/8.31001e+06 floor 0/0 Pic 2.45763 | all E_L1 0.1399 | aff N 508 E_L1 0.1958 bias -0.0585 Fr_L1 0.1963 | aff_far E_L1 0.2289 | umbra N 28 E/J0 0.0273 (ex 0.0016)
T2cs_shdC_FIN        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.31001e+06 (0.148668) realis 1.6439e+07/8.31001e+06 floor 0/0 Pic 2.45763 | all E_L1 0.1399 | aff N 508 E_L1 0.1958 bias -0.0585 Fr_L1 0.1963 | aff_far E_L1 0.2289 | umbra N 28 E/J0 0.0273 (ex 0.0016)
T2cs_shdC_P1         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 9.6571e+06 (0.157397) realis 1.83039e+07/9.6571e+06 floor 0/0 Pic 2.45206 | all E_L1 0.1435 | aff N 508 E_L1 0.2047 bias -0.0704 Fr_L1 0.2046 | aff_far E_L1 0.2374 | umbra N 28 E/J0 0.0250 (ex 0.0016)
T2cs_shdC_P18        cyc 30 NA 1 NP 2 NC 0 F 0 st/BE 2.700000e+01/3.000000e+00 stfb 1 fclip 9.11576e+06 (0.133114) realis 2.32197e+07/9.1158e+06 floor 33/0 Pic 2.7541 | all E_L1 0.1399 | aff N 508 E_L1 0.1958 bias -0.0585 Fr_L1 0.1963 | aff_far E_L1 0.2289 | umbra N 28 E/J0 0.0273 (ex 0.0016)
T2cs_shdC_P18H       cyc 30 NA 1 NP 2 NC 0 F 0 st/BE 2.700000e+01/3.000000e+00 stfb 1 fclip 9.11576e+06 (0.133114) realis 2.32197e+07/9.1158e+06 floor 33/0 Pic 2.7541 | all E_L1 0.1399 | aff N 508 E_L1 0.1958 bias -0.0585 Fr_L1 0.1963 | aff_far E_L1 0.2289 | umbra N 28 E/J0 0.0273 (ex 0.0016)
T2cs_shdC_X1         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.31001e+06 (0.148668) realis 1.6439e+07/8.31001e+06 floor 0/0 Pic 2.45763 | all E_L1 0.1399 | aff N 508 E_L1 0.1958 bias -0.0585 Fr_L1 0.1963 | aff_far E_L1 0.2289 | umbra N 28 E/J0 0.0273 (ex 0.0016)
T2sp_shd_D6          cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 827675 (0.818787) realis 1.68144e+06/827683 floor 0/0 Pic 2.35593 | all E_L1 0.1614 | aff N 519 E_L1 0.1537 bias +0.0240 Fr_L1 0.1558 | aff_far E_L1 0.1778 | umbra N 27 E/J0 0.0446 (ex 0.0073)
T2sp_shd_D6f         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 827675 (0.818787) realis 1.68144e+06/827683 floor 0/0 Pic 2.35593 | all E_L1 0.1614 | aff N 519 E_L1 0.1537 bias +0.0240 Fr_L1 0.1558 | aff_far E_L1 0.1778 | umbra N 27 E/J0 0.0446 (ex 0.0073)
T2sp_shd_D6i         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 827675 (0.818787) realis 1.68144e+06/827683 floor 0/0 Pic 2.35593 | all E_L1 0.1614 | aff N 519 E_L1 0.1537 bias +0.0240 Fr_L1 0.1558 | aff_far E_L1 0.1778 | umbra N 27 E/J0 0.0446 (ex 0.0073)
T2sp_shd_FIN         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 827675 (0.818787) realis 1.68144e+06/827683 floor 0/0 Pic 2.35593 | all E_L1 0.1614 | aff N 519 E_L1 0.1537 bias +0.0240 Fr_L1 0.1558 | aff_far E_L1 0.1778 | umbra N 27 E/J0 0.0446 (ex 0.0073)
T2sp_shd_P1          cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 984166 (0.851229) realis 1.9024e+06/984174 floor 0/0 Pic 2.43836 | all E_L1 0.1538 | aff N 519 E_L1 0.1496 bias +0.0291 Fr_L1 0.1503 | aff_far E_L1 0.1738 | umbra N 27 E/J0 0.0450 (ex 0.0073)
T2sp_shd_P18         cyc 30 NA 2 NP 2 NC 0 F 0 st/BE 2.500000e+01/5.000000e+00 stfb 2 fclip 903512 (0.764931) realis 2.34681e+06/903544 floor 24/0 Pic 2.77049 | all E_L1 0.1614 | aff N 519 E_L1 0.1537 bias +0.0240 Fr_L1 0.1558 | aff_far E_L1 0.1778 | umbra N 27 E/J0 0.0446 (ex 0.0073)
T2sp_shd_P18H        cyc 30 NA 2 NP 2 NC 0 F 0 st/BE 2.500000e+01/5.000000e+00 stfb 2 fclip 903512 (0.764931) realis 2.34681e+06/903544 floor 24/0 Pic 2.77049 | all E_L1 0.1614 | aff N 519 E_L1 0.1537 bias +0.0240 Fr_L1 0.1558 | aff_far E_L1 0.1778 | umbra N 27 E/J0 0.0446 (ex 0.0073)
T2sp_shd_X1          cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 827675 (0.818787) realis 1.68144e+06/827683 floor 0/0 Pic 2.35593 | all E_L1 0.1614 | aff N 519 E_L1 0.1537 bias +0.0240 Fr_L1 0.1558 | aff_far E_L1 0.1778 | umbra N 27 E/J0 0.0446 (ex 0.0073)
T3cs_hot5V_D6        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.13444e+06 (0.0709192) realis 2.41346e+07/8.13444e+06 floor 0/0 Pic 2.66102 | all E_L1 0.1571 | aff N 1413 E_L1 0.1593 bias +0.0188 Fr_L1 0.1560 | aff_far E_L1 0.2327 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1040
T3cs_hot5V_D6f       cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.13444e+06 (0.0709192) realis 2.41346e+07/8.13444e+06 floor 0/0 Pic 2.66102 | all E_L1 0.1571 | aff N 1413 E_L1 0.1593 bias +0.0188 Fr_L1 0.1560 | aff_far E_L1 0.2327 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1040
T3cs_hot5V_D6i       cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.13444e+06 (0.0709192) realis 2.41346e+07/8.13444e+06 floor 0/0 Pic 2.66102 | all E_L1 0.1571 | aff N 1413 E_L1 0.1593 bias +0.0188 Fr_L1 0.1560 | aff_far E_L1 0.2327 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1040
T3cs_hot5V_FIN       cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.13444e+06 (0.0709192) realis 2.41344e+07/8.13444e+06 floor 0/0 Pic 2.66102 | all E_L1 0.1571 | aff N 1413 E_L1 0.1593 bias +0.0188 Fr_L1 0.1560 | aff_far E_L1 0.2327 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1040
T3cs_hot5V_P1        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 9.45774e+06 (0.0746801) realis 3.43201e+07/9.45774e+06 floor 0/0 Pic 3.46575 | all E_L1 0.1572 | aff N 1413 E_L1 0.1547 bias +0.0174 Fr_L1 0.1509 | aff_far E_L1 0.2301 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1079
T3cs_hot5V_P18       cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.13444e+06 (0.0709192) realis 2.41346e+07/8.13444e+06 floor 0/0 Pic 2.66102 | all E_L1 0.1571 | aff N 1413 E_L1 0.1593 bias +0.0188 Fr_L1 0.1560 | aff_far E_L1 0.2327 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1040
T3cs_hot5V_P18H      cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.13444e+06 (0.0709192) realis 2.41346e+07/8.13444e+06 floor 0/0 Pic 2.66102 | all E_L1 0.1571 | aff N 1413 E_L1 0.1593 bias +0.0188 Fr_L1 0.1560 | aff_far E_L1 0.2327 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1040
T3cs_hot5V_X1        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 8.13446e+06 (0.0709188) realis 2.41347e+07/8.13446e+06 floor 0/0 Pic 2.66102 | all E_L1 0.1571 | aff N 1413 E_L1 0.1592 bias +0.0187 Fr_L1 0.1558 | aff_far E_L1 0.2325 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1045
T3sp_hot5_D6         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 747573 (0.0811779) realis 1.6201e+06/747573 floor 0/0 Pic 2.13559 | all E_L1 0.1596 | aff N 1127 E_L1 0.1406 bias +0.0666 Fr_L1 0.1358 | aff_far E_L1 0.1804 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1008
T3sp_hot5_D6f        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 747573 (0.0811779) realis 1.6201e+06/747573 floor 0/0 Pic 2.13559 | all E_L1 0.1596 | aff N 1127 E_L1 0.1406 bias +0.0666 Fr_L1 0.1358 | aff_far E_L1 0.1804 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1008
T3sp_hot5_D6i        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 747573 (0.0811779) realis 1.6201e+06/747573 floor 0/0 Pic 2.13559 | all E_L1 0.1596 | aff N 1127 E_L1 0.1406 bias +0.0666 Fr_L1 0.1358 | aff_far E_L1 0.1804 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1008
T3sp_hot5_FIN        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 747567 (0.0811779) realis 1.6201e+06/747567 floor 0/0 Pic 2.13559 | all E_L1 0.1596 | aff N 1127 E_L1 0.1406 bias +0.0666 Fr_L1 0.1358 | aff_far E_L1 0.1804 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1007
T3sp_hot5_P1         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 900353 (0.0834127) realis 1.99954e+06/900353 floor 0/0 Pic 2.42466 | all E_L1 0.1526 | aff N 1127 E_L1 0.1357 bias +0.0656 Fr_L1 0.1313 | aff_far E_L1 0.1731 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1023
T3sp_hot5_P18        cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 747573 (0.0811779) realis 1.6201e+06/747573 floor 0/0 Pic 2.13559 | all E_L1 0.1596 | aff N 1127 E_L1 0.1406 bias +0.0666 Fr_L1 0.1358 | aff_far E_L1 0.1804 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1008
T3sp_hot5_P18H       cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 747573 (0.0811779) realis 1.6201e+06/747573 floor 0/0 Pic 2.13559 | all E_L1 0.1596 | aff N 1127 E_L1 0.1406 bias +0.0666 Fr_L1 0.1358 | aff_far E_L1 0.1804 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1008
T3sp_hot5_X1         cyc 30 NA 0 NP 0 NC 0 F 0 st/BE 2.900000e+01/1.000000e+00 stfb 0 fclip 747564 (0.0811782) realis 1.6201e+06/747564 floor 0/0 Pic 2.13559 | all E_L1 0.1596 | aff N 1127 E_L1 0.1406 bias +0.0666 Fr_L1 0.1358 | aff_far E_L1 0.1804 | umbra N 0 E/J0 nan (ex nan) | inside E_L1 0.1007
```

## (a) Wave grid, q = prod-1118 + EXTRA_HG (180 cases): **PASS. q = d6int BITWISE on all 180 cases** (cmp_wbs dev 0 in every closure / P / tau)
Lists fin_P18H_* = fin_D6i_* + the six EXTRA_HG tokens (see the header); run dirs runs_finP18H_*; smoke N16 vetf 18/18, edd 18/18; full on orion051.
**Counts: P18H_edd 90/90, P18H_vetf 90/90, 0 FAILED, 0 FATAL.**
- **No key fires on the waves:** in all 180 logs face_mmat limited x1 faces = 0 (the Zalesak cap cannot act); implicit_vimp_lag_local lag-path
  triggers 0; time2_adm_poslim admitted stages 0; time2_dtr_poslim (eint <= 0 BE solves) 0. time2_stage_local_be=false changes nothing here either:
  no stage fallbacks occur on the waves.
- So q vs FIN = d6int vs FIN (22b8e580d): the P3.33 tau10 fix holds (vetf p rho / E 1.99 / 2.00); the P3.33 tau0.01 wbsolve flag is unchanged
  (L1 E-mrm 3.564e-4 vs FIN 3.572e-4); the only > 1 % flag vs FIN is P3.33 tau10 N256 L1 rho (+1.2 % vetf, +1.4 % edd). No order flag.

### cmp_p18h_d6i_vetf.md
```
| clos | P | tau | N | g^4 | dev rho | dev E | |dw|/|w| | max dev / g^4 | note |
|---|---|---|---|---|---|---|---|---|---|
| vetf | 0.03333 | 0.01 | 16 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 9.469e-02 W 9.469e-02 |
| vetf | 0.03333 | 0.01 | 32 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 3.781e-02 W 3.781e-02 |
| vetf | 0.03333 | 0.01 | 64 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.059e-02 W 1.059e-02 |
| vetf | 0.03333 | 0.01 | 128 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.778e-03 W 2.778e-03 |
| vetf | 0.03333 | 0.01 | 256 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 6.939e-04 W 6.939e-04 |
| vetf | 0.03333 | 0.1 | 16 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 9.144e-02 W 9.144e-02 |
| vetf | 0.03333 | 0.1 | 32 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 3.616e-02 W 3.616e-02 |
| vetf | 0.03333 | 0.1 | 64 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.004e-02 W 1.004e-02 |
| vetf | 0.03333 | 0.1 | 128 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.519e-03 W 2.519e-03 |
| vetf | 0.03333 | 0.1 | 256 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 6.005e-04 W 6.005e-04 |
| vetf | 0.03333 | 1 | 16 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 6.923e-02 W 6.923e-02 |
| vetf | 0.03333 | 1 | 32 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.384e-02 W 2.384e-02 |
| vetf | 0.03333 | 1 | 64 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 6.330e-03 W 6.330e-03 |
| vetf | 0.03333 | 1 | 128 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.522e-03 W 1.522e-03 |
| vetf | 0.03333 | 1 | 256 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 3.424e-04 W 3.424e-04 |
| vetf | 0.03333 | 10 | 16 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 6.070e-02 W 6.070e-02 |
| vetf | 0.03333 | 10 | 32 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.638e-02 W 1.638e-02 |
| vetf | 0.03333 | 10 | 64 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 4.085e-03 W 4.085e-03 |
| vetf | 0.03333 | 10 | 128 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 9.839e-04 W 9.839e-04 |
| vetf | 0.03333 | 10 | 256 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.360e-04 W 2.360e-04 |
| vetf | 0.03333 | 100 | 16 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 8.425e-02 W 8.425e-02 |
| vetf | 0.03333 | 100 | 32 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.870e-02 W 2.870e-02 |
| vetf | 0.03333 | 100 | 64 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 6.757e-03 W 6.757e-03 |
| vetf | 0.03333 | 100 | 128 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.376e-03 W 1.376e-03 |
| vetf | 0.03333 | 100 | 256 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.870e-04 W 2.870e-04 |
| vetf | 0.03333 | 1000 | 16 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 8.856e-02 W 8.856e-02 |
| vetf | 0.03333 | 1000 | 32 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 3.509e-02 W 3.509e-02 |
| vetf | 0.03333 | 1000 | 64 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 9.449e-03 W 9.449e-03 |
| vetf | 0.03333 | 1000 | 128 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.264e-03 W 2.264e-03 |
| vetf | 0.03333 | 1000 | 256 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 4.888e-04 W 4.888e-04 |
| vetf | 0.3333 | 0.01 | 16 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 9.101e-02 W 9.101e-02 |
| vetf | 0.3333 | 0.01 | 32 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 3.605e-02 W 3.605e-02 |
| vetf | 0.3333 | 0.01 | 64 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 9.990e-03 W 9.990e-03 |
| vetf | 0.3333 | 0.01 | 128 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.502e-03 W 2.502e-03 |
| vetf | 0.3333 | 0.01 | 256 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 5.956e-04 W 5.956e-04 |
| vetf | 0.3333 | 0.1 | 16 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 6.367e-02 W 6.367e-02 |
| vetf | 0.3333 | 0.1 | 32 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.191e-02 W 2.191e-02 |
| vetf | 0.3333 | 0.1 | 64 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 5.623e-03 W 5.623e-03 |
| vetf | 0.3333 | 0.1 | 128 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.330e-03 W 1.330e-03 |
| vetf | 0.3333 | 0.1 | 256 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.956e-04 W 2.956e-04 |
| vetf | 0.3333 | 1 | 16 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 4.032e-02 W 4.032e-02 |
| vetf | 0.3333 | 1 | 32 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 9.650e-03 W 9.650e-03 |
| vetf | 0.3333 | 1 | 64 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.416e-03 W 2.416e-03 |
| vetf | 0.3333 | 1 | 128 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 6.090e-04 W 6.090e-04 |
| vetf | 0.3333 | 1 | 256 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.507e-04 W 1.507e-04 |
| vetf | 0.3333 | 10 | 16 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.448e-02 W 2.448e-02 |
| vetf | 0.3333 | 10 | 32 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 6.135e-03 W 6.135e-03 |
| vetf | 0.3333 | 10 | 64 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.423e-03 W 1.423e-03 |
| vetf | 0.3333 | 10 | 128 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 3.278e-04 W 3.278e-04 |
| vetf | 0.3333 | 10 | 256 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 8.184e-05 W 8.184e-05 |
| vetf | 0.3333 | 100 | 16 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 4.637e-02 W 4.637e-02 |
| vetf | 0.3333 | 100 | 32 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.154e-02 W 1.154e-02 |
| vetf | 0.3333 | 100 | 64 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.019e-03 W 2.019e-03 |
| vetf | 0.3333 | 100 | 128 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 4.069e-04 W 4.069e-04 |
| vetf | 0.3333 | 100 | 256 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 9.210e-05 W 9.210e-05 |
| vetf | 0.3333 | 1000 | 16 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 6.408e-02 W 6.408e-02 |
| vetf | 0.3333 | 1000 | 32 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.237e-02 W 2.237e-02 |
| vetf | 0.3333 | 1000 | 64 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 4.928e-03 W 4.928e-03 |
| vetf | 0.3333 | 1000 | 128 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 9.091e-04 W 9.091e-04 |
| vetf | 0.3333 | 1000 | 256 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.484e-04 W 1.484e-04 |
| vetf | 3.333 | 0.01 | 16 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 6.102e-02 W 6.102e-02 |
| vetf | 3.333 | 0.01 | 32 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.114e-02 W 2.114e-02 |
| vetf | 3.333 | 0.01 | 64 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 5.370e-03 W 5.370e-03 |
| vetf | 3.333 | 0.01 | 128 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.243e-03 W 1.243e-03 |
| vetf | 3.333 | 0.01 | 256 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.720e-04 W 2.720e-04 |
| vetf | 3.333 | 0.1 | 16 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 4.642e-02 W 4.642e-02 |
| vetf | 3.333 | 0.1 | 32 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.074e-02 W 1.074e-02 |
| vetf | 3.333 | 0.1 | 64 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.653e-03 W 2.653e-03 |
| vetf | 3.333 | 0.1 | 128 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 6.664e-04 W 6.664e-04 |
| vetf | 3.333 | 0.1 | 256 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.657e-04 W 1.657e-04 |
| vetf | 3.333 | 1 | 16 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 3.547e-02 W 3.547e-02 |
| vetf | 3.333 | 1 | 32 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.142e-02 W 1.142e-02 |
| vetf | 3.333 | 1 | 64 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.597e-03 W 2.597e-03 |
| vetf | 3.333 | 1 | 128 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 5.649e-04 W 5.649e-04 |
| vetf | 3.333 | 1 | 256 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.301e-04 W 1.301e-04 |
| vetf | 3.333 | 10 | 16 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.183e-02 W 1.183e-02 |
| vetf | 3.333 | 10 | 32 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.896e-03 W 2.896e-03 |
| vetf | 3.333 | 10 | 64 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 7.354e-04 W 7.354e-04 |
| vetf | 3.333 | 10 | 128 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.860e-04 W 1.860e-04 |
| vetf | 3.333 | 10 | 256 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 4.693e-05 W 4.693e-05 |
| vetf | 3.333 | 100 | 16 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 6.138e-02 W 6.138e-02 |
| vetf | 3.333 | 100 | 32 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.534e-02 W 1.534e-02 |
| vetf | 3.333 | 100 | 64 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 3.762e-03 W 3.762e-03 |
| vetf | 3.333 | 100 | 128 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 9.255e-04 W 9.255e-04 |
| vetf | 3.333 | 100 | 256 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.283e-04 W 2.283e-04 |
| vetf | 3.333 | 1000 | 16 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 6.987e-02 W 6.987e-02 |
| vetf | 3.333 | 1000 | 32 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.820e-02 W 1.820e-02 |
| vetf | 3.333 | 1000 | 64 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 4.660e-03 W 4.660e-03 |
| vetf | 3.333 | 1000 | 128 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.169e-03 W 1.169e-03 |
| vetf | 3.333 | 1000 | 256 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.858e-04 W 2.858e-04 |

per (closure, P): max over tau, N of (dev/g^4, dev)
('vetf', 0.03333333333333333) dev/g^4 0  dev 0.000e+00  g^4 3.160e-12
('vetf', 0.3333333333333333) dev/g^4 0  dev 0.000e+00  g^4 3.160e-08
('vetf', 3.3333333333333335) dev/g^4 0  dev 0.000e+00  g^4 3.160e-04
```

### cmp_p18h_d6i_edd.md
```
| clos | P | tau | N | g^4 | dev rho | dev E | |dw|/|w| | max dev / g^4 | note |
|---|---|---|---|---|---|---|---|---|---|
| edd | 0.03333 | 0.01 | 16 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 9.469e-02 W 9.469e-02 |
| edd | 0.03333 | 0.01 | 32 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 3.781e-02 W 3.781e-02 |
| edd | 0.03333 | 0.01 | 64 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.059e-02 W 1.059e-02 |
| edd | 0.03333 | 0.01 | 128 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.778e-03 W 2.778e-03 |
| edd | 0.03333 | 0.01 | 256 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 6.939e-04 W 6.939e-04 |
| edd | 0.03333 | 0.1 | 16 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 9.140e-02 W 9.140e-02 |
| edd | 0.03333 | 0.1 | 32 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 3.614e-02 W 3.614e-02 |
| edd | 0.03333 | 0.1 | 64 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.004e-02 W 1.004e-02 |
| edd | 0.03333 | 0.1 | 128 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.517e-03 W 2.517e-03 |
| edd | 0.03333 | 0.1 | 256 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 5.995e-04 W 5.995e-04 |
| edd | 0.03333 | 1 | 16 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 6.480e-02 W 6.480e-02 |
| edd | 0.03333 | 1 | 32 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.235e-02 W 2.235e-02 |
| edd | 0.03333 | 1 | 64 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 5.843e-03 W 5.843e-03 |
| edd | 0.03333 | 1 | 128 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.397e-03 W 1.397e-03 |
| edd | 0.03333 | 1 | 256 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 3.157e-04 W 3.157e-04 |
| edd | 0.03333 | 10 | 16 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 5.848e-02 W 5.848e-02 |
| edd | 0.03333 | 10 | 32 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.494e-02 W 1.494e-02 |
| edd | 0.03333 | 10 | 64 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 3.780e-03 W 3.780e-03 |
| edd | 0.03333 | 10 | 128 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 9.240e-04 W 9.240e-04 |
| edd | 0.03333 | 10 | 256 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.224e-04 W 2.224e-04 |
| edd | 0.03333 | 100 | 16 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 8.429e-02 W 8.429e-02 |
| edd | 0.03333 | 100 | 32 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.854e-02 W 2.854e-02 |
| edd | 0.03333 | 100 | 64 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 6.644e-03 W 6.644e-03 |
| edd | 0.03333 | 100 | 128 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.335e-03 W 1.335e-03 |
| edd | 0.03333 | 100 | 256 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.792e-04 W 2.792e-04 |
| edd | 0.03333 | 1000 | 16 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 8.858e-02 W 8.858e-02 |
| edd | 0.03333 | 1000 | 32 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 3.509e-02 W 3.509e-02 |
| edd | 0.03333 | 1000 | 64 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 9.451e-03 W 9.451e-03 |
| edd | 0.03333 | 1000 | 128 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.265e-03 W 2.265e-03 |
| edd | 0.03333 | 1000 | 256 | 3.16e-12 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 4.886e-04 W 4.886e-04 |
| edd | 0.3333 | 0.01 | 16 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 9.101e-02 W 9.101e-02 |
| edd | 0.3333 | 0.01 | 32 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 3.605e-02 W 3.605e-02 |
| edd | 0.3333 | 0.01 | 64 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 9.990e-03 W 9.990e-03 |
| edd | 0.3333 | 0.01 | 128 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.502e-03 W 2.502e-03 |
| edd | 0.3333 | 0.01 | 256 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 5.957e-04 W 5.957e-04 |
| edd | 0.3333 | 0.1 | 16 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 6.348e-02 W 6.348e-02 |
| edd | 0.3333 | 0.1 | 32 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.183e-02 W 2.183e-02 |
| edd | 0.3333 | 0.1 | 64 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 5.597e-03 W 5.597e-03 |
| edd | 0.3333 | 0.1 | 128 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.324e-03 W 1.324e-03 |
| edd | 0.3333 | 0.1 | 256 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.932e-04 W 2.932e-04 |
| edd | 0.3333 | 1 | 16 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 4.352e-02 W 4.352e-02 |
| edd | 0.3333 | 1 | 32 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.032e-02 W 1.032e-02 |
| edd | 0.3333 | 1 | 64 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.618e-03 W 2.618e-03 |
| edd | 0.3333 | 1 | 128 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 6.621e-04 W 6.621e-04 |
| edd | 0.3333 | 1 | 256 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.649e-04 W 1.649e-04 |
| edd | 0.3333 | 10 | 16 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.370e-02 W 2.370e-02 |
| edd | 0.3333 | 10 | 32 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 6.023e-03 W 6.023e-03 |
| edd | 0.3333 | 10 | 64 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.391e-03 W 1.391e-03 |
| edd | 0.3333 | 10 | 128 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 3.160e-04 W 3.160e-04 |
| edd | 0.3333 | 10 | 256 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 7.808e-05 W 7.808e-05 |
| edd | 0.3333 | 100 | 16 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 4.659e-02 W 4.659e-02 |
| edd | 0.3333 | 100 | 32 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.148e-02 W 1.148e-02 |
| edd | 0.3333 | 100 | 64 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.028e-03 W 2.028e-03 |
| edd | 0.3333 | 100 | 128 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 4.143e-04 W 4.143e-04 |
| edd | 0.3333 | 100 | 256 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 9.294e-05 W 9.294e-05 |
| edd | 0.3333 | 1000 | 16 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 6.416e-02 W 6.416e-02 |
| edd | 0.3333 | 1000 | 32 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.240e-02 W 2.240e-02 |
| edd | 0.3333 | 1000 | 64 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 4.955e-03 W 4.955e-03 |
| edd | 0.3333 | 1000 | 128 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 9.197e-04 W 9.197e-04 |
| edd | 0.3333 | 1000 | 256 | 3.16e-08 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.524e-04 W 1.524e-04 |
| edd | 3.333 | 0.01 | 16 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 6.102e-02 W 6.102e-02 |
| edd | 3.333 | 0.01 | 32 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.114e-02 W 2.114e-02 |
| edd | 3.333 | 0.01 | 64 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 5.369e-03 W 5.369e-03 |
| edd | 3.333 | 0.01 | 128 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.243e-03 W 1.243e-03 |
| edd | 3.333 | 0.01 | 256 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.719e-04 W 2.719e-04 |
| edd | 3.333 | 0.1 | 16 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 4.656e-02 W 4.656e-02 |
| edd | 3.333 | 0.1 | 32 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.078e-02 W 1.078e-02 |
| edd | 3.333 | 0.1 | 64 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.662e-03 W 2.662e-03 |
| edd | 3.333 | 0.1 | 128 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 6.692e-04 W 6.692e-04 |
| edd | 3.333 | 0.1 | 256 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.665e-04 W 1.665e-04 |
| edd | 3.333 | 1 | 16 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 3.615e-02 W 3.615e-02 |
| edd | 3.333 | 1 | 32 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.155e-02 W 1.155e-02 |
| edd | 3.333 | 1 | 64 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.641e-03 W 2.641e-03 |
| edd | 3.333 | 1 | 128 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 5.782e-04 W 5.782e-04 |
| edd | 3.333 | 1 | 256 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.337e-04 W 1.337e-04 |
| edd | 3.333 | 10 | 16 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 9.008e-03 W 9.008e-03 |
| edd | 3.333 | 10 | 32 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.264e-03 W 2.264e-03 |
| edd | 3.333 | 10 | 64 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 5.919e-04 W 5.919e-04 |
| edd | 3.333 | 10 | 128 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.526e-04 W 1.526e-04 |
| edd | 3.333 | 10 | 256 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 3.948e-05 W 3.948e-05 |
| edd | 3.333 | 100 | 16 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 5.903e-02 W 5.903e-02 |
| edd | 3.333 | 100 | 32 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.461e-02 W 1.461e-02 |
| edd | 3.333 | 100 | 64 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 3.602e-03 W 3.602e-03 |
| edd | 3.333 | 100 | 128 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 8.974e-04 W 8.974e-04 |
| edd | 3.333 | 100 | 256 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.244e-04 W 2.244e-04 |
| edd | 3.333 | 1000 | 16 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 6.949e-02 W 6.949e-02 |
| edd | 3.333 | 1000 | 32 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.801e-02 W 1.801e-02 |
| edd | 3.333 | 1000 | 64 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 4.578e-03 W 4.578e-03 |
| edd | 3.333 | 1000 | 128 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 1.130e-03 W 1.130e-03 |
| edd | 3.333 | 1000 | 256 | 3.16e-04 | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0 | L1rho B 2.719e-04 W 2.719e-04 |

per (closure, P): max over tau, N of (dev/g^4, dev)
('edd', 0.03333333333333333) dev/g^4 0  dev 0.000e+00  g^4 3.160e-12
('edd', 0.3333333333333333) dev/g^4 0  dev 0.000e+00  g^4 3.160e-08
('edd', 3.3333333333333335) dev/g^4 0  dev 0.000e+00  g^4 3.160e-04
```

### cmp_p18h_fin_vetf.md
```
| clos | P | tau | N | g^4 | dev rho | dev E | |dw|/|w| | max dev / g^4 | note |
|---|---|---|---|---|---|---|---|---|---|
| vetf | 0.03333 | 0.01 | 16 | 3.16e-12 | 8.88e-10 | 3.14e-06 | 2.56e-11 | 9.9e+05 | L1rho B 9.469e-02 W 9.469e-02 |
| vetf | 0.03333 | 0.01 | 32 | 3.16e-12 | 2.00e-09 | 1.65e-06 | 5.27e-11 | 5.2e+05 | L1rho B 3.781e-02 W 3.781e-02 |
| vetf | 0.03333 | 0.01 | 64 | 3.16e-12 | 2.22e-09 | 2.07e-06 | 1.64e-11 | 6.6e+05 | L1rho B 1.059e-02 W 1.059e-02 |
| vetf | 0.03333 | 0.01 | 128 | 3.16e-12 | 6.11e-09 | 3.03e-06 | 9.96e-11 | 9.6e+05 | L1rho B 2.778e-03 W 2.778e-03 |
| vetf | 0.03333 | 0.01 | 256 | 3.16e-12 | 7.77e-09 | 5.05e-06 | 1.21e-10 | 1.6e+06 | L1rho B 6.939e-04 W 6.939e-04 |
| vetf | 0.03333 | 0.1 | 16 | 3.16e-12 | 1.11e-09 | 5.90e-08 | 5.65e-11 | 1.9e+04 | L1rho B 9.144e-02 W 9.144e-02 |
| vetf | 0.03333 | 0.1 | 32 | 3.16e-12 | 1.78e-09 | 7.38e-08 | 4.64e-11 | 2.3e+04 | L1rho B 3.616e-02 W 3.616e-02 |
| vetf | 0.03333 | 0.1 | 64 | 3.16e-12 | 3.11e-09 | 1.48e-07 | 6.95e-11 | 4.7e+04 | L1rho B 1.004e-02 W 1.004e-02 |
| vetf | 0.03333 | 0.1 | 128 | 3.16e-12 | 5.55e-09 | 1.44e-07 | 8.93e-12 | 4.6e+04 | L1rho B 2.519e-03 W 2.519e-03 |
| vetf | 0.03333 | 0.1 | 256 | 3.16e-12 | 7.55e-09 | 2.62e-07 | 3.40e-11 | 8.3e+04 | L1rho B 6.005e-04 W 6.005e-04 |
| vetf | 0.03333 | 1 | 16 | 3.16e-12 | 8.88e-10 | 2.20e-09 | 1.23e-10 | 7e+02 | L1rho B 6.923e-02 W 6.923e-02 |
| vetf | 0.03333 | 1 | 32 | 3.16e-12 | 1.33e-09 | 2.80e-09 | 5.92e-11 | 8.9e+02 | L1rho B 2.384e-02 W 2.384e-02 |
| vetf | 0.03333 | 1 | 64 | 3.16e-12 | 3.11e-09 | 3.60e-09 | 4.38e-11 | 1.1e+03 | L1rho B 6.330e-03 W 6.330e-03 |
| vetf | 0.03333 | 1 | 128 | 3.16e-12 | 4.44e-09 | 6.00e-09 | 4.64e-11 | 1.9e+03 | L1rho B 1.522e-03 W 1.522e-03 |
| vetf | 0.03333 | 1 | 256 | 3.16e-12 | 5.88e-09 | 6.20e-09 | 5.60e-11 | 2e+03 | L1rho B 3.424e-04 W 3.424e-04 |
| vetf | 0.03333 | 10 | 16 | 3.16e-12 | 4.44e-10 | 7.52e-10 | 6.55e-11 | 2.4e+02 | L1rho B 6.070e-02 W 6.070e-02 |
| vetf | 0.03333 | 10 | 32 | 3.16e-12 | 7.77e-10 | 1.38e-09 | 2.52e-11 | 4.4e+02 | L1rho B 1.638e-02 W 1.638e-02 |
| vetf | 0.03333 | 10 | 64 | 3.16e-12 | 2.44e-09 | 9.84e-09 | 4.34e-10 | 3.1e+03 | L1rho B 4.085e-03 W 4.085e-03 |
| vetf | 0.03333 | 10 | 128 | 3.16e-12 | 1.33e-09 | 2.13e-09 | 5.16e-11 | 6.7e+02 | L1rho B 9.839e-04 W 9.839e-04 |
| vetf | 0.03333 | 10 | 256 | 3.16e-12 | 2.89e-09 | 1.38e-09 | 4.82e-11 | 9.1e+02 | L1rho B 2.360e-04 W 2.360e-04 |
| vetf | 0.03333 | 100 | 16 | 3.16e-12 | 2.33e-09 | 7.57e-09 | 1.53e-10 | 2.4e+03 | L1rho B 8.425e-02 W 8.425e-02 |
| vetf | 0.03333 | 100 | 32 | 3.16e-12 | 1.44e-09 | 1.80e-09 | 1.13e-10 | 5.7e+02 | L1rho B 2.870e-02 W 2.870e-02 |
| vetf | 0.03333 | 100 | 64 | 3.16e-12 | 1.17e-08 | 1.05e-08 | 1.71e-11 | 3.7e+03 | L1rho B 6.757e-03 W 6.757e-03 |
| vetf | 0.03333 | 100 | 128 | 3.16e-12 | 5.22e-09 | 2.71e-08 | 1.50e-10 | 8.6e+03 | L1rho B 1.376e-03 W 1.376e-03 |
| vetf | 0.03333 | 100 | 256 | 3.16e-12 | 1.78e-09 | 3.73e-09 | 8.61e-11 | 1.2e+03 | L1rho B 2.870e-04 W 2.870e-04 |
| vetf | 0.03333 | 1000 | 16 | 3.16e-12 | 3.55e-09 | 1.09e-08 | 9.00e-11 | 3.5e+03 | L1rho B 8.856e-02 W 8.856e-02 |
| vetf | 0.03333 | 1000 | 32 | 3.16e-12 | 5.77e-09 | 6.22e-09 | 2.59e-10 | 2e+03 | L1rho B 3.509e-02 W 3.509e-02 |
| vetf | 0.03333 | 1000 | 64 | 3.16e-12 | 7.16e-08 | 9.11e-08 | 3.81e-10 | 2.9e+04 | L1rho B 9.449e-03 W 9.449e-03 |
| vetf | 0.03333 | 1000 | 128 | 3.16e-12 | 2.35e-08 | 2.65e-08 | 6.61e-11 | 8.4e+03 | L1rho B 2.264e-03 W 2.264e-03 |
| vetf | 0.03333 | 1000 | 256 | 3.16e-12 | 5.51e-08 | 5.72e-08 | 1.13e-10 | 1.8e+04 | L1rho B 4.888e-04 W 4.888e-04 |
| vetf | 0.3333 | 0.01 | 16 | 3.16e-08 | 8.88e-10 | 1.21e-05 | 4.14e-11 | 3.8e+02 | L1rho B 9.101e-02 W 9.101e-02 |
| vetf | 0.3333 | 0.01 | 32 | 3.16e-08 | 1.89e-09 | 2.21e-06 | 1.08e-10 | 70 | L1rho B 3.605e-02 W 3.605e-02 |
| vetf | 0.3333 | 0.01 | 64 | 3.16e-08 | 2.66e-09 | 2.13e-06 | 3.20e-11 | 67 | L1rho B 9.990e-03 W 9.990e-03 |
| vetf | 0.3333 | 0.01 | 128 | 3.16e-08 | 4.66e-09 | 3.06e-06 | 6.81e-11 | 97 | L1rho B 2.502e-03 W 2.502e-03 |
| vetf | 0.3333 | 0.01 | 256 | 3.16e-08 | 7.33e-09 | 3.91e-06 | 8.98e-11 | 1.2e+02 | L1rho B 5.956e-04 W 5.956e-04 |
| vetf | 0.3333 | 0.1 | 16 | 3.16e-08 | 1.44e-09 | 2.64e-07 | 7.93e-11 | 8.4 | L1rho B 6.367e-02 W 6.367e-02 |
| vetf | 0.3333 | 0.1 | 32 | 3.16e-08 | 2.66e-09 | 1.89e-07 | 5.62e-11 | 6 | L1rho B 2.191e-02 W 2.191e-02 |
| vetf | 0.3333 | 0.1 | 64 | 3.16e-08 | 2.44e-09 | 1.11e-07 | 6.78e-11 | 3.5 | L1rho B 5.623e-03 W 5.623e-03 |
| vetf | 0.3333 | 0.1 | 128 | 3.16e-08 | 3.33e-09 | 2.52e-07 | 1.47e-11 | 8 | L1rho B 1.330e-03 W 1.330e-03 |
| vetf | 0.3333 | 0.1 | 256 | 3.16e-08 | 5.77e-09 | 2.10e-07 | 1.06e-10 | 6.6 | L1rho B 2.956e-04 W 2.956e-04 |
| vetf | 0.3333 | 1 | 16 | 3.16e-08 | 2.44e-09 | 1.24e-08 | 9.72e-10 | 0.39 | L1rho B 4.032e-02 W 4.032e-02 |
| vetf | 0.3333 | 1 | 32 | 3.16e-08 | 2.22e-09 | 7.75e-09 | 6.26e-10 | 0.25 | L1rho B 9.650e-03 W 9.650e-03 |
| vetf | 0.3333 | 1 | 64 | 3.16e-08 | 1.78e-09 | 1.67e-08 | 3.83e-10 | 0.53 | L1rho B 2.416e-03 W 2.416e-03 |
| vetf | 0.3333 | 1 | 128 | 3.16e-08 | 1.78e-09 | 1.63e-08 | 1.72e-10 | 0.52 | L1rho B 6.090e-04 W 6.090e-04 |
| vetf | 0.3333 | 1 | 256 | 3.16e-08 | 2.00e-09 | 1.47e-08 | 1.17e-10 | 0.47 | L1rho B 1.507e-04 W 1.507e-04 |
| vetf | 0.3333 | 10 | 16 | 3.16e-08 | 5.77e-09 | 1.03e-08 | 3.65e-09 | 0.33 | L1rho B 2.448e-02 W 2.448e-02 |
| vetf | 0.3333 | 10 | 32 | 3.16e-08 | 2.33e-09 | 4.98e-09 | 1.32e-09 | 0.16 | L1rho B 6.135e-03 W 6.135e-03 |
| vetf | 0.3333 | 10 | 64 | 3.16e-08 | 3.11e-09 | 3.36e-08 | 1.43e-09 | 1.1 | L1rho B 1.423e-03 W 1.423e-03 |
| vetf | 0.3333 | 10 | 128 | 3.16e-08 | 3.77e-09 | 7.01e-09 | 1.91e-09 | 0.22 | L1rho B 3.278e-04 W 3.278e-04 |
| vetf | 0.3333 | 10 | 256 | 3.16e-08 | 2.00e-09 | 2.53e-09 | 5.22e-10 | 0.08 | L1rho B 8.184e-05 W 8.184e-05 |
| vetf | 0.3333 | 100 | 16 | 3.16e-08 | 2.63e-08 | 2.12e-08 | 5.06e-09 | 0.83 | L1rho B 4.637e-02 W 4.637e-02 |
| vetf | 0.3333 | 100 | 32 | 3.16e-08 | 3.33e-08 | 3.72e-08 | 3.06e-09 | 1.2 | L1rho B 1.154e-02 W 1.154e-02 |
| vetf | 0.3333 | 100 | 64 | 3.16e-08 | 2.22e-09 | 2.75e-09 | 2.82e-10 | 0.087 | L1rho B 2.019e-03 W 2.019e-03 |
| vetf | 0.3333 | 100 | 128 | 3.16e-08 | 9.10e-09 | 4.00e-08 | 3.30e-10 | 1.3 | L1rho B 4.069e-04 W 4.069e-04 |
| vetf | 0.3333 | 100 | 256 | 3.16e-08 | 5.33e-09 | 9.25e-09 | 7.21e-10 | 0.29 | L1rho B 9.210e-05 W 9.210e-05 |
| vetf | 0.3333 | 1000 | 16 | 3.16e-08 | 6.51e-08 | 5.53e-08 | 1.07e-08 | 2.1 | L1rho B 6.408e-02 W 6.408e-02 |
| vetf | 0.3333 | 1000 | 32 | 3.16e-08 | 4.84e-08 | 5.17e-08 | 3.67e-09 | 1.6 | L1rho B 2.237e-02 W 2.237e-02 |
| vetf | 0.3333 | 1000 | 64 | 3.16e-08 | 6.51e-08 | 5.42e-08 | 9.79e-10 | 2.1 | L1rho B 4.928e-03 W 4.928e-03 |
| vetf | 0.3333 | 1000 | 128 | 3.16e-08 | 1.20e-08 | 1.15e-08 | 1.40e-10 | 0.38 | L1rho B 9.091e-04 W 9.091e-04 |
| vetf | 0.3333 | 1000 | 256 | 3.16e-08 | 3.24e-08 | 3.78e-08 | 1.37e-10 | 1.2 | L1rho B 1.484e-04 W 1.484e-04 |
| vetf | 3.333 | 0.01 | 16 | 3.16e-04 | 1.05e-07 | 7.10e-06 | 2.75e-08 | 0.022 | L1rho B 6.102e-02 W 6.102e-02 |
| vetf | 3.333 | 0.01 | 32 | 3.16e-04 | 6.53e-08 | 1.14e-05 | 1.60e-08 | 0.036 | L1rho B 2.114e-02 W 2.114e-02 |
| vetf | 3.333 | 0.01 | 64 | 3.16e-04 | 3.42e-08 | 6.25e-06 | 8.34e-09 | 0.02 | L1rho B 5.370e-03 W 5.370e-03 |
| vetf | 3.333 | 0.01 | 128 | 3.16e-04 | 2.02e-08 | 4.97e-06 | 4.21e-09 | 0.016 | L1rho B 1.243e-03 W 1.243e-03 |
| vetf | 3.333 | 0.01 | 256 | 3.16e-04 | 1.20e-08 | 4.40e-06 | 2.10e-09 | 0.014 | L1rho B 2.720e-04 W 2.720e-04 |
| vetf | 3.333 | 0.1 | 16 | 3.16e-04 | 1.12e-06 | 1.75e-06 | 3.65e-07 | 0.0055 | L1rho B 4.642e-02 W 4.642e-02 |
| vetf | 3.333 | 0.1 | 32 | 3.16e-04 | 7.48e-07 | 7.27e-06 | 2.22e-07 | 0.023 | L1rho B 1.073e-02 W 1.074e-02 |
| vetf | 3.333 | 0.1 | 64 | 3.16e-04 | 3.96e-07 | 4.34e-06 | 1.17e-07 | 0.014 | L1rho B 2.653e-03 W 2.653e-03 |
| vetf | 3.333 | 0.1 | 128 | 3.16e-04 | 2.02e-07 | 2.38e-06 | 5.93e-08 | 0.0075 | L1rho B 6.662e-04 W 6.664e-04 |
| vetf | 3.333 | 0.1 | 256 | 3.16e-04 | 1.02e-07 | 1.50e-06 | 2.97e-08 | 0.0048 | L1rho B 1.656e-04 W 1.657e-04 |
| vetf | 3.333 | 1 | 16 | 3.16e-04 | 6.54e-06 | 2.21e-06 | 2.54e-06 | 0.021 | L1rho B 3.546e-02 W 3.547e-02 |
| vetf | 3.333 | 1 | 32 | 3.16e-04 | 4.89e-06 | 9.89e-07 | 1.89e-06 | 0.015 | L1rho B 1.142e-02 W 1.142e-02 |
| vetf | 3.333 | 1 | 64 | 3.16e-04 | 2.78e-06 | 5.74e-07 | 1.10e-06 | 0.0088 | L1rho B 2.595e-03 W 2.597e-03 |
| vetf | 3.333 | 1 | 128 | 3.16e-04 | 1.48e-06 | 2.33e-07 | 5.85e-07 | 0.0047 | L1rho B 5.641e-04 W 5.649e-04 |
| vetf | 3.333 | 1 | 256 | 3.16e-04 | 7.64e-07 | 1.36e-07 | 3.01e-07 | 0.0024 | L1rho B 1.297e-04 W 1.301e-04 |
| vetf | 3.333 | 10 | 16 | 3.16e-04 | 1.14e-06 | 6.05e-07 | 1.37e-06 | 0.0036 | L1rho B 1.183e-02 W 1.183e-02 |
| vetf | 3.333 | 10 | 32 | 3.16e-04 | 2.63e-06 | 1.96e-06 | 3.43e-06 | 0.0083 | L1rho B 2.894e-03 W 2.896e-03 |
| vetf | 3.333 | 10 | 64 | 3.16e-04 | 2.31e-06 | 1.92e-06 | 3.11e-06 | 0.0073 | L1rho B 7.340e-04 W 7.354e-04 |
| vetf | 3.333 | 10 | 128 | 3.16e-04 | 1.52e-06 | 1.31e-06 | 2.05e-06 | 0.0048 | L1rho B 1.851e-04 W 1.860e-04 |
| vetf | 3.333 | 10 | 256 | 3.16e-04 | 8.56e-07 | 7.57e-07 | 1.17e-06 | 0.0027 | L1rho B 4.638e-05 W 4.693e-05 |
| vetf | 3.333 | 100 | 16 | 3.16e-04 | 2.45e-06 | 1.91e-06 | 4.87e-07 | 0.0078 | L1rho B 6.138e-02 W 6.138e-02 |
| vetf | 3.333 | 100 | 32 | 3.16e-04 | 4.68e-07 | 4.34e-07 | 1.22e-07 | 0.0015 | L1rho B 1.534e-02 W 1.534e-02 |
| vetf | 3.333 | 100 | 64 | 3.16e-04 | 1.08e-08 | 1.72e-08 | 6.15e-11 | 5.5e-05 | L1rho B 3.762e-03 W 3.762e-03 |
| vetf | 3.333 | 100 | 128 | 3.16e-04 | 3.31e-08 | 3.16e-08 | 9.69e-09 | 0.0001 | L1rho B 9.255e-04 W 9.255e-04 |
| vetf | 3.333 | 100 | 256 | 3.16e-04 | 9.00e-08 | 8.53e-08 | 2.59e-08 | 0.00028 | L1rho B 2.283e-04 W 2.283e-04 |
| vetf | 3.333 | 1000 | 16 | 3.16e-04 | 1.27e-05 | 9.50e-06 | 1.79e-06 | 0.04 | L1rho B 6.988e-02 W 6.987e-02 |
| vetf | 3.333 | 1000 | 32 | 3.16e-04 | 8.32e-06 | 7.56e-06 | 1.26e-06 | 0.026 | L1rho B 1.820e-02 W 1.820e-02 |
| vetf | 3.333 | 1000 | 64 | 3.16e-04 | 1.42e-06 | 1.40e-06 | 3.30e-07 | 0.0045 | L1rho B 4.660e-03 W 4.660e-03 |
| vetf | 3.333 | 1000 | 128 | 3.16e-04 | 1.96e-07 | 1.93e-07 | 4.82e-08 | 0.00062 | L1rho B 1.169e-03 W 1.169e-03 |
| vetf | 3.333 | 1000 | 256 | 3.16e-04 | 5.13e-08 | 4.53e-08 | 3.90e-09 | 0.00016 | L1rho B 2.859e-04 W 2.858e-04 |

per (closure, P): max over tau, N of (dev/g^4, dev)
('vetf', 0.03333333333333333) dev/g^4 1.6e+06  dev 5.049e-06  g^4 3.160e-12
('vetf', 0.3333333333333333) dev/g^4 383  dev 1.212e-05  g^4 3.160e-08
('vetf', 3.3333333333333335) dev/g^4 0.0403  dev 1.272e-05  g^4 3.160e-04
```

### cmp_p18h_fin_edd.md
```
| clos | P | tau | N | g^4 | dev rho | dev E | |dw|/|w| | max dev / g^4 | note |
|---|---|---|---|---|---|---|---|---|---|
| edd | 0.03333 | 0.01 | 16 | 3.16e-12 | 1.11e-09 | 1.55e-06 | 2.10e-11 | 4.9e+05 | L1rho B 9.469e-02 W 9.469e-02 |
| edd | 0.03333 | 0.01 | 32 | 3.16e-12 | 2.89e-09 | 1.02e-06 | 4.37e-11 | 3.2e+05 | L1rho B 3.781e-02 W 3.781e-02 |
| edd | 0.03333 | 0.01 | 64 | 3.16e-12 | 2.66e-09 | 1.87e-06 | 7.62e-11 | 5.9e+05 | L1rho B 1.059e-02 W 1.059e-02 |
| edd | 0.03333 | 0.01 | 128 | 3.16e-12 | 6.00e-09 | 2.89e-06 | 5.94e-11 | 9.1e+05 | L1rho B 2.778e-03 W 2.778e-03 |
| edd | 0.03333 | 0.01 | 256 | 3.16e-12 | 7.55e-09 | 3.90e-06 | 4.19e-11 | 1.2e+06 | L1rho B 6.939e-04 W 6.939e-04 |
| edd | 0.03333 | 0.1 | 16 | 3.16e-12 | 1.11e-09 | 1.23e-07 | 4.49e-11 | 3.9e+04 | L1rho B 9.140e-02 W 9.140e-02 |
| edd | 0.03333 | 0.1 | 32 | 3.16e-12 | 2.89e-09 | 1.39e-07 | 1.26e-11 | 4.4e+04 | L1rho B 3.614e-02 W 3.614e-02 |
| edd | 0.03333 | 0.1 | 64 | 3.16e-12 | 3.33e-09 | 1.65e-07 | 1.19e-10 | 5.2e+04 | L1rho B 1.004e-02 W 1.004e-02 |
| edd | 0.03333 | 0.1 | 128 | 3.16e-12 | 6.00e-09 | 2.45e-07 | 7.14e-11 | 7.8e+04 | L1rho B 2.517e-03 W 2.517e-03 |
| edd | 0.03333 | 0.1 | 256 | 3.16e-12 | 7.55e-09 | 2.72e-07 | 3.84e-11 | 8.6e+04 | L1rho B 5.995e-04 W 5.995e-04 |
| edd | 0.03333 | 1 | 16 | 3.16e-12 | 1.33e-09 | 7.07e-09 | 1.04e-10 | 2.2e+03 | L1rho B 6.480e-02 W 6.480e-02 |
| edd | 0.03333 | 1 | 32 | 3.16e-12 | 1.78e-09 | 4.86e-09 | 4.78e-11 | 1.5e+03 | L1rho B 2.235e-02 W 2.235e-02 |
| edd | 0.03333 | 1 | 64 | 3.16e-12 | 2.22e-09 | 7.07e-09 | 3.59e-11 | 2.2e+03 | L1rho B 5.843e-03 W 5.843e-03 |
| edd | 0.03333 | 1 | 128 | 3.16e-12 | 3.11e-09 | 7.51e-09 | 2.82e-11 | 2.4e+03 | L1rho B 1.397e-03 W 1.397e-03 |
| edd | 0.03333 | 1 | 256 | 3.16e-12 | 5.77e-09 | 1.24e-08 | 7.34e-11 | 3.9e+03 | L1rho B 3.157e-04 W 3.157e-04 |
| edd | 0.03333 | 10 | 16 | 3.16e-12 | 5.55e-10 | 1.71e-09 | 6.86e-11 | 5.4e+02 | L1rho B 5.848e-02 W 5.848e-02 |
| edd | 0.03333 | 10 | 32 | 3.16e-12 | 1.11e-09 | 1.08e-09 | 4.87e-11 | 3.5e+02 | L1rho B 1.494e-02 W 1.494e-02 |
| edd | 0.03333 | 10 | 64 | 3.16e-12 | 2.78e-09 | 1.14e-08 | 5.08e-10 | 3.6e+03 | L1rho B 3.780e-03 W 3.780e-03 |
| edd | 0.03333 | 10 | 128 | 3.16e-12 | 1.78e-09 | 1.46e-09 | 5.58e-11 | 5.6e+02 | L1rho B 9.240e-04 W 9.240e-04 |
| edd | 0.03333 | 10 | 256 | 3.16e-12 | 2.44e-09 | 9.49e-10 | 1.17e-11 | 7.7e+02 | L1rho B 2.224e-04 W 2.224e-04 |
| edd | 0.03333 | 100 | 16 | 3.16e-12 | 2.11e-09 | 8.09e-09 | 1.17e-10 | 2.6e+03 | L1rho B 8.429e-02 W 8.429e-02 |
| edd | 0.03333 | 100 | 32 | 3.16e-12 | 1.33e-09 | 2.33e-09 | 7.17e-12 | 7.4e+02 | L1rho B 2.854e-02 W 2.854e-02 |
| edd | 0.03333 | 100 | 64 | 3.16e-12 | 1.71e-08 | 4.18e-08 | 2.01e-10 | 1.3e+04 | L1rho B 6.644e-03 W 6.644e-03 |
| edd | 0.03333 | 100 | 128 | 3.16e-12 | 1.44e-08 | 4.76e-08 | 3.79e-10 | 1.5e+04 | L1rho B 1.335e-03 W 1.335e-03 |
| edd | 0.03333 | 100 | 256 | 3.16e-12 | 1.55e-09 | 5.06e-09 | 4.07e-11 | 1.6e+03 | L1rho B 2.792e-04 W 2.792e-04 |
| edd | 0.03333 | 1000 | 16 | 3.16e-12 | 2.11e-09 | 1.06e-08 | 5.65e-11 | 3.4e+03 | L1rho B 8.858e-02 W 8.858e-02 |
| edd | 0.03333 | 1000 | 32 | 3.16e-12 | 9.33e-09 | 6.51e-09 | 4.34e-10 | 3e+03 | L1rho B 3.509e-02 W 3.509e-02 |
| edd | 0.03333 | 1000 | 64 | 3.16e-12 | 6.67e-08 | 8.37e-08 | 3.50e-10 | 2.6e+04 | L1rho B 9.451e-03 W 9.451e-03 |
| edd | 0.03333 | 1000 | 128 | 3.16e-12 | 2.29e-08 | 1.92e-08 | 7.71e-11 | 7.2e+03 | L1rho B 2.265e-03 W 2.265e-03 |
| edd | 0.03333 | 1000 | 256 | 3.16e-12 | 5.35e-08 | 5.96e-08 | 1.29e-10 | 1.9e+04 | L1rho B 4.886e-04 W 4.886e-04 |
| edd | 0.3333 | 0.01 | 16 | 3.16e-08 | 1.11e-09 | 4.88e-06 | 3.74e-11 | 1.5e+02 | L1rho B 9.101e-02 W 9.101e-02 |
| edd | 0.3333 | 0.01 | 32 | 3.16e-08 | 2.44e-09 | 1.71e-06 | 6.98e-11 | 54 | L1rho B 3.605e-02 W 3.605e-02 |
| edd | 0.3333 | 0.01 | 64 | 3.16e-08 | 3.55e-09 | 2.57e-06 | 4.42e-11 | 81 | L1rho B 9.990e-03 W 9.990e-03 |
| edd | 0.3333 | 0.01 | 128 | 3.16e-08 | 4.77e-09 | 2.14e-06 | 5.67e-11 | 68 | L1rho B 2.502e-03 W 2.502e-03 |
| edd | 0.3333 | 0.01 | 256 | 3.16e-08 | 7.33e-09 | 4.19e-06 | 7.39e-11 | 1.3e+02 | L1rho B 5.957e-04 W 5.957e-04 |
| edd | 0.3333 | 0.1 | 16 | 3.16e-08 | 1.33e-09 | 4.44e-08 | 5.61e-11 | 1.4 | L1rho B 6.348e-02 W 6.348e-02 |
| edd | 0.3333 | 0.1 | 32 | 3.16e-08 | 1.67e-09 | 2.04e-07 | 6.10e-11 | 6.5 | L1rho B 2.183e-02 W 2.183e-02 |
| edd | 0.3333 | 0.1 | 64 | 3.16e-08 | 2.11e-09 | 1.51e-07 | 3.66e-11 | 4.8 | L1rho B 5.597e-03 W 5.597e-03 |
| edd | 0.3333 | 0.1 | 128 | 3.16e-08 | 3.55e-09 | 2.13e-07 | 1.90e-11 | 6.7 | L1rho B 1.324e-03 W 1.324e-03 |
| edd | 0.3333 | 0.1 | 256 | 3.16e-08 | 4.88e-09 | 2.67e-07 | 2.85e-11 | 8.4 | L1rho B 2.932e-04 W 2.932e-04 |
| edd | 0.3333 | 1 | 16 | 3.16e-08 | 2.55e-09 | 3.95e-08 | 8.33e-10 | 1.2 | L1rho B 4.352e-02 W 4.352e-02 |
| edd | 0.3333 | 1 | 32 | 3.16e-08 | 2.44e-09 | 2.69e-08 | 6.15e-10 | 0.85 | L1rho B 1.032e-02 W 1.032e-02 |
| edd | 0.3333 | 1 | 64 | 3.16e-08 | 1.89e-09 | 3.59e-08 | 4.02e-10 | 1.1 | L1rho B 2.618e-03 W 2.618e-03 |
| edd | 0.3333 | 1 | 128 | 3.16e-08 | 1.78e-09 | 3.50e-08 | 1.96e-10 | 1.1 | L1rho B 6.621e-04 W 6.621e-04 |
| edd | 0.3333 | 1 | 256 | 3.16e-08 | 2.66e-09 | 4.31e-08 | 6.78e-11 | 1.4 | L1rho B 1.649e-04 W 1.649e-04 |
| edd | 0.3333 | 10 | 16 | 3.16e-08 | 6.44e-09 | 1.38e-08 | 4.19e-09 | 0.44 | L1rho B 2.370e-02 W 2.370e-02 |
| edd | 0.3333 | 10 | 32 | 3.16e-08 | 2.44e-09 | 4.08e-09 | 1.34e-09 | 0.13 | L1rho B 6.023e-03 W 6.023e-03 |
| edd | 0.3333 | 10 | 64 | 3.16e-08 | 2.44e-09 | 3.23e-08 | 1.32e-09 | 1 | L1rho B 1.391e-03 W 1.391e-03 |
| edd | 0.3333 | 10 | 128 | 3.16e-08 | 3.77e-09 | 6.93e-09 | 2.16e-09 | 0.22 | L1rho B 3.160e-04 W 3.160e-04 |
| edd | 0.3333 | 10 | 256 | 3.16e-08 | 2.00e-09 | 3.02e-09 | 7.17e-10 | 0.096 | L1rho B 7.808e-05 W 7.808e-05 |
| edd | 0.3333 | 100 | 16 | 3.16e-08 | 2.82e-08 | 2.11e-08 | 5.33e-09 | 0.89 | L1rho B 4.659e-02 W 4.659e-02 |
| edd | 0.3333 | 100 | 32 | 3.16e-08 | 3.60e-08 | 4.00e-08 | 3.05e-09 | 1.3 | L1rho B 1.148e-02 W 1.148e-02 |
| edd | 0.3333 | 100 | 64 | 3.16e-08 | 2.22e-09 | 3.41e-09 | 1.77e-10 | 0.11 | L1rho B 2.028e-03 W 2.028e-03 |
| edd | 0.3333 | 100 | 128 | 3.16e-08 | 6.00e-09 | 2.35e-08 | 4.42e-10 | 0.74 | L1rho B 4.143e-04 W 4.143e-04 |
| edd | 0.3333 | 100 | 256 | 3.16e-08 | 4.00e-09 | 1.11e-08 | 6.67e-10 | 0.35 | L1rho B 9.294e-05 W 9.294e-05 |
| edd | 0.3333 | 1000 | 16 | 3.16e-08 | 6.48e-08 | 5.53e-08 | 1.07e-08 | 2.1 | L1rho B 6.416e-02 W 6.416e-02 |
| edd | 0.3333 | 1000 | 32 | 3.16e-08 | 4.77e-08 | 5.04e-08 | 3.69e-09 | 1.6 | L1rho B 2.240e-02 W 2.240e-02 |
| edd | 0.3333 | 1000 | 64 | 3.16e-08 | 6.31e-08 | 5.68e-08 | 9.16e-10 | 2 | L1rho B 4.955e-03 W 4.955e-03 |
| edd | 0.3333 | 1000 | 128 | 3.16e-08 | 1.17e-08 | 1.16e-08 | 1.68e-10 | 0.37 | L1rho B 9.197e-04 W 9.197e-04 |
| edd | 0.3333 | 1000 | 256 | 3.16e-08 | 2.62e-08 | 3.32e-08 | 2.72e-10 | 1.1 | L1rho B 1.524e-04 W 1.524e-04 |
| edd | 3.333 | 0.01 | 16 | 3.16e-04 | 1.04e-07 | 7.11e-06 | 2.76e-08 | 0.022 | L1rho B 6.102e-02 W 6.102e-02 |
| edd | 3.333 | 0.01 | 32 | 3.16e-04 | 6.59e-08 | 1.01e-05 | 1.60e-08 | 0.032 | L1rho B 2.114e-02 W 2.114e-02 |
| edd | 3.333 | 0.01 | 64 | 3.16e-04 | 3.40e-08 | 6.25e-06 | 8.41e-09 | 0.02 | L1rho B 5.369e-03 W 5.369e-03 |
| edd | 3.333 | 0.01 | 128 | 3.16e-04 | 1.81e-08 | 4.59e-06 | 4.15e-09 | 0.015 | L1rho B 1.243e-03 W 1.243e-03 |
| edd | 3.333 | 0.01 | 256 | 3.16e-04 | 1.18e-08 | 4.95e-06 | 2.16e-09 | 0.016 | L1rho B 2.719e-04 W 2.719e-04 |
| edd | 3.333 | 0.1 | 16 | 3.16e-04 | 1.13e-06 | 2.28e-06 | 3.65e-07 | 0.0072 | L1rho B 4.656e-02 W 4.656e-02 |
| edd | 3.333 | 0.1 | 32 | 3.16e-04 | 7.50e-07 | 1.34e-05 | 2.22e-07 | 0.042 | L1rho B 1.078e-02 W 1.078e-02 |
| edd | 3.333 | 0.1 | 64 | 3.16e-04 | 3.98e-07 | 7.65e-06 | 1.17e-07 | 0.024 | L1rho B 2.662e-03 W 2.662e-03 |
| edd | 3.333 | 0.1 | 128 | 3.16e-04 | 2.03e-07 | 3.42e-06 | 5.92e-08 | 0.011 | L1rho B 6.691e-04 W 6.692e-04 |
| edd | 3.333 | 0.1 | 256 | 3.16e-04 | 1.03e-07 | 2.07e-06 | 2.97e-08 | 0.0065 | L1rho B 1.664e-04 W 1.665e-04 |
| edd | 3.333 | 1 | 16 | 3.16e-04 | 6.56e-06 | 2.18e-06 | 2.48e-06 | 0.021 | L1rho B 3.615e-02 W 3.615e-02 |
| edd | 3.333 | 1 | 32 | 3.16e-04 | 4.90e-06 | 1.95e-06 | 1.84e-06 | 0.016 | L1rho B 1.155e-02 W 1.155e-02 |
| edd | 3.333 | 1 | 64 | 3.16e-04 | 2.78e-06 | 1.28e-06 | 1.07e-06 | 0.0088 | L1rho B 2.640e-03 W 2.641e-03 |
| edd | 3.333 | 1 | 128 | 3.16e-04 | 1.48e-06 | 6.24e-07 | 5.72e-07 | 0.0047 | L1rho B 5.774e-04 W 5.782e-04 |
| edd | 3.333 | 1 | 256 | 3.16e-04 | 7.64e-07 | 3.12e-07 | 2.94e-07 | 0.0024 | L1rho B 1.333e-04 W 1.337e-04 |
| edd | 3.333 | 10 | 16 | 3.16e-04 | 1.29e-06 | 6.62e-07 | 1.77e-06 | 0.0041 | L1rho B 9.007e-03 W 9.008e-03 |
| edd | 3.333 | 10 | 32 | 3.16e-04 | 2.78e-06 | 2.00e-06 | 4.13e-06 | 0.0088 | L1rho B 2.262e-03 W 2.264e-03 |
| edd | 3.333 | 10 | 64 | 3.16e-04 | 2.38e-06 | 1.91e-06 | 3.66e-06 | 0.0075 | L1rho B 5.904e-04 W 5.919e-04 |
| edd | 3.333 | 10 | 128 | 3.16e-04 | 1.55e-06 | 1.29e-06 | 2.40e-06 | 0.0049 | L1rho B 1.516e-04 W 1.526e-04 |
| edd | 3.333 | 10 | 256 | 3.16e-04 | 8.72e-07 | 7.40e-07 | 1.36e-06 | 0.0028 | L1rho B 3.893e-05 W 3.948e-05 |
| edd | 3.333 | 100 | 16 | 3.16e-04 | 2.52e-06 | 1.97e-06 | 4.86e-07 | 0.008 | L1rho B 5.903e-02 W 5.903e-02 |
| edd | 3.333 | 100 | 32 | 3.16e-04 | 4.62e-07 | 4.31e-07 | 1.22e-07 | 0.0015 | L1rho B 1.461e-02 W 1.461e-02 |
| edd | 3.333 | 100 | 64 | 3.16e-04 | 1.31e-08 | 2.44e-08 | 1.62e-10 | 7.7e-05 | L1rho B 3.602e-03 W 3.602e-03 |
| edd | 3.333 | 100 | 128 | 3.16e-04 | 4.13e-08 | 4.78e-08 | 1.12e-08 | 0.00015 | L1rho B 8.974e-04 W 8.974e-04 |
| edd | 3.333 | 100 | 256 | 3.16e-04 | 8.93e-08 | 8.56e-08 | 2.58e-08 | 0.00028 | L1rho B 2.243e-04 W 2.244e-04 |
| edd | 3.333 | 1000 | 16 | 3.16e-04 | 1.28e-05 | 9.52e-06 | 1.79e-06 | 0.04 | L1rho B 6.949e-02 W 6.949e-02 |
| edd | 3.333 | 1000 | 32 | 3.16e-04 | 8.34e-06 | 7.58e-06 | 1.26e-06 | 0.026 | L1rho B 1.802e-02 W 1.801e-02 |
| edd | 3.333 | 1000 | 64 | 3.16e-04 | 1.42e-06 | 1.40e-06 | 3.29e-07 | 0.0045 | L1rho B 4.578e-03 W 4.578e-03 |
| edd | 3.333 | 1000 | 128 | 3.16e-04 | 2.05e-07 | 1.91e-07 | 4.75e-08 | 0.00065 | L1rho B 1.131e-03 W 1.130e-03 |
| edd | 3.333 | 1000 | 256 | 3.16e-04 | 1.07e-07 | 8.90e-08 | 4.59e-09 | 0.00034 | L1rho B 2.719e-04 W 2.719e-04 |

per (closure, P): max over tau, N of (dev/g^4, dev)
('edd', 0.03333333333333333) dev/g^4 1.24e+06  dev 3.903e-06  g^4 3.160e-12
('edd', 0.3333333333333333) dev/g^4 154  dev 4.878e-06  g^4 3.160e-08
('edd', 3.3333333333333335) dev/g^4 0.0425  dev 1.342e-05  g^4 3.160e-04
```

### analyze.py via view_p18h_vetf (PROD = FIN, FINAL = q)
```
| P | tau | arm | p rho | p E raw | p E mean-rm | p E raw 64->128 | L1 rho 256 | L1 E 256 | L1 E-mrm 256 | mean off 256 (A|dE|) | done N |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.0333 | 0.01 | PROD | 2.00 | 1.91 | 1.94 | 1.89 | 6.939e-04 | 2.800e-04 | 2.745e-04 | 7.946e-05 | 16,32,64,128,256 |
| 0.0333 | 0.01 | FINAL | 2.00 | 1.91 | 1.94 | 1.89 | 6.939e-04 | 2.800e-04 | 2.747e-04 | 7.898e-05 | 16,32,64,128,256 |
| 0.0333 | 0.1 | PROD | 2.07 | 1.95 | 1.95 | 1.94 | 6.005e-04 | 2.917e-04 | 2.908e-04 | 3.455e-05 | 16,32,64,128,256 |
| 0.0333 | 0.1 | FINAL | 2.07 | 1.95 | 1.95 | 1.94 | 6.005e-04 | 2.917e-04 | 2.908e-04 | 3.454e-05 | 16,32,64,128,256 |
| 0.0333 | 1 | PROD | 2.15 | 2.01 | 2.01 | 2.01 | 3.424e-04 | 2.177e-04 | 2.177e-04 | 2.709e-06 | 16,32,64,128,256 |
| 0.0333 | 1 | FINAL | 2.15 | 2.01 | 2.01 | 2.01 | 3.424e-04 | 2.177e-04 | 2.177e-04 | 2.710e-06 | 16,32,64,128,256 |
| 0.0333 | 10 | PROD | 2.06 | 2.04 | 2.04 | 2.05 | 2.360e-04 | 1.870e-04 | 1.870e-04 | 6.300e-07 | 16,32,64,128,256 |
| 0.0333 | 10 | FINAL | 2.06 | 2.04 | 2.04 | 2.05 | 2.360e-04 | 1.870e-04 | 1.870e-04 | 6.302e-07 | 16,32,64,128,256 |
| 0.0333 | 100 | PROD | 2.26 | 2.26 | 2.26 | 2.29 | 2.870e-04 | 2.798e-04 | 2.798e-04 | 3.871e-07 | 16,32,64,128,256 |
| 0.0333 | 100 | FINAL | 2.26 | 2.26 | 2.26 | 2.29 | 2.870e-04 | 2.798e-04 | 2.798e-04 | 3.872e-07 | 16,32,64,128,256 |
| 0.0333 | 1000 | PROD | 2.21 | 2.21 | 2.21 | 2.06 | 4.888e-04 | 4.860e-04 | 4.860e-04 | 3.606e-07 | 16,32,64,128,256 |
| 0.0333 | 1000 | FINAL | 2.21 | 2.21 | 2.21 | 2.06 | 4.888e-04 | 4.860e-04 | 4.860e-04 | 3.604e-07 | 16,32,64,128,256 |
| 0.3333 | 0.01 | PROD | 2.07 | 1.91 | 1.93 | 1.89 | 5.956e-04 | 2.745e-04 | 2.706e-04 | 6.251e-05 | 16,32,64,128,256 |
| 0.3333 | 0.01 | FINAL | 2.07 | 1.91 | 1.93 | 1.89 | 5.956e-04 | 2.746e-04 | 2.708e-04 | 6.302e-05 | 16,32,64,128,256 |
| 0.3333 | 0.1 | PROD | 2.17 | 1.97 | 1.97 | 1.95 | 2.956e-04 | 2.462e-04 | 2.459e-04 | 1.668e-05 | 16,32,64,128,256 |
| 0.3333 | 0.1 | FINAL | 2.17 | 1.97 | 1.97 | 1.95 | 2.956e-04 | 2.462e-04 | 2.459e-04 | 1.668e-05 | 16,32,64,128,256 |
| 0.3333 | 1 | PROD | 2.01 | 2.00 | 2.00 | 1.97 | 1.507e-04 | 1.600e-04 | 1.600e-04 | 6.475e-07 | 16,32,64,128,256 |
| 0.3333 | 1 | FINAL | 2.01 | 2.00 | 2.00 | 1.97 | 1.507e-04 | 1.600e-04 | 1.600e-04 | 6.452e-07 | 16,32,64,128,256 |
| 0.3333 | 10 | PROD | 2.00 | 2.00 | 2.00 | 2.02 | 8.184e-05 | 5.844e-05 | 5.844e-05 | 2.773e-07 | 16,32,64,128,256 |
| 0.3333 | 10 | FINAL | 2.00 | 2.00 | 2.00 | 2.02 | 8.184e-05 | 5.844e-05 | 5.844e-05 | 2.776e-07 | 16,32,64,128,256 |
| 0.3333 | 100 | PROD | 2.14 | 2.14 | 2.14 | 2.28 | 9.210e-05 | 8.361e-05 | 8.361e-05 | 1.575e-07 | 16,32,64,128,256 |
| 0.3333 | 100 | FINAL | 2.14 | 2.14 | 2.14 | 2.28 | 9.210e-05 | 8.361e-05 | 8.361e-05 | 1.572e-07 | 16,32,64,128,256 |
| 0.3333 | 1000 | PROD | 2.61 | 2.63 | 2.63 | 2.46 | 1.484e-04 | 1.391e-04 | 1.391e-04 | 9.709e-08 | 16,32,64,128,256 |
| 0.3333 | 1000 | FINAL | 2.61 | 2.63 | 2.63 | 2.46 | 1.484e-04 | 1.391e-04 | 1.391e-04 | 9.696e-08 | 16,32,64,128,256 |
| 3.333 | 0.01 | PROD | 2.19 | 1.96 | 1.96 | 1.90 | 2.720e-04 | 3.577e-04 | 3.572e-04 | 3.039e-05 | 16,32,64,128,256 |
| 3.333 | 0.01 | FINAL | 2.19 | 1.96 | 1.97 | 1.90 | 2.720e-04 | 3.569e-04 | 3.564e-04 | 3.103e-05 | 16,32,64,128,256 |
| 3.333 | 0.1 | PROD | 2.01 | 1.98 | 1.98 | 1.91 | 1.656e-04 | 3.122e-04 | 3.122e-04 | 1.308e-06 | 16,32,64,128,256 |
| 3.333 | 0.1 | FINAL | 2.01 | 1.98 | 1.98 | 1.91 | 1.657e-04 | 3.116e-04 | 3.116e-04 | 1.412e-06 | 16,32,64,128,256 |
| 3.333 | 1 | PROD | 2.12 | 1.91 | 1.91 | 1.81 | 1.297e-04 | 1.258e-04 | 1.258e-04 | -1.130e-07 | 16,32,64,128,256 |
| 3.333 | 1 | FINAL | 2.12 | 1.91 | 1.91 | 1.81 | 1.301e-04 | 1.257e-04 | 1.257e-04 | -1.218e-07 | 16,32,64,128,256 |
| 3.333 | 10 | PROD | 2.00 | 2.01 | 2.01 | 2.00 | 4.638e-05 | 6.230e-05 | 6.230e-05 | 9.356e-08 | 16,32,64,128,256 |
| 3.333 | 10 | FINAL | 1.99 | 2.00 | 2.00 | 2.00 | 4.693e-05 | 6.273e-05 | 6.273e-05 | 9.343e-08 | 16,32,64,128,256 |
| 3.333 | 100 | PROD | 2.02 | 2.04 | 2.04 | 2.06 | 2.283e-04 | 2.383e-04 | 2.383e-04 | 3.755e-08 | 16,32,64,128,256 |
| 3.333 | 100 | FINAL | 2.02 | 2.04 | 2.04 | 2.06 | 2.283e-04 | 2.384e-04 | 2.384e-04 | 3.755e-08 | 16,32,64,128,256 |
| 3.333 | 1000 | PROD | 2.03 | 2.04 | 2.04 | 2.04 | 2.859e-04 | 2.965e-04 | 2.965e-04 | 1.247e-08 | 16,32,64,128,256 |
| 3.333 | 1000 | FINAL | 2.03 | 2.04 | 2.04 | 2.04 | 2.858e-04 | 2.965e-04 | 2.965e-04 | 1.312e-08 | 16,32,64,128,256 |

FINAL/PROD L1 at N256 (rho, E):
  P 0.0333 tau 0.01: 1.000 1.000
  P 0.0333 tau 0.1: 1.000 1.000
  P 0.0333 tau 1: 1.000 1.000
  P 0.0333 tau 10: 1.000 1.000
  P 0.0333 tau 100: 1.000 1.000
  P 0.0333 tau 1000: 1.000 1.000
  P 0.3333 tau 0.01: 1.000 1.000
  P 0.3333 tau 0.1: 1.000 1.000
  P 0.3333 tau 1: 1.000 1.000
  P 0.3333 tau 10: 1.000 1.000
  P 0.3333 tau 100: 1.000 1.000
  P 0.3333 tau 1000: 1.000 1.000
  P 3.333 tau 0.01: 1.000 0.998
  P 3.333 tau 0.1: 1.000 0.998
  P 3.333 tau 1: 1.003 1.000
  P 3.333 tau 10: 1.012 1.007
  P 3.333 tau 100: 1.000 1.000
  P 3.333 tau 1000: 1.000 1.000
```

### analyze_edd.py via view_p18h_edd (PROD = FIN, FINAL = q)
```
| P | tau | arm | p rho 16-32-64-128-256 | p E 16..256 | p E mean-rm 128->256 | L1 rho 256 | L1 E 256 | |dw|/|w| 256 vs edd | vs eddc | done N |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.0333 | 0.01 | PROD | 1.32 1.84 1.93 2.00 | 2.41 1.73 1.88 1.88 | 1.93 | 6.939e-04 | 2.769e-04 | 6.93e-05 | 7.06e-05 | 16,32,64,128,256 |
| 0.0333 | 0.01 | FINAL | 1.32 1.84 1.93 2.00 | 2.41 1.73 1.88 1.88 | 1.93 | 6.939e-04 | 2.771e-04 | 6.93e-05 | 7.06e-05 | 16,32,64,128,256 |
| 0.0333 | 0.1 | PROD | 1.34 1.85 2.00 2.07 | 2.21 1.85 1.91 1.93 | 1.94 | 5.995e-04 | 2.666e-04 | 7.08e-05 | 1.88e-04 | 16,32,64,128,256 |
| 0.0333 | 0.1 | FINAL | 1.34 1.85 2.00 2.07 | 2.21 1.85 1.91 1.93 | 1.94 | 5.995e-04 | 2.667e-04 | 7.08e-05 | 1.88e-04 | 16,32,64,128,256 |
| 0.0333 | 1 | PROD | 1.54 1.94 2.06 2.15 | 2.04 2.00 2.01 1.99 | 1.99 | 3.157e-04 | 2.049e-04 | 7.59e-05 | 1.68e-03 | 16,32,64,128,256 |
| 0.0333 | 1 | FINAL | 1.54 1.94 2.06 2.15 | 2.04 2.00 2.01 1.99 | 1.99 | 3.157e-04 | 2.049e-04 | 7.59e-05 | 1.68e-03 | 16,32,64,128,256 |
| 0.0333 | 10 | PROD | 1.97 1.98 2.03 2.05 | 2.09 2.05 2.03 2.02 | 2.02 | 2.224e-04 | 1.768e-04 | 7.00e-05 | 1.79e-03 | 16,32,64,128,256 |
| 0.0333 | 10 | FINAL | 1.97 1.98 2.03 2.05 | 2.09 2.05 2.03 2.02 | 2.02 | 2.224e-04 | 1.768e-04 | 7.00e-05 | 1.79e-03 | 16,32,64,128,256 |
| 0.0333 | 100 | PROD | 1.56 2.10 2.32 2.26 | 1.58 2.10 2.31 2.25 | 2.25 | 2.792e-04 | 2.735e-04 | 6.49e-05 | 6.27e-04 | 16,32,64,128,256 |
| 0.0333 | 100 | FINAL | 1.56 2.10 2.32 2.26 | 1.58 2.10 2.31 2.25 | 2.25 | 2.792e-04 | 2.735e-04 | 6.49e-05 | 6.27e-04 | 16,32,64,128,256 |
| 0.0333 | 1000 | PROD | 1.34 1.89 2.06 2.21 | 1.34 1.89 2.06 2.21 | 2.21 | 4.886e-04 | 4.856e-04 | 6.24e-05 | 6.04e-04 | 16,32,64,128,256 |
| 0.0333 | 1000 | FINAL | 1.34 1.89 2.06 2.21 | 1.34 1.89 2.06 2.21 | 2.21 | 4.886e-04 | 4.856e-04 | 6.24e-05 | 6.04e-04 | 16,32,64,128,256 |
| 0.3333 | 0.01 | PROD | 1.34 1.85 2.00 2.07 | 2.45 1.69 1.87 1.88 | 1.93 | 5.957e-04 | 2.696e-04 | 6.99e-05 | 1.89e-04 | 16,32,64,128,256 |
| 0.3333 | 0.01 | FINAL | 1.34 1.85 2.00 2.07 | 2.45 1.69 1.87 1.88 | 1.93 | 5.957e-04 | 2.695e-04 | 6.99e-05 | 1.89e-04 | 16,32,64,128,256 |
| 0.3333 | 0.1 | PROD | 1.54 1.96 2.08 2.17 | 1.69 1.66 1.90 1.96 | 1.96 | 2.932e-04 | 2.587e-04 | 7.57e-05 | 1.80e-03 | 16,32,64,128,256 |
| 0.3333 | 0.1 | FINAL | 1.54 1.96 2.08 2.17 | 1.69 1.66 1.90 1.96 | 1.96 | 2.932e-04 | 2.587e-04 | 7.57e-05 | 1.80e-03 | 16,32,64,128,256 |
| 0.3333 | 1 | PROD | 2.08 1.98 1.98 2.01 | 2.06 1.88 1.97 1.99 | 1.99 | 1.649e-04 | 1.811e-04 | 8.43e-05 | 5.69e-03 | 16,32,64,128,256 |
| 0.3333 | 1 | FINAL | 2.08 1.98 1.98 2.01 | 2.06 1.88 1.97 1.99 | 1.99 | 1.649e-04 | 1.811e-04 | 8.43e-05 | 5.69e-03 | 16,32,64,128,256 |
| 0.3333 | 10 | PROD | 1.98 2.11 2.14 2.02 | 2.07 2.07 2.02 2.00 | 2.00 | 7.808e-05 | 5.545e-05 | 7.14e-05 | 9.68e-03 | 16,32,64,128,256 |
| 0.3333 | 10 | FINAL | 1.98 2.11 2.14 2.02 | 2.07 2.07 2.02 2.00 | 2.00 | 7.808e-05 | 5.545e-05 | 7.14e-05 | 9.68e-03 | 16,32,64,128,256 |
| 0.3333 | 100 | PROD | 2.02 2.50 2.29 2.16 | 2.09 2.50 2.28 2.15 | 2.15 | 9.294e-05 | 8.497e-05 | 2.04e-05 | 6.67e-03 | 16,32,64,128,256 |
| 0.3333 | 100 | FINAL | 2.02 2.50 2.29 2.16 | 2.09 2.50 2.28 2.15 | 2.15 | 9.294e-05 | 8.497e-05 | 2.04e-05 | 6.67e-03 | 16,32,64,128,256 |
| 0.3333 | 1000 | PROD | 1.52 2.18 2.43 2.59 | 1.50 2.15 2.44 2.61 | 2.61 | 1.524e-04 | 1.433e-04 | 1.78e-05 | 6.63e-03 | 16,32,64,128,256 |
| 0.3333 | 1000 | FINAL | 1.52 2.18 2.43 2.59 | 1.50 2.15 2.44 2.61 | 2.61 | 1.524e-04 | 1.433e-04 | 1.78e-05 | 6.63e-03 | 16,32,64,128,256 |
| 3.3333 | 0.01 | PROD | 1.53 1.98 2.11 2.19 | 1.27 1.44 1.89 1.96 | 1.97 | 2.719e-04 | 3.522e-04 | 6.57e-05 | 1.80e-03 | 16,32,64,128,256 |
| 3.3333 | 0.01 | FINAL | 1.53 1.98 2.11 2.19 | 1.27 1.44 1.89 1.96 | 1.97 | 2.719e-04 | 3.516e-04 | 6.57e-05 | 1.80e-03 | 16,32,64,128,256 |
| 3.3333 | 0.1 | PROD | 2.11 2.02 1.99 2.01 | 0.86 1.70 1.93 1.99 | 1.99 | 1.664e-04 | 3.977e-04 | 7.40e-05 | 5.49e-03 | 16,32,64,128,256 |
| 3.3333 | 0.1 | FINAL | 2.11 2.02 1.99 2.01 | 0.86 1.70 1.93 1.99 | 1.99 | 1.665e-04 | 3.967e-04 | 7.40e-05 | 5.49e-03 | 16,32,64,128,256 |
| 3.3333 | 1 | PROD | 1.65 2.13 2.19 2.11 | 2.23 1.90 1.96 2.00 | 2.00 | 1.333e-04 | 1.225e-04 | 7.46e-05 | 6.82e-03 | 16,32,64,128,256 |
| 3.3333 | 1 | FINAL | 1.65 2.13 2.19 2.11 | 2.23 1.90 1.96 2.00 | 2.00 | 1.337e-04 | 1.224e-04 | 7.48e-05 | 6.82e-03 | 16,32,64,128,256 |
| 3.3333 | 10 | PROD | 1.99 1.94 1.96 1.96 | 1.89 1.94 1.97 1.98 | 1.98 | 3.893e-05 | 5.607e-05 | 9.41e-05 | 8.38e-02 | 16,32,64,128,256 |
| 3.3333 | 10 | FINAL | 1.99 1.94 1.96 1.95 | 1.89 1.93 1.97 1.98 | 1.98 | 3.948e-05 | 5.647e-05 | 9.52e-05 | 8.38e-02 | 16,32,64,128,256 |
| 3.3333 | 100 | PROD | 2.01 2.02 2.00 2.00 | 2.03 2.05 2.04 2.02 | 2.02 | 2.243e-04 | 2.344e-04 | 7.16e-05 | 6.48e-02 | 16,32,64,128,256 |
| 3.3333 | 100 | FINAL | 2.01 2.02 2.00 2.00 | 2.03 2.05 2.04 2.02 | 2.02 | 2.244e-04 | 2.345e-04 | 7.16e-05 | 6.48e-02 | 16,32,64,128,256 |
| 3.3333 | 1000 | PROD | 1.95 1.98 2.02 2.06 | 1.95 2.03 2.06 2.07 | 2.07 | 2.719e-04 | 2.821e-04 | 7.04e-05 | 6.47e-02 | 16,32,64,128,256 |
| 3.3333 | 1000 | FINAL | 1.95 1.98 2.02 2.06 | 1.95 2.03 2.06 2.07 | 2.07 | 2.719e-04 | 2.821e-04 | 7.04e-05 | 6.47e-02 | 16,32,64,128,256 |
```

## (c) Raven IDORT set, q = prod-1118 + EXTRA_HG (one_q_orion.sh): **correctness PASS on 160/160 rows; q != d6int on 126 rows**
main 41 / rsx 12 / spxx 100 / OMP 7 (OMP via sched.py with NC 10 and the build_omp binary, W=omp_q, as the TASK). Nodes orion052-055 from
12:30 PST 11-13; walls: spxx 4671 s, OMP 12847 s (rows run one after another), main 21778 s, rsx 30783 s. Smoke: list_smoke 5/5 q rc=0 (plain EXTRA,
run before the UPDATE; the full set uses EXTRA_HG). The HG keys are in every q deck (inserted by splitkeys; checked in runs/*/q/in.athinput).

- **Gates: all 160 q rows rc=0, FATAL 0, NON-CONVERGED 0, DT_RETRY 0.** cmp_o.sh with q as the D6 arm (view dir p -> q) vs FIN:
  **0 FLAG_RC, 0 FLAG_NC, 0 FLAG_PICARD, 0 MISSING.**
- **q vs d6int (p), cmpout per row: 34 BITWISE, 126 DIFF.** This is NOT the expected "q = d6int bitwise". Where it comes from (q logs):
  - **vimp_lag_local fires in 0 rows; dtr_poslim in 0; adm_poslim in 2 rows.** So the HG positivity keys are almost inert on this set.
  - The Zalesak cap can act in 29 DIFF rows (q max rounds > 4: 5, 9, up to 75 uncapped in p-type rows).
  - **97 DIFF rows have no cap action (rounds <= 4), no adm / lag / dtr firing:** there the only remaining change from the D6 deck is
    **time2_stage_local_be = false** (decks_d6 have true; EXTRA turns it off). It is the dominant source of q != d6int, e.g. all spx rsw / rlat
    families, nltek1-8 (Picard mean p 3.03-3.12 -> q 2.26-2.85, max 4-5 -> 3), cb/cbg, eq2, pipe. face_mmat limits faces in 41 q rows
    (counts differ from p where stage_local_be changes the stage path).
  - cmpout's maxrel saturates at 2 (sign-changing near-zero fields), so the size of each change is better read from ndiff, from the Picard change
    in the table below, and from your analyzers. I have no IDORT row analyzers on Orion: the accuracy rule (metric <= 1.05 x FIN) is **not
    evaluated** for (c). The q outputs sit next to p and f in /orion/ptmp/jinma/idortd6int_1113/runs/<row>/q and omp_q/runs/<row>/q.
- **Correction:** I first parsed these log counters by position and briefly had face_mmat-faces mislabelled as lag-local (41 rows). The numbers
  above use label-based extraction. That error was never posted before this section.

### q vs p per row (q_vs_p_rows.md)
| row | q-vs-p | ndiff | maxrel | Picard mean p -> q | face_mmat limited faces q (p) | max Zalesak rounds q | lag_local | adm_poslim | dtr_poslim | candidate key(s) |
|---|---|---|---|---|---|---|---|---|---|---|
| omp_q/runs/cart/ddf/FIN3_d | BITWISE | - | - | 3.981 -> 3.981 | 0 (0) | 0 | 0 | 0 | 0 | - |
| omp_q/runs/cart/ddn/FIN3_d | BITWISE | - | - | 3.368 -> 3.368 | 0 (0) | 0 | 0 | 0 | 0 | - |
| omp_q/runs/cart/rs_m1.2/FIN3_n512 | BITWISE | - | - | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | - |
| omp_q/runs/cart/rs_m2/FIN3_n1024 | DIFF | 12 | 2.000e+00 | 3.000 -> 3.000 | 94596180 (1211532) | 5 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| omp_q/runs/cart/rs_m2/FIN3_n512 | DIFF | 12 | 1.999e+00 | 3.000 -> 3.000 | 5866066 (269265) | 5 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| omp_q/runs/cart/rs_m3/FIN3_n1024 | DIFF | 12 | 1.999e+00 | 3.000 -> 3.000 | 122903700 (709752) | 5 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| omp_q/runs/cart/rs_m3/FIN3_n512 | DIFF | 12 | 2.000e+00 | 3.457 -> 3.456 | 197483 (190256) | 5 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| runs/cart/cb/FIN3_r12 | DIFF | 1 | 1.169e-02 | 2.513 -> 2.615 | 288 (288) | 69 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| runs/cart/cb/FIN3_r4c0.5 | BITWISE | - | - | 2.267 -> 2.267 | 0 (0) | 0 | 0 | 0 | 0 | - |
| runs/cart/cb/FIN3_r4c32 | BITWISE | - | - | 2.702 -> 2.702 | 256 (256) | 61 | 0 | 0 | 0 | - |
| runs/cart/cb/FIN3_r4 | DIFF | 1 | 1.541e+00 | 2.800 -> 2.900 | 324 (303) | 66 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| runs/cart/cbg/FIN3_g8x32 | DIFF | 2 | 1.999e+00 | 2.949 -> 3.000 | 94623 (32004) | 75 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| runs/cart/ddf/FIN3_d | BITWISE | - | - | 3.981 -> 3.981 | 0 (0) | 0 | 0 | 0 | 0 | - |
| runs/cart/ddiff/FIN3_d | BITWISE | - | - | 3.981 -> 3.981 | 0 (0) | 0 | 0 | 0 | 0 | - |
| runs/cart/ddn/FIN3_d | BITWISE | - | - | 3.368 -> 3.368 | 0 (0) | 0 | 0 | 0 | 0 | - |
| runs/cart/eq1/FIN3_fine | BITWISE | - | - | 3.005 -> 3.005 | 0 (0) | 0 | 0 | 0 | 0 | - |
| runs/cart/eq2/FIN3_fine | DIFF | 80 | 2.000e+00 | 3.015 -> 3.003 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/cart/eq3t/FIN3_t | BITWISE | - | - | 6.185 -> 6.185 | 0 (0) | 0 | 0 | 0 | 0 | - |
| runs/cart/hcyl/FIN3_g8x32 | BITWISE | - | - | 2.173 -> 2.173 | 0 (0) | 0 | 0 | 0 | 0 | - |
| runs/cart/marsh/FIN3_c10 | BITWISE | - | - | 3.011 -> 3.011 | 21229 (21229) | 1 | 0 | 0 | 0 | - |
| runs/cart/marsh/FIN3_c1 | BITWISE | - | - | 3.000 -> 3.000 | 13 (13) | 1 | 0 | 0 | 0 | - |
| runs/cart/nltek1/FIN3_d | DIFF | 1 | 1.989e+00 | 3.050 -> 2.513 | 279502 (104283) | 5 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| runs/cart/nltek2/FIN3_d | DIFF | 1 | 1.999e+00 | 3.050 -> 2.613 | 312094 (110702) | 5 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| runs/cart/nltek3/FIN3_d | DIFF | 1 | 1.996e+00 | 3.118 -> 2.731 | 356207 (171006) | 5 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| runs/cart/nltek4/FIN3_d | DIFF | 1 | 1.998e+00 | 3.025 -> 2.849 | 434077 (179368) | 5 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| runs/cart/nltek5/FIN3_d | DIFF | 1 | 2.000e+00 | 3.042 -> 2.639 | 470097 (292224) | 5 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| runs/cart/nltek6/FIN3_d | DIFF | 1 | 1.996e+00 | 3.084 -> 2.454 | 501083 (336883) | 5 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| runs/cart/nltek7/FIN3_d | DIFF | 1 | 1.998e+00 | 3.034 -> 2.261 | 522726 (248304) | 5 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| runs/cart/nltek8/FIN3_d | DIFF | 1 | 2.000e+00 | 3.059 -> 2.412 | 534578 (304523) | 5 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| runs/cart/oneo/FIN3_b3k3 | BITWISE | - | - | 2.759 -> 2.759 | 35 (35) | 9 | 0 | 0 | 0 | - |
| runs/cart/oneo/FIN3_b3k5 | BITWISE | - | - | 2.823 -> 2.823 | 1404 (1404) | 9 | 0 | 0 | 0 | - |
| runs/cart/oneo/FIN3_b3k7 | BITWISE | - | - | 2.873 -> 2.873 | 5779 (5779) | 9 | 0 | 0 | 0 | - |
| runs/cart/pipe/FIN3_c2200 | DIFF | 296 | 2.000e+00 | 3.608 -> 3.545 | 40694820 (1642688) | 4096 | 0 | 6 | 0 | Zalesak cap (rounds>4) + time2_adm_poslim |
| runs/cart/pipe/FIN3_c220 | DIFF | 300 | 2.000e+00 | 2.993 -> 3.005 | 186683800 (1008061) | 1083 | 0 | 1 | 0 | Zalesak cap (rounds>4) + time2_adm_poslim |
| runs/cart/rs_m1.2/FIN3_n512 | BITWISE | - | - | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | - |
| runs/cart/rs_m2/FIN3_n1024 | DIFF | 12 | 2.000e+00 | 3.000 -> 3.000 | 94596180 (1211532) | 5 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| runs/cart/rs_m2/FIN3_n512 | DIFF | 12 | 1.999e+00 | 3.000 -> 3.000 | 5866066 (269265) | 5 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| runs/cart/rs_m3/FIN3_n1024 | DIFF | 12 | 1.999e+00 | 3.000 -> 3.000 | 122903700 (709752) | 5 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| runs/cart/rs_m3/FIN3_n512 | DIFF | 12 | 2.000e+00 | 3.457 -> 3.456 | 197483 (190256) | 5 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| runs/cart/rsx_m1.05/P_n512_c0.4 | BITWISE | - | - | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | - |
| runs/cart/rsx_m1.2/P_n512_c0.4 | BITWISE | - | - | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | - |
| runs/cart/rsx_m27/P_n1024_c0.4 | BITWISE | - | - | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | - |
| runs/cart/rsx_m2/P_n512_c0.4 | DIFF | 6 | 2.000e+00 | 3.000 -> 3.000 | 663446 (136448) | 5 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| runs/cart/rsx_m3/P_n512_c0.4 | DIFF | 6 | 1.997e+00 | 3.430 -> 3.429 | 198243 (190820) | 5 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| runs/cart/rsx_m50/P_n1024_c0.4 | BITWISE | - | - | 4.000 -> 4.000 | 0 (0) | 0 | 0 | 0 | 0 | - |
| runs/cart/rsx_m5/P_n1024_c0.4 | BITWISE | - | - | 3.061 -> 3.061 | 0 (0) | 0 | 0 | 0 | 0 | - |
| runs/cart/rsx_m5/P_n2048_c0.2 | BITWISE | - | - | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | - |
| runs/cart/rsx_m5/P_n2048_c0.4 | BITWISE | - | - | 3.995 -> 3.995 | 0 (0) | 0 | 0 | 0 | 0 | - |
| runs/cart/sdifk/FIN3_d | BITWISE | - | - | 2.000 -> 2.000 | 0 (0) | 1 | 0 | 0 | 0 | - |
| runs/cart/shd/FIN3_r12g2 | DIFF | 2 | 2.000e+00 | 4.105 -> 4.053 | 105147 (85934) | 171 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| runs/cart/shd/FIN3_r36g2 | DIFF | 2 | 2.000e+00 | 4.000 -> 4.000 | 98325 (82309) | 163 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| runs/cart/shdg/FIN3_g8x32 | DIFF | 2 | 2.000e+00 | 4.000 -> 4.000 | 103487 (80181) | 191 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| runs/cart/slabc/FIN3_n128 | BITWISE | - | - | 2.102 -> 2.102 | 0 (0) | 0 | 0 | 0 | 0 | - |
| runs/cart/slabc/FIN3_n256 | BITWISE | - | - | 2.094 -> 2.094 | 0 (0) | 0 | 0 | 0 | 0 | - |
| runs/cart/slabc/FIN3_n32 | BITWISE | - | - | 2.039 -> 2.039 | 0 (0) | 0 | 0 | 0 | 0 | - |
| runs/cart/slabc/FIN3_n512 | BITWISE | - | - | 2.118 -> 2.118 | 0 (0) | 0 | 0 | 0 | 0 | - |
| runs/cart/slabc/FIN3_n64 | BITWISE | - | - | 2.063 -> 2.063 | 0 (0) | 0 | 0 | 0 | 0 | - |
| runs/sp/rsx_m27/P_n1024_c0.4 | DIFF | 6 | 1.943e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/sp/rsx_m2/P_n512_c0.4 | DIFF | 6 | 1.999e+00 | 3.000 -> 3.000 | 1279831 (1279869) | 3 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/sp/rsx_m5/P_n1024_c0.4 | DIFF | 6 | 2.000e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/marsh_r1_u/L128 | DIFF | 12 | 1.978e+00 | 3.000 -> 3.000 | 84 (84) | 1 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/marsh_r1_u/L256 | DIFF | 12 | 1.989e+00 | 3.000 -> 3.000 | 350 (350) | 1 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/marsh_r1_u/L32 | DIFF | 12 | 1.786e+00 | 3.000 -> 3.000 | 29 (29) | 1 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/marsh_r1_u/L512 | DIFF | 12 | 1.993e+00 | 3.000 -> 3.000 | 271 (271) | 1 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/marsh_r1_u/L64 | DIFF | 12 | 1.990e+00 | 3.000 -> 3.000 | 22 (22) | 1 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/patmT_L16_s/L128 | DIFF | 3 | 1.991e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/patmT_L16_s/L256 | DIFF | 3 | 1.981e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/patmT_L16_s/L32 | DIFF | 3 | 1.974e+00 | 3.000 -> 3.000 | 3155 (3348) | 5 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| runs/spx/patmT_L16_s/L512 | DIFF | 3 | 1.990e+00 | 3.013 -> 3.006 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/patmT_L16_s/L64 | DIFF | 3 | 1.989e+00 | 3.000 -> 3.000 | 6926 (6902) | 5 | 0 | 0 | 0 | Zalesak cap (rounds>4) |
| runs/spx/rlat/L128 | DIFF | 19 | 1.258e-07 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rlat/L16 | DIFF | 7 | 1.032e-07 | 3.098 -> 3.098 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rlat/L256 | DIFF | 21 | 1.142e-07 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rlat/L32 | DIFF | 7 | 3.279e-07 | 3.029 -> 3.029 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rlat/L64 | DIFF | 14 | 1.183e-07 | 3.006 -> 3.006 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rlat_om5/L128 | DIFF | 24 | 8.293e-03 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rlat_om5/L16 | DIFF | 24 | 4.731e-05 | 3.106 -> 3.106 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rlat_om5/L256 | DIFF | 24 | 1.531e-03 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rlat_om5/L32 | DIFF | 24 | 4.451e-04 | 3.029 -> 3.029 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rlat_om5/L64 | DIFF | 24 | 3.312e-04 | 3.008 -> 3.008 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rlatt/L128 | DIFF | 22 | 2.166e-07 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rlatt/L16 | DIFF | 18 | 1.054e-06 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rlatt/L256 | DIFF | 16 | 4.413e-07 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rlatt/L32 | DIFF | 19 | 1.159e-07 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rlatt/L64 | DIFF | 12 | 1.104e-07 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rlatt_om5/L128 | DIFF | 24 | 1.016e-01 | 3.003 -> 3.003 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rlatt_om5/L16 | DIFF | 24 | 1.182e-02 | 3.049 -> 3.049 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rlatt_om5/L256 | DIFF | 24 | 9.002e-01 | 3.001 -> 3.001 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rlatt_om5/L32 | DIFF | 24 | 1.706e-02 | 3.021 -> 3.021 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rlatt_om5/L64 | DIFF | 24 | 2.792e-01 | 3.010 -> 3.010 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rsw_C/L128 | DIFF | 3 | 1.372e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rsw_C/L256 | DIFF | 3 | 1.550e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rsw_C/L32 | DIFF | 3 | 3.060e-02 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rsw_C/L512 | DIFF | 3 | 1.985e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rsw_C/L64 | DIFF | 3 | 9.610e-02 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rsw/L128 | DIFF | 3 | 1.914e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rsw/L256 | DIFF | 3 | 2.000e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rsw/L32 | DIFF | 3 | 9.200e-01 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rsw/L512 | DIFF | 3 | 2.000e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rsw/L64 | DIFF | 3 | 9.610e-02 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswt_C/L128 | DIFF | 3 | 1.990e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswt_C/L256 | DIFF | 3 | 1.996e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswt_C/L32 | DIFF | 3 | 7.156e-01 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswt_C/L512 | DIFF | 3 | 1.999e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswt_C/L64 | DIFF | 3 | 1.999e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswt_dt/L16 | DIFF | 3 | 1.999e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswt_dt/L1 | DIFF | 3 | 1.998e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswt_dt/L2 | DIFF | 3 | 1.999e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswt_dt/L4 | DIFF | 3 | 1.996e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswt_dt/L8 | DIFF | 3 | 1.998e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswt/L128 | DIFF | 3 | 2.000e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswt/L256 | DIFF | 3 | 1.998e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswt/L32 | DIFF | 3 | 1.926e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswt/L512 | DIFF | 3 | 1.999e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswt/L64 | DIFF | 3 | 1.999e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswtv_C/L128 | DIFF | 3 | 1.992e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswtv_C/L256 | DIFF | 3 | 1.995e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswtv_C/L32 | DIFF | 3 | 1.870e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswtv_C/L512 | DIFF | 3 | 1.993e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswtv_C/L64 | DIFF | 3 | 1.972e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswtv_F_C/L128 | DIFF | 3 | 1.986e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswtv_F_C/L256 | DIFF | 3 | 1.992e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswtv_F_C/L32 | DIFF | 3 | 1.837e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswtv_F_C/L512 | DIFF | 3 | 1.997e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswtv_F_C/L64 | DIFF | 3 | 1.995e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswtv_F/L128 | DIFF | 3 | 1.988e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswtv_F/L256 | DIFF | 3 | 2.000e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswtv_F/L32 | DIFF | 3 | 1.966e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswtv_F/L512 | DIFF | 3 | 1.995e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswtv_F/L64 | DIFF | 3 | 1.995e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswtv/L128 | DIFF | 3 | 1.955e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswtv/L256 | DIFF | 3 | 1.997e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswtv/L32 | DIFF | 3 | 1.988e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswtv/L512 | DIFF | 3 | 1.999e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswtv/L64 | DIFF | 3 | 1.972e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswv_C/L128 | DIFF | 3 | 1.921e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswv_C/L256 | DIFF | 3 | 1.995e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswv_C/L32 | DIFF | 3 | 5.063e-01 | 3.016 -> 3.016 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswv_C/L512 | DIFF | 3 | 1.983e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswv_C/L64 | DIFF | 3 | 1.831e+00 | 3.004 -> 3.004 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswv_F_C/L128 | DIFF | 3 | 6.076e-01 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswv_F_C/L256 | DIFF | 3 | 1.972e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswv_F_C/L32 | DIFF | 3 | 9.261e-02 | 3.016 -> 3.016 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswv_F_C/L512 | DIFF | 3 | 1.998e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswv_F_C/L64 | DIFF | 3 | 1.826e+00 | 3.004 -> 3.004 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswv_F/L128 | DIFF | 3 | 1.954e+00 | 3.004 -> 3.004 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswv_F/L256 | DIFF | 3 | 1.997e+00 | 3.004 -> 3.004 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswv_F/L32 | DIFF | 3 | 1.079e-02 | 3.004 -> 3.004 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswv_F/L512 | DIFF | 3 | 1.995e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswv_F/L64 | DIFF | 3 | 1.826e+00 | 3.004 -> 3.004 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswv/L128 | DIFF | 3 | 1.651e+00 | 3.004 -> 3.004 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswv/L256 | DIFF | 3 | 1.995e+00 | 3.004 -> 3.004 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswv/L32 | DIFF | 3 | 3.504e-01 | 3.004 -> 3.004 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswv/L512 | DIFF | 3 | 1.996e+00 | 3.000 -> 3.000 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/spx/rswv/L64 | DIFF | 3 | 1.831e+00 | 3.004 -> 3.004 | 0 (0) | 0 | 0 | 0 | 0 | time2_stage_local_be=false (only candidate) |
| runs/xtra/g5_k12800 | BITWISE | - | - | 2.000 -> 2.000 | 0 (0) | 0 | 0 | 0 | 0 | - |
| runs/xtra/g5_k1280 | BITWISE | - | - | 2.000 -> 2.000 | 0 (0) | 0 | 0 | 0 | 0 | - |
| runs/xtra/g5_k12.8 | BITWISE | - | - | 2.000 -> 2.000 | 0 (0) | 0 | 0 | 0 | 0 | - |
| runs/xtra/g5_k128 | BITWISE | - | - | 2.000 -> 2.000 | 0 (0) | 0 | 0 | 0 | 0 | - |
| runs/xtra/shd432 | DIFF | 1 | 2.000e+00 | 2.089 -> 2.250 | 199656 (146642) | 6 | 0 | 0 | 0 | Zalesak cap (rounds>4) |

SUMMARY (160 rows): BITWISE 34, DIFF 126. Counters over all q rows: face_mmat limited faces > 0 in 41; vimp_lag_local > 0 in 0; adm_poslim > 0 in 2; dtr_poslim > 0 in 0.
DIFF rows with Zalesak rounds > 4 (cap can act): 29; with adm_poslim firing: 2; with lag_local firing: 0; DIFF rows where none of these act (only time2_stage_local_be=false is left): 97

### cmp_o.sh with q as D6 arm vs FIN (view_qf; whole)
```
== (1) rows: D6 (p) | FIN (f)
runs/cart/cb/FIN3_r12 | D6 rc=0 fatal=0 dtr=0 mean=2.615385e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=65 | FIN rc=0 fatal=0 dtr=0 mean=2.512821e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=64 |
runs/cart/cb/FIN3_r4 | D6 rc=0 fatal=0 dtr=0 mean=2.900000e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 t=64 | FIN rc=0 fatal=0 dtr=0 mean=2.800000e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 t=70 |
runs/cart/cb/FIN3_r4c0.5 | D6 rc=0 fatal=0 dtr=0 mean=2.266940e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=270 | FIN rc=0 fatal=0 dtr=0 mean=2.266940e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=281 |
runs/cart/cb/FIN3_r4c32 | D6 rc=0 fatal=0 dtr=0 mean=2.702128e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=29 | FIN rc=0 fatal=0 dtr=0 mean=2.702128e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=30 |
runs/cart/cbg/FIN3_g8x32 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=42 | FIN rc=0 fatal=0 dtr=0 mean=2.948718e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=45 |
runs/cart/ddf/FIN3_d | D6 rc=0 fatal=0 dtr=0 mean=3.981228e+00 max=1.000000e+01 NON-CONVERGED=0.000000e+00 t=4207 | FIN rc=0 fatal=0 dtr=0 mean=4.000000e+00 max=1.000000e+01 NON-CONVERGED=0.000000e+00 t=4271 |
runs/cart/ddiff/FIN3_d | D6 rc=0 fatal=0 dtr=0 mean=3.981228e+00 max=1.000000e+01 NON-CONVERGED=0.000000e+00 t=4210 | FIN rc=0 fatal=0 dtr=0 mean=4.000000e+00 max=1.000000e+01 NON-CONVERGED=0.000000e+00 t=4269 |
runs/cart/ddn/FIN3_d | D6 rc=0 fatal=0 dtr=0 mean=3.368488e+00 max=9.000000e+00 NON-CONVERGED=0.000000e+00 t=3812 | FIN rc=0 fatal=0 dtr=0 mean=3.368488e+00 max=9.000000e+00 NON-CONVERGED=0.000000e+00 t=3877 |
runs/cart/eq1/FIN3_fine | D6 rc=0 fatal=0 dtr=0 mean=3.005013e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=27 | FIN rc=0 fatal=0 dtr=0 mean=3.005013e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=28 |
runs/cart/eq2/FIN3_fine | D6 rc=0 fatal=0 dtr=0 mean=3.002506e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=29 | FIN rc=0 fatal=0 dtr=0 mean=3.015038e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=34 |
runs/cart/eq3t/FIN3_t | D6 rc=0 fatal=0 dtr=0 mean=6.185022e+00 max=9.000000e+00 NON-CONVERGED=0.000000e+00 t=53 | FIN rc=0 fatal=0 dtr=0 mean=6.198238e+00 max=1.000000e+01 NON-CONVERGED=0.000000e+00 t=57 |
runs/cart/hcyl/FIN3_g8x32 | D6 rc=0 fatal=0 dtr=0 mean=2.173228e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=278 | FIN rc=0 fatal=0 dtr=0 mean=2.173228e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=290 |
runs/cart/marsh/FIN3_c1 | D6 rc=0 fatal=0 dtr=0 mean=3.000097e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=995 | FIN rc=0 fatal=0 dtr=0 mean=3.000097e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=1021 |
runs/cart/marsh/FIN3_c10 | D6 rc=0 fatal=0 dtr=0 mean=3.010636e+00 max=7.000000e+00 NON-CONVERGED=0.000000e+00 t=105 | FIN rc=0 fatal=0 dtr=0 mean=3.010636e+00 max=7.000000e+00 NON-CONVERGED=0.000000e+00 t=111 |
runs/cart/nltek1/FIN3_d | D6 rc=0 fatal=0 dtr=0 mean=2.512605e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=46 | FIN rc=0 fatal=0 dtr=0 mean=3.067227e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=38 |
runs/cart/nltek2/FIN3_d | D6 rc=0 fatal=0 dtr=0 mean=2.613445e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=54 | FIN rc=0 fatal=0 dtr=0 mean=3.067227e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=46 |
runs/cart/nltek3/FIN3_d | D6 rc=0 fatal=0 dtr=0 mean=2.731092e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=53 | FIN rc=0 fatal=0 dtr=0 mean=3.142857e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=48 |
runs/cart/nltek4/FIN3_d | D6 rc=0 fatal=0 dtr=0 mean=2.848739e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=53 | FIN rc=0 fatal=0 dtr=0 mean=3.100840e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=49 |
runs/cart/nltek5/FIN3_d | D6 rc=0 fatal=0 dtr=0 mean=2.638655e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=51 | FIN rc=0 fatal=0 dtr=0 mean=3.117647e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=47 |
runs/cart/nltek6/FIN3_d | D6 rc=0 fatal=0 dtr=0 mean=2.453782e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=46 | FIN rc=0 fatal=0 dtr=0 mean=3.058824e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=47 |
runs/cart/nltek7/FIN3_d | D6 rc=0 fatal=0 dtr=0 mean=2.260504e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=43 | FIN rc=0 fatal=0 dtr=0 mean=3.042017e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=50 |
runs/cart/nltek8/FIN3_d | D6 rc=0 fatal=0 dtr=0 mean=2.411765e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=47 | FIN rc=0 fatal=0 dtr=0 mean=3.117647e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=48 |
runs/cart/oneo/FIN3_b3k3 | D6 rc=0 fatal=0 dtr=0 mean=2.759494e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 t=400 | FIN rc=0 fatal=0 dtr=0 mean=2.759494e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 t=416 |
runs/cart/oneo/FIN3_b3k5 | D6 rc=0 fatal=0 dtr=0 mean=2.822785e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 t=374 | FIN rc=0 fatal=0 dtr=0 mean=2.822785e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 t=389 |
runs/cart/oneo/FIN3_b3k7 | D6 rc=0 fatal=0 dtr=0 mean=2.873418e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 t=477 | FIN rc=0 fatal=0 dtr=0 mean=2.873418e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 t=501 |
runs/cart/pipe/FIN3_c220 | D6 rc=0 fatal=0 dtr=0 mean=3.004990e+00 max=3.300000e+01 NON-CONVERGED=0.000000e+00 t=21775 | FIN rc=0 fatal=0 dtr=0 mean=2.992982e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 t=20128 |
runs/cart/pipe/FIN3_c2200 | D6 rc=0 fatal=0 dtr=0 mean=3.545455e+00 max=1.000000e+01 NON-CONVERGED=0.000000e+00 t=3666 | FIN rc=0 fatal=0 dtr=0 mean=3.326633e+00 max=1.800000e+01 NON-CONVERGED=0.000000e+00 t=7893 |
runs/cart/rs_m1.2/FIN3_n512 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=969 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=991 |
runs/cart/rs_m2/FIN3_n1024 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=3187 | FIN rc=0 fatal=0 dtr=0 mean=3.003010e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=2662 |
runs/cart/rs_m2/FIN3_n512 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=1037 | FIN rc=0 fatal=0 dtr=0 mean=3.010581e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=933 |
runs/cart/rs_m3/FIN3_n1024 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=2535 | FIN rc=0 fatal=0 dtr=0 mean=3.000093e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=2253 |
runs/cart/rs_m3/FIN3_n512 | D6 rc=0 fatal=0 dtr=0 mean=3.456350e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=909 | FIN rc=0 fatal=0 dtr=0 mean=3.015957e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=830 |
runs/cart/rsx_m1.05/P_n512_c0.4 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=1574 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=1672 |
runs/cart/rsx_m1.2/P_n512_c0.4 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=1100 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=1073 |
runs/cart/rsx_m27/P_n1024_c0.4 | D6 rc=0 fatal=0 dtr=0 mean=3.000250e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 t=8904 | FIN rc=0 fatal=0 dtr=0 mean=3.000250e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 t=9056 |
runs/cart/rsx_m2/P_n512_c0.4 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=1057 | FIN rc=0 fatal=0 dtr=0 mean=3.007828e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=967 |
runs/cart/rsx_m3/P_n512_c0.4 | D6 rc=0 fatal=0 dtr=0 mean=3.429075e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=913 | FIN rc=0 fatal=0 dtr=0 mean=3.003804e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=794 |
runs/cart/rsx_m50/P_n1024_c0.4 | D6 rc=0 fatal=0 dtr=0 mean=4.000095e+00 max=6.000000e+00 NON-CONVERGED=0.000000e+00 t=4238 | FIN rc=0 fatal=0 dtr=0 mean=3.724557e+00 max=6.000000e+00 NON-CONVERGED=0.000000e+00 t=3876 |
runs/cart/rsx_m5/P_n1024_c0.4 | D6 rc=0 fatal=0 dtr=0 mean=3.060809e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=4273 | FIN rc=0 fatal=0 dtr=0 mean=3.076059e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=4337 |
runs/cart/rsx_m5/P_n2048_c0.2 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=12653 | FIN rc=0 fatal=0 dtr=0 mean=3.009586e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=12265 |
runs/cart/rsx_m5/P_n2048_c0.4 | D6 rc=0 fatal=0 dtr=0 mean=3.995498e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=11681 | FIN rc=0 fatal=0 dtr=0 mean=3.493331e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=9758 |
runs/cart/sdifk/FIN3_d | D6 rc=0 fatal=0 dtr=0 mean=2.000000e+00 max=2.000000e+00 NON-CONVERGED=0.000000e+00 t=47 | FIN rc=0 fatal=0 dtr=0 mean=2.000000e+00 max=2.000000e+00 NON-CONVERGED=0.000000e+00 t=52 |
runs/cart/shd/FIN3_r12g2 | D6 rc=0 fatal=0 dtr=0 mean=4.052632e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 t=85 | FIN rc=0 fatal=0 dtr=0 mean=4.000000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=110 |
runs/cart/shd/FIN3_r36g2 | D6 rc=0 fatal=0 dtr=0 mean=4.000000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=94 | FIN rc=0 fatal=0 dtr=0 mean=4.000000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=122 |
runs/cart/shdg/FIN3_g8x32 | D6 rc=0 fatal=0 dtr=0 mean=4.000000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=170 | FIN rc=0 fatal=0 dtr=0 mean=4.000000e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=174 |
runs/cart/slabc/FIN3_n128 | D6 rc=0 fatal=0 dtr=0 mean=2.102362e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=8 | FIN rc=0 fatal=0 dtr=0 mean=2.094488e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=14 |
runs/cart/slabc/FIN3_n256 | D6 rc=0 fatal=0 dtr=0 mean=2.094488e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=13 | FIN rc=0 fatal=0 dtr=0 mean=2.094488e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=17 |
runs/cart/slabc/FIN3_n32 | D6 rc=0 fatal=0 dtr=0 mean=2.039370e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=7 | FIN rc=0 fatal=0 dtr=0 mean=2.039370e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=10 |
runs/cart/slabc/FIN3_n512 | D6 rc=0 fatal=0 dtr=0 mean=2.118110e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=23 | FIN rc=0 fatal=0 dtr=0 mean=2.118110e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=26 |
runs/cart/slabc/FIN3_n64 | D6 rc=0 fatal=0 dtr=0 mean=2.062992e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=8 | FIN rc=0 fatal=0 dtr=0 mean=2.062992e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=11 |
runs_omp/cart/ddf/FIN3_d | D6 rc=0 fatal=0 dtr=0 mean=3.981228e+00 max=1.000000e+01 NON-CONVERGED=0.000000e+00 t=2327 | FIN rc=0 fatal=0 dtr=0 mean=4.000000e+00 max=1.000000e+01 NON-CONVERGED=0.000000e+00 t=1893 |
runs_omp/cart/ddn/FIN3_d | D6 rc=0 fatal=0 dtr=0 mean=3.368488e+00 max=9.000000e+00 NON-CONVERGED=0.000000e+00 t=2335 | FIN rc=0 fatal=0 dtr=0 mean=3.368488e+00 max=9.000000e+00 NON-CONVERGED=0.000000e+00 t=2018 |
runs_omp/cart/rs_m1.2/FIN3_n512 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=1135 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=976 |
runs_omp/cart/rs_m2/FIN3_n1024 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=2671 | FIN rc=0 fatal=0 dtr=0 mean=3.003010e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=1744 |
runs_omp/cart/rs_m2/FIN3_n512 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=1360 | FIN rc=0 fatal=0 dtr=0 mean=3.010581e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=878 |
runs_omp/cart/rs_m3/FIN3_n1024 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=1912 | FIN rc=0 fatal=0 dtr=0 mean=3.000093e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=1589 |
runs_omp/cart/rs_m3/FIN3_n512 | D6 rc=0 fatal=0 dtr=0 mean=3.456350e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=1097 | FIN rc=0 fatal=0 dtr=0 mean=3.015957e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=779 |
runs/sp/rsx_m27/P_n1024_c0.4 | D6 rc=0 fatal=0 dtr=0 mean=3.000225e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 t=26614 | FIN rc=0 fatal=0 dtr=0 mean=3.000225e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 t=26484 |
runs/sp/rsx_m2/P_n512_c0.4 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=8732 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=8659 |
runs/sp/rsx_m5/P_n1024_c0.4 | D6 rc=0 fatal=0 dtr=0 mean=3.000048e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=30779 | FIN rc=0 fatal=0 dtr=0 mean=3.000048e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=30518 |
runs/spx/marsh_r1_u/L128 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=92 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=104 |
runs/spx/marsh_r1_u/L256 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=184 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=197 |
runs/spx/marsh_r1_u/L32 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=26 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=29 |
runs/spx/marsh_r1_u/L512 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=368 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=400 |
runs/spx/marsh_r1_u/L64 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=48 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=53 |
runs/spx/patmT_L16_s/L128 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=396 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=424 |
runs/spx/patmT_L16_s/L256 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=778 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=812 |
runs/spx/patmT_L16_s/L32 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=124 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=128 |
runs/spx/patmT_L16_s/L512 | D6 rc=0 fatal=0 dtr=0 mean=3.006289e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=1507 | FIN rc=0 fatal=0 dtr=0 mean=3.012579e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=1544 |
runs/spx/patmT_L16_s/L64 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=214 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=225 |
runs/spx/rlat/L128 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=166 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=191 |
runs/spx/rlat/L16 | D6 rc=0 fatal=0 dtr=0 mean=3.097561e+00 max=8.000000e+00 NON-CONVERGED=0.000000e+00 t=4 | FIN rc=0 fatal=0 dtr=0 mean=3.097561e+00 max=8.000000e+00 NON-CONVERGED=0.000000e+00 t=5 |
runs/spx/rlat/L256 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=658 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=697 |
runs/spx/rlat/L32 | D6 rc=0 fatal=0 dtr=0 mean=3.029289e+00 max=6.000000e+00 NON-CONVERGED=0.000000e+00 t=13 | FIN rc=0 fatal=0 dtr=0 mean=3.025105e+00 max=6.000000e+00 NON-CONVERGED=0.000000e+00 t=13 |
runs/spx/rlat/L64 | D6 rc=0 fatal=0 dtr=0 mean=3.006263e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=42 | FIN rc=0 fatal=0 dtr=0 mean=3.004175e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=48 |
runs/spx/rlat_om5/L128 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=165 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=187 |
runs/spx/rlat_om5/L16 | D6 rc=0 fatal=0 dtr=0 mean=3.105691e+00 max=9.000000e+00 NON-CONVERGED=0.000000e+00 t=5 | FIN rc=0 fatal=0 dtr=0 mean=3.105691e+00 max=9.000000e+00 NON-CONVERGED=0.000000e+00 t=5 |
runs/spx/rlat_om5/L256 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=659 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=709 |
runs/spx/rlat_om5/L32 | D6 rc=0 fatal=0 dtr=0 mean=3.029289e+00 max=6.000000e+00 NON-CONVERGED=0.000000e+00 t=13 | FIN rc=0 fatal=0 dtr=0 mean=3.025105e+00 max=6.000000e+00 NON-CONVERGED=0.000000e+00 t=15 |
runs/spx/rlat_om5/L64 | D6 rc=0 fatal=0 dtr=0 mean=3.008351e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 t=42 | FIN rc=0 fatal=0 dtr=0 mean=3.006263e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 t=47 |
runs/spx/rlatt/L128 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=741 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=771 |
runs/spx/rlatt/L16 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=13 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=16 |
runs/spx/rlatt/L256 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=2996 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=3016 |
runs/spx/rlatt/L32 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=49 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=51 |
runs/spx/rlatt/L64 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=185 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=197 |
runs/spx/rlatt_om5/L128 | D6 rc=0 fatal=0 dtr=0 mean=3.003135e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=744 | FIN rc=0 fatal=0 dtr=0 mean=3.003135e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=774 |
runs/spx/rlatt_om5/L16 | D6 rc=0 fatal=0 dtr=0 mean=3.048780e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 t=14 | FIN rc=0 fatal=0 dtr=0 mean=3.048780e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 t=16 |
runs/spx/rlatt_om5/L256 | D6 rc=0 fatal=0 dtr=0 mean=3.001045e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=2993 | FIN rc=0 fatal=0 dtr=0 mean=3.001045e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=3031 |
runs/spx/rlatt_om5/L32 | D6 rc=0 fatal=0 dtr=0 mean=3.020921e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 t=48 | FIN rc=0 fatal=0 dtr=0 mean=3.020921e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 t=52 |
runs/spx/rlatt_om5/L64 | D6 rc=0 fatal=0 dtr=0 mean=3.010438e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 t=185 | FIN rc=0 fatal=0 dtr=0 mean=3.010438e+00 max=5.000000e+00 NON-CONVERGED=0.000000e+00 t=200 |
runs/spx/rsw_C/L128 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=19 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=24 |
runs/spx/rsw_C/L256 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=76 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=84 |
runs/spx/rsw_C/L32 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=4 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=3 |
runs/spx/rsw_C/L512 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=295 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=323 |
runs/spx/rsw_C/L64 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=6 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=8 |
runs/spx/rsw/L128 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=12 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=13 |
runs/spx/rsw/L256 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=19 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=22 |
runs/spx/rsw/L32 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=4 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=5 |
runs/spx/rsw/L512 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=66 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=73 |
runs/spx/rsw/L64 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=6 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=7 |
runs/spx/rswt_C/L128 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=101 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=116 |
runs/spx/rswt_C/L256 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=393 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=435 |
runs/spx/rswt_C/L32 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=8 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=9 |
runs/spx/rswt_C/L512 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=1549 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=1595 |
runs/spx/rswt_C/L64 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=27 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=29 |
runs/spx/rswt_dt/L1 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=100 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=110 |
runs/spx/rswt_dt/L16 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=1555 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=1604 |
runs/spx/rswt_dt/L2 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=198 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=224 |
runs/spx/rswt_dt/L4 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=397 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=444 |
runs/spx/rswt_dt/L8 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=787 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=822 |
runs/spx/rswt/L128 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=51 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=58 |
runs/spx/rswt/L256 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=103 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=116 |
runs/spx/rswt/L32 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=15 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=17 |
runs/spx/rswt/L512 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=201 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=217 |
runs/spx/rswt/L64 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=27 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=29 |
runs/spx/rswtv_C/L128 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=103 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=114 |
runs/spx/rswtv_C/L256 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=399 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=433 |
runs/spx/rswtv_C/L32 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=9 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=8 |
runs/spx/rswtv_C/L512 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=1555 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=1585 |
runs/spx/rswtv_C/L64 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=27 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=30 |
runs/spx/rswtv_F_C/L128 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=101 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=111 |
runs/spx/rswtv_F_C/L256 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=394 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=422 |
runs/spx/rswtv_F_C/L32 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=9 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=9 |
runs/spx/rswtv_F_C/L512 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=1554 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=1601 |
runs/spx/rswtv_F_C/L64 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=28 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=31 |
runs/spx/rswtv_F/L128 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=52 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=59 |
runs/spx/rswtv_F/L256 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=99 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=112 |
runs/spx/rswtv_F/L32 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=16 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=16 |
runs/spx/rswtv_F/L512 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=200 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=213 |
runs/spx/rswtv_F/L64 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=27 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=30 |
runs/spx/rswtv/L128 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=52 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=56 |
runs/spx/rswtv/L256 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=105 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=114 |
runs/spx/rswtv/L32 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=15 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=16 |
runs/spx/rswtv/L512 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=202 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=212 |
runs/spx/rswtv/L64 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=27 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=31 |
runs/spx/rswv_C/L128 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=22 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=24 |
runs/spx/rswv_C/L256 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=80 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=88 |
runs/spx/rswv_C/L32 | D6 rc=0 fatal=0 dtr=0 mean=3.015748e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=2 | FIN rc=0 fatal=0 dtr=0 mean=3.015748e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=2 |
runs/spx/rswv_C/L512 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=289 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=320 |
runs/spx/rswv_C/L64 | D6 rc=0 fatal=0 dtr=0 mean=3.003922e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=7 | FIN rc=0 fatal=0 dtr=0 mean=3.007843e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=7 |
runs/spx/rswv_F_C/L128 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=20 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=24 |
runs/spx/rswv_F_C/L256 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=75 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=78 |
runs/spx/rswv_F_C/L32 | D6 rc=0 fatal=0 dtr=0 mean=3.015748e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=3 | FIN rc=0 fatal=0 dtr=0 mean=3.015748e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=3 |
runs/spx/rswv_F_C/L512 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=288 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=322 |
runs/spx/rswv_F_C/L64 | D6 rc=0 fatal=0 dtr=0 mean=3.003922e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=7 | FIN rc=0 fatal=0 dtr=0 mean=3.007843e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=7 |
runs/spx/rswv_F/L128 | D6 rc=0 fatal=0 dtr=0 mean=3.003922e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=11 | FIN rc=0 fatal=0 dtr=0 mean=3.007843e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=13 |
runs/spx/rswv_F/L256 | D6 rc=0 fatal=0 dtr=0 mean=3.003922e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=20 | FIN rc=0 fatal=0 dtr=0 mean=3.007843e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=22 |
runs/spx/rswv_F/L32 | D6 rc=0 fatal=0 dtr=0 mean=3.003922e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=5 | FIN rc=0 fatal=0 dtr=0 mean=3.007843e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=5 |
runs/spx/rswv_F/L512 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=72 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=80 |
runs/spx/rswv_F/L64 | D6 rc=0 fatal=0 dtr=0 mean=3.003922e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=6 | FIN rc=0 fatal=0 dtr=0 mean=3.007843e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=6 |
runs/spx/rswv/L128 | D6 rc=0 fatal=0 dtr=0 mean=3.003922e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=11 | FIN rc=0 fatal=0 dtr=0 mean=3.007843e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=12 |
runs/spx/rswv/L256 | D6 rc=0 fatal=0 dtr=0 mean=3.003922e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=20 | FIN rc=0 fatal=0 dtr=0 mean=3.007843e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=22 |
runs/spx/rswv/L32 | D6 rc=0 fatal=0 dtr=0 mean=3.003922e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=4 | FIN rc=0 fatal=0 dtr=0 mean=3.007843e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=4 |
runs/spx/rswv/L512 | D6 rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=72 | FIN rc=0 fatal=0 dtr=0 mean=3.000000e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=78 |
runs/spx/rswv/L64 | D6 rc=0 fatal=0 dtr=0 mean=3.003922e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=7 | FIN rc=0 fatal=0 dtr=0 mean=3.007843e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=7 |
runs/xtra/g5_k12.8 | D6 rc=0 fatal=0 dtr=0 mean=2.000000e+00 max=2.000000e+00 NON-CONVERGED=0.000000e+00 t=99 | FIN rc=0 fatal=0 dtr=0 mean=2.000000e+00 max=2.000000e+00 NON-CONVERGED=0.000000e+00 t=110 |
runs/xtra/g5_k128 | D6 rc=0 fatal=0 dtr=0 mean=2.000000e+00 max=2.000000e+00 NON-CONVERGED=0.000000e+00 t=10 | FIN rc=0 fatal=0 dtr=0 mean=2.007874e+00 max=3.000000e+00 NON-CONVERGED=0.000000e+00 t=13 |
runs/xtra/g5_k1280 | D6 rc=0 fatal=0 dtr=0 mean=2.000000e+00 max=2.000000e+00 NON-CONVERGED=0.000000e+00 t=46 | FIN rc=0 fatal=0 dtr=0 mean=2.000000e+00 max=2.000000e+00 NON-CONVERGED=0.000000e+00 t=54 |
runs/xtra/g5_k12800 | D6 rc=0 fatal=0 dtr=0 mean=2.000000e+00 max=2.000000e+00 NON-CONVERGED=0.000000e+00 t=398 | FIN rc=0 fatal=0 dtr=0 mean=2.000000e+00 max=2.000000e+00 NON-CONVERGED=0.000000e+00 t=459 |
runs/xtra/shd432 | D6 rc=0 fatal=0 dtr=0 mean=2.250000e+00 max=9.000000e+00 NON-CONVERGED=0.000000e+00 t=4668 | FIN rc=0 fatal=0 dtr=0 mean=2.088608e+00 max=4.000000e+00 NON-CONVERGED=0.000000e+00 t=4620 |
== (2) outputs FIN -> D6 (cmpout.py)
runs/cart/cb/FIN3_r12 DIFF 4 0 0 ndiff=1 first=bin/cb.m1.00001.bin maxrel=1.169e-02
runs/cart/cb/FIN3_r4 DIFF 4 0 0 ndiff=1 first=bin/cb.m1.00001.bin maxrel=1.541e+00
runs/cart/cb/FIN3_r4c0.5 BITWISE 4 0 0
runs/cart/cb/FIN3_r4c32 BITWISE 4 0 0
runs/cart/cbg/FIN3_g8x32 DIFF 4 0 0 ndiff=2 first=bin/cbg.m1.00001.bin maxrel=1.999e+00
runs/cart/ddf/FIN3_d DIFF 0 10 0 ndiff=8 first=tab/ddf.hydro_w.00001.tab maxrel=2.000e+00
runs/cart/ddiff/FIN3_d DIFF 0 5 0 ndiff=4 first=tab/ddiff.m1.00001.tab maxrel=2.000e+00
runs/cart/ddn/FIN3_d DIFF 0 10 0 ndiff=4 first=tab/ddn.m1.00001.tab maxrel=3.223e-07
runs/cart/eq1/FIN3_fine DIFF 82 0 0 ndiff=80 first=bin/eq1.hydro_w.00001.bin maxrel=2.000e+00
runs/cart/eq2/FIN3_fine DIFF 82 0 0 ndiff=80 first=bin/eq2.hydro_w.00001.bin maxrel=2.000e+00
runs/cart/eq3t/FIN3_t DIFF 0 402 0 ndiff=400 first=tab/eq3.hydro_w.00001.tab maxrel=2.000e+00
runs/cart/hcyl/FIN3_g8x32 DIFF 9 0 0 ndiff=2 first=bin/shd.m1.00001.bin maxrel=7.732e-08
runs/cart/marsh/FIN3_c1 DIFF 0 4 0 ndiff=2 first=tab/m1_marshak.hydro_w.00001.tab maxrel=1.986e+00
runs/cart/marsh/FIN3_c10 DIFF 0 4 0 ndiff=2 first=tab/m1_marshak.hydro_w.00001.tab maxrel=1.971e+00
runs/cart/nltek1/FIN3_d DIFF 0 2 0 ndiff=1 first=tab/nlte1.m1.00001.tab maxrel=1.998e+00
runs/cart/nltek2/FIN3_d DIFF 0 2 0 ndiff=1 first=tab/nlte2.m1.00001.tab maxrel=2.000e+00
runs/cart/nltek3/FIN3_d DIFF 0 2 0 ndiff=1 first=tab/nlte3.m1.00001.tab maxrel=2.000e+00
runs/cart/nltek4/FIN3_d DIFF 0 2 0 ndiff=1 first=tab/nlte4.m1.00001.tab maxrel=1.997e+00
runs/cart/nltek5/FIN3_d DIFF 0 2 0 ndiff=1 first=tab/nlte5.m1.00001.tab maxrel=2.000e+00
runs/cart/nltek6/FIN3_d DIFF 0 2 0 ndiff=1 first=tab/nlte6.m1.00001.tab maxrel=1.999e+00
runs/cart/nltek7/FIN3_d DIFF 0 2 0 ndiff=1 first=tab/nlte7.m1.00001.tab maxrel=1.999e+00
runs/cart/nltek8/FIN3_d DIFF 0 2 0 ndiff=1 first=tab/nlte8.m1.00001.tab maxrel=2.000e+00
runs/cart/oneo/FIN3_b3k3 BITWISE 18 0 0
runs/cart/oneo/FIN3_b3k5 BITWISE 18 0 0
runs/cart/oneo/FIN3_b3k7 BITWISE 18 0 0
runs/cart/pipe/FIN3_c220 DIFF 303 0 0 ndiff=300 first=bin/pipe.hydro_w.00001.bin maxrel=2.000e+00
runs/cart/pipe/FIN3_c2200 DIFF 303 0 0 ndiff=299 first=bin/pipe.hydro_w.00001.bin maxrel=2.000e+00
runs/cart/rs_m1.2/FIN3_n512 DIFF 0 14 0 ndiff=12 first=tab/m1_radshock.hydro_w.00001.tab maxrel=2.000e+00
runs/cart/rs_m2/FIN3_n1024 DIFF 0 14 0 ndiff=12 first=tab/m1_radshock.hydro_w.00001.tab maxrel=1.999e+00
runs/cart/rs_m2/FIN3_n512 DIFF 0 14 0 ndiff=12 first=tab/m1_radshock.hydro_w.00001.tab maxrel=2.000e+00
runs/cart/rs_m3/FIN3_n1024 DIFF 0 14 0 ndiff=12 first=tab/m1_radshock.hydro_w.00001.tab maxrel=2.000e+00
runs/cart/rs_m3/FIN3_n512 DIFF 0 14 0 ndiff=12 first=tab/m1_radshock.hydro_w.00001.tab maxrel=2.000e+00
runs/cart/rsx_m1.05/P_n512_c0.4 DIFF 0 8 0 ndiff=6 first=tab/m1_radshock.hydro_w.00001.tab maxrel=2.000e+00
runs/cart/rsx_m1.2/P_n512_c0.4 DIFF 0 8 0 ndiff=6 first=tab/m1_radshock.hydro_w.00001.tab maxrel=2.000e+00
runs/cart/rsx_m27/P_n1024_c0.4 DIFF 0 8 0 ndiff=6 first=tab/m1_radshock.hydro_w.00001.tab maxrel=1.999e+00
runs/cart/rsx_m2/P_n512_c0.4 DIFF 0 8 0 ndiff=6 first=tab/m1_radshock.hydro_w.00001.tab maxrel=1.999e+00
runs/cart/rsx_m3/P_n512_c0.4 DIFF 0 8 0 ndiff=6 first=tab/m1_radshock.hydro_w.00001.tab maxrel=2.000e+00
runs/cart/rsx_m50/P_n1024_c0.4 DIFF 0 8 0 ndiff=6 first=tab/m1_radshock.hydro_w.00001.tab maxrel=1.999e+00
runs/cart/rsx_m5/P_n1024_c0.4 DIFF 0 8 0 ndiff=6 first=tab/m1_radshock.hydro_w.00001.tab maxrel=1.998e+00
runs/cart/rsx_m5/P_n2048_c0.2 DIFF 0 8 0 ndiff=6 first=tab/m1_radshock.hydro_w.00001.tab maxrel=2.000e+00
runs/cart/rsx_m5/P_n2048_c0.4 DIFF 0 8 0 ndiff=6 first=tab/m1_radshock.hydro_w.00001.tab maxrel=2.000e+00
runs/cart/sdifk/FIN3_d BITWISE 0 4 0
runs/cart/shd/FIN3_r12g2 DIFF 4 0 0 ndiff=2 first=bin/shd.m1.00001.bin maxrel=2.000e+00
runs/cart/shd/FIN3_r36g2 DIFF 4 0 0 ndiff=2 first=bin/shd.m1.00001.bin maxrel=2.000e+00
runs/cart/shdg/FIN3_g8x32 DIFF 4 0 0 ndiff=2 first=bin/shdg.m1.00001.bin maxrel=2.000e+00
runs/cart/slabc/FIN3_n128 DIFF 6 0 0 ndiff=1 first=bin/slab.m1.00001.bin maxrel=1.994e+00
runs/cart/slabc/FIN3_n256 DIFF 6 0 0 ndiff=1 first=bin/slab.m1.00001.bin maxrel=1.997e+00
runs/cart/slabc/FIN3_n32 BITWISE 9 0 0
runs/cart/slabc/FIN3_n512 DIFF 6 0 0 ndiff=1 first=bin/slab.m1.00001.bin maxrel=1.998e+00
runs/cart/slabc/FIN3_n64 DIFF 9 0 0 ndiff=2 first=bin/slab.m1.00001.bin maxrel=2.000e+00
runs_omp/cart/ddf/FIN3_d DIFF 0 10 0 ndiff=8 first=tab/ddf.hydro_w.00001.tab maxrel=2.000e+00
runs_omp/cart/ddn/FIN3_d DIFF 0 10 0 ndiff=4 first=tab/ddn.m1.00001.tab maxrel=3.223e-07
runs_omp/cart/rs_m1.2/FIN3_n512 DIFF 0 14 0 ndiff=12 first=tab/m1_radshock.hydro_w.00001.tab maxrel=2.000e+00
runs_omp/cart/rs_m2/FIN3_n1024 DIFF 0 14 0 ndiff=12 first=tab/m1_radshock.hydro_w.00001.tab maxrel=1.999e+00
runs_omp/cart/rs_m2/FIN3_n512 DIFF 0 14 0 ndiff=12 first=tab/m1_radshock.hydro_w.00001.tab maxrel=2.000e+00
runs_omp/cart/rs_m3/FIN3_n1024 DIFF 0 14 0 ndiff=12 first=tab/m1_radshock.hydro_w.00001.tab maxrel=2.000e+00
runs_omp/cart/rs_m3/FIN3_n512 DIFF 0 14 0 ndiff=12 first=tab/m1_radshock.hydro_w.00001.tab maxrel=2.000e+00
runs/sp/rsx_m27/P_n1024_c0.4 DIFF 0 8 0 ndiff=6 first=tab/m1_radshock.hydro_w.00001.tab maxrel=2.000e+00
runs/sp/rsx_m2/P_n512_c0.4 DIFF 0 8 0 ndiff=6 first=tab/m1_radshock.hydro_w.00001.tab maxrel=2.000e+00
runs/sp/rsx_m5/P_n1024_c0.4 DIFF 0 8 0 ndiff=6 first=tab/m1_radshock.hydro_w.00001.tab maxrel=2.000e+00
runs/spx/marsh_r1_u/L128 DIFF 0 15 0 ndiff=12 first=tab/sp.hydro_u.00001.tab maxrel=1.001e+00
runs/spx/marsh_r1_u/L256 DIFF 0 15 0 ndiff=12 first=tab/sp.hydro_u.00001.tab maxrel=1.987e+00
runs/spx/marsh_r1_u/L32 DIFF 0 15 0 ndiff=12 first=tab/sp.hydro_u.00001.tab maxrel=1.976e+00
runs/spx/marsh_r1_u/L512 DIFF 0 15 0 ndiff=12 first=tab/sp.hydro_u.00001.tab maxrel=1.205e+00
runs/spx/marsh_r1_u/L64 DIFF 0 15 0 ndiff=12 first=tab/sp.hydro_u.00001.tab maxrel=1.951e+00
runs/spx/patmT_L16_s/L128 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.996e+00
runs/spx/patmT_L16_s/L256 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.982e+00
runs/spx/patmT_L16_s/L32 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.984e+00
runs/spx/patmT_L16_s/L512 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.997e+00
runs/spx/patmT_L16_s/L64 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.998e+00
runs/spx/rlat/L128 DIFF 27 0 0 ndiff=24 first=bin/sp.hydro_u.00001.bin maxrel=7.926e-03
runs/spx/rlat/L16 DIFF 27 0 0 ndiff=24 first=bin/sp.hydro_u.00001.bin maxrel=9.895e-02
runs/spx/rlat/L256 DIFF 27 0 0 ndiff=24 first=bin/sp.hydro_u.00001.bin maxrel=1.377e-03
runs/spx/rlat/L32 DIFF 27 0 0 ndiff=24 first=bin/sp.hydro_u.00001.bin maxrel=3.417e-01
runs/spx/rlat/L64 DIFF 27 0 0 ndiff=24 first=bin/sp.hydro_u.00001.bin maxrel=1.917e-01
runs/spx/rlat_om5/L128 DIFF 27 0 0 ndiff=24 first=bin/sp.hydro_u.00001.bin maxrel=5.214e-03
runs/spx/rlat_om5/L16 DIFF 27 0 0 ndiff=24 first=bin/sp.hydro_u.00001.bin maxrel=3.577e-02
runs/spx/rlat_om5/L256 DIFF 27 0 0 ndiff=24 first=bin/sp.hydro_u.00001.bin maxrel=5.085e-03
runs/spx/rlat_om5/L32 DIFF 27 0 0 ndiff=24 first=bin/sp.hydro_u.00001.bin maxrel=8.447e-03
runs/spx/rlat_om5/L64 DIFF 27 0 0 ndiff=24 first=bin/sp.hydro_u.00001.bin maxrel=2.614e-03
runs/spx/rlatt/L128 DIFF 27 0 0 ndiff=24 first=bin/sp.hydro_u.00001.bin maxrel=1.687e-04
runs/spx/rlatt/L16 DIFF 27 0 0 ndiff=24 first=bin/sp.hydro_u.00001.bin maxrel=1.233e-04
runs/spx/rlatt/L256 DIFF 27 0 0 ndiff=24 first=bin/sp.hydro_u.00001.bin maxrel=4.016e-03
runs/spx/rlatt/L32 DIFF 27 0 0 ndiff=24 first=bin/sp.hydro_u.00001.bin maxrel=3.953e-03
runs/spx/rlatt/L64 DIFF 27 0 0 ndiff=24 first=bin/sp.hydro_u.00001.bin maxrel=3.862e-04
runs/spx/rlatt_om5/L128 DIFF 27 0 0 ndiff=24 first=bin/sp.hydro_u.00001.bin maxrel=1.016e-01
runs/spx/rlatt_om5/L16 DIFF 27 0 0 ndiff=24 first=bin/sp.hydro_u.00001.bin maxrel=1.151e-02
runs/spx/rlatt_om5/L256 DIFF 27 0 0 ndiff=24 first=bin/sp.hydro_u.00001.bin maxrel=8.993e-01
runs/spx/rlatt_om5/L32 DIFF 27 0 0 ndiff=24 first=bin/sp.hydro_u.00001.bin maxrel=1.661e-02
runs/spx/rlatt_om5/L64 DIFF 27 0 0 ndiff=24 first=bin/sp.hydro_u.00001.bin maxrel=2.796e-01
runs/spx/rsw_C/L128 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.994e+00
runs/spx/rsw_C/L256 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.982e+00
runs/spx/rsw_C/L32 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.889e+00
runs/spx/rsw_C/L512 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.997e+00
runs/spx/rsw_C/L64 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.952e+00
runs/spx/rsw/L128 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.999e+00
runs/spx/rsw/L256 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.992e+00
runs/spx/rsw/L32 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.992e+00
runs/spx/rsw/L512 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=2.000e+00
runs/spx/rsw/L64 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.952e+00
runs/spx/rswt_C/L128 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.996e+00
runs/spx/rswt_C/L256 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.995e+00
runs/spx/rswt_C/L32 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=2.000e+00
runs/spx/rswt_C/L512 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.997e+00
runs/spx/rswt_C/L64 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.997e+00
runs/spx/rswt_dt/L1 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.994e+00
runs/spx/rswt_dt/L16 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=2.000e+00
runs/spx/rswt_dt/L2 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.994e+00
runs/spx/rswt_dt/L4 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.995e+00
runs/spx/rswt_dt/L8 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.999e+00
runs/spx/rswt/L128 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.999e+00
runs/spx/rswt/L256 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.994e+00
runs/spx/rswt/L32 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.989e+00
runs/spx/rswt/L512 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.999e+00
runs/spx/rswt/L64 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.997e+00
runs/spx/rswtv_C/L128 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.998e+00
runs/spx/rswtv_C/L256 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.999e+00
runs/spx/rswtv_C/L32 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.975e+00
runs/spx/rswtv_C/L512 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.998e+00
runs/spx/rswtv_C/L64 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.983e+00
runs/spx/rswtv_F_C/L128 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.994e+00
runs/spx/rswtv_F_C/L256 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.996e+00
runs/spx/rswtv_F_C/L32 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.940e+00
runs/spx/rswtv_F_C/L512 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=2.000e+00
runs/spx/rswtv_F_C/L64 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.991e+00
runs/spx/rswtv_F/L128 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.991e+00
runs/spx/rswtv_F/L256 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.985e+00
runs/spx/rswtv_F/L32 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.994e+00
runs/spx/rswtv_F/L512 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.999e+00
runs/spx/rswtv_F/L64 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.991e+00
runs/spx/rswtv/L128 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.977e+00
runs/spx/rswtv/L256 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.996e+00
runs/spx/rswtv/L32 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.972e+00
runs/spx/rswtv/L512 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.990e+00
runs/spx/rswtv/L64 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.983e+00
runs/spx/rswv_C/L128 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.996e+00
runs/spx/rswv_C/L256 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.993e+00
runs/spx/rswv_C/L32 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.990e+00
runs/spx/rswv_C/L512 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.998e+00
runs/spx/rswv_C/L64 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.989e+00
runs/spx/rswv_F_C/L128 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.994e+00
runs/spx/rswv_F_C/L256 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.996e+00
runs/spx/rswv_F_C/L32 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.941e+00
runs/spx/rswv_F_C/L512 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.998e+00
runs/spx/rswv_F_C/L64 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.915e+00
runs/spx/rswv_F/L128 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.997e+00
runs/spx/rswv_F/L256 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.992e+00
runs/spx/rswv_F/L32 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.972e+00
runs/spx/rswv_F/L512 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=2.000e+00
runs/spx/rswv_F/L64 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.915e+00
runs/spx/rswv/L128 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.995e+00
runs/spx/rswv/L256 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.996e+00
runs/spx/rswv/L32 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.900e+00
runs/spx/rswv/L512 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.992e+00
runs/spx/rswv/L64 DIFF 0 6 0 ndiff=3 first=tab/sp.hydro_u.00001.tab maxrel=1.989e+00
runs/xtra/g5_k12.8 DIFF 21 0 0 ndiff=20 first=bin/m1_thick_pulse.m1.00001.bin maxrel=2.000e+00
runs/xtra/g5_k128 DIFF 21 0 0 ndiff=20 first=bin/m1_thick_pulse.m1.00001.bin maxrel=2.000e+00
runs/xtra/g5_k1280 BITWISE 21 0 0
runs/xtra/g5_k12800 BITWISE 21 0 0
runs/xtra/shd432 DIFF 2 0 0 ndiff=1 first=bin/sp.m1.00001.bin maxrel=2.000e+00
```

## Verdict idortp1118 (prod-1118 1d1972805 + EXTRA_HG)
- (a) waves: **PASS**, q = d6int bitwise (no key fires).
- (b) 5 sp/cs cells: **PASS** (metrics = D6i to 4 digits; P18 = P18H; the shadow cells differ from D6i at 1e-8 to 5e-7 rel).
- (c) Raven set: **correctness PASS 160/160** (rc / FATAL / NC / DT_RETRY / Picard vs FIN); **q != d6int on 126 rows, mostly from
  time2_stage_local_be=false** (97 rows with no other acting key); HG keys: adm_poslim fires in 2 rows, lag_local and dtr_poslim never.
  The accuracy-metric rule needs your analyzers on the kept outputs.
