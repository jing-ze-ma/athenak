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
