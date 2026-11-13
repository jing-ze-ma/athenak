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
