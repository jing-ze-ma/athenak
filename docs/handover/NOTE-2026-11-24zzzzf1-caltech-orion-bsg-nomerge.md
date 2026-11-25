# NOTE 2026-11-24zzzzf1 caltech -> orion: bsg without merged cells (your 24zzzzf1 / zzzzf3 / zzzzf7)

## NM0 did not start: question
- NM0 (c4 3ae130e9, 02-09, 11-24 19:46): the c1 production deck without its three merge_lat lines stops at start-up:
```
### FATAL ERROR in /resnick/groups/carnegie_poc/jingze/wt_orion1121/c4/src/pgen/he_star_vet.cpp at line 3525
he_star_vet: problem/he_ic_profile_seed_group needs hydro/merge_lat_rc
```
- `problem/he_ic_profile_seed_group` in the bsg deck needs `hydro/merge_lat_rc`. Deck lines concerned:
```
338:he_ic_profile_seed_group = true   # lhxseed-1106: one deviate per merged group, as his coarse cells
```
- **How should NM0 be started?** For example merge_lat_f = 1,1,1 with rc kept (as NM1), or another value of the seed key, or the seed key removed: each changes the start state, so I edit nothing until you say. MG0 (control, fresh) ran clean: 100 cycles, CW / FATAL / DT_RETRY 0.
- K1 (G): dropped (zzzzf9), no repeat; its one result is in the zzzzd1 NOTE. rst cycles of the c1 bsg: 00001 c3575, 00002 c7300, 00003 c10883, 00004 c14605; none at c5000-5500, MG1c skipped.

## NM0 (fresh, c4, merge_lat_f = 1,1,1 with rc kept, zzzzg1) and MG0 (as production)
- NM0 is run as NM0b: the c1 production deck with merge_lat_f = 1,1,1, dt_report = 10, final_restart = false; MG0: the same deck with merge_lat_f = 2,1,2. The seed noise of the initial state differs between the two (the seed group is set by the merge layout), so their early rows differ by construction.
- NM0b: 2026-11-24 20:19:36 rc=0 (srun) ended: Terminating on cycle limit wall=116 s CW 0 FATAL 0 inf 0 DT_RETRY 0 stage fallback 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- MG0: 2026-11-24 19:27:49 rc=0 (srun) ended: Terminating on cycle limit wall=129 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
/resnick/groups/carnegie_poc/jingze/b1_1124/nm0.py:33: RuntimeWarning: invalid value encountered in true_divide
  r = a[m, n[c]] / np.interp(t, b[:, 0], b[:, D[y][2][c]]) - 1
| arm | rows | dt min/mean/max | s/cyc (last two elapsed prints) | its t2s2 (sum, solves) | t2s3 | CW/FATAL/inf/DT_RETRY | fallback/BE | floor_pay sum | nfb (last MLATDIAG) | wall per 1000 s model | L_top/L_in end | KE_int end |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| NM0b | 101 | 82.253/83.294/83.816 | 0.7754 (c50-100) | 17.8, 3.03 | 22.9, 3.62 | 0/0/0/0 | 0/1 | 0.000e+00 | 0 | 8.47 | 1.062339 | 7.24744e+40 |
| MG0 | 101 | 164.644/167.310/168.458 | 0.8070 (c50-100) | 20.5, 3.26 | 25.5, 3.94 | 0/0/0/0 | 0/1 | 0.000e+00 | 0 | 4.45 | 1.040158 | 9.26537e+40 |
- NM0b: MLAT lines: ['MLAT merged lateral cells: 3 shell(s): x1v < 2.43495e+12 f 1; x1v < 3.82635e+12 f 1; x1v < 5.5656e+12 f 1; mom rigid, rad fine', 'MLAT init: 0 non-uniform gas groups, 0 non-uniform radiation groups -> state kept']
  DTREP 22 lines; binders {('3', '33,130,3', '3'): 10, ('0', '5,3,3', '3'): 5, ('0', '70,3,3', '3'): 3, ('3', '130,130,3', '3'): 2}; last r 1.40233e+12 rho 5.64582e-05 v 15795.6 ceff 1.77996e+07 h 4.97371e+09 cfl*dt_hyd 83.754; min cfl*dt_hyd 82.4384 max 83.754
- MG0: MLAT lines: ['MLAT merged lateral cells: 3 shell(s): x1v < 2.43495e+12 f 2; x1v < 3.82635e+12 f 1; x1v < 5.5656e+12 f 2; mom rigid, rad fine', 'MLAT init: 0 non-uniform gas groups, 0 non-uniform radiation groups -> state kept']
  DTREP 22 lines; binders {('0', '33,3,3', '3'): 11, ('0', '53,3,3', '3'): 3, ('0', '35,3,3', '3'): 3, ('0', '42,3,3', '3'): 2}; last r 1.402326e+12 rho 5.729379e-05 v -4.191858e+04 ceff 1.768455e+07 h 4.973708e+09 cfl*dt_hyd 1.683485e+02; min cfl*dt_hyd 165.038 max 168.3485
- NM0b vs MG0 at equal model time (101 rows of NM0b; interpolated):
  - L_top rel: row 1 -2.784e-01, row 5 -9.261e-01, row 10 -1.061e-01, last -3.876e-04; max rows 1-10 9.551e-01, max rows 11-end 7.765e-02
  - KE_int rel: row 1 8.683e-01, row 5 8.181e-01, row 10 5.054e-02, last -1.705e-01; max rows 1-10 1.813e+00, max rows 11-end 1.705e-01
  - Lz_tot rel: row 1 -2.044e+00, row 5 -2.041e+00, row 10 -2.090e+00, last -6.115e-01; max rows 1-10 2.090e+00, max rows 11-end 3.156e+02
  - Etot rel: row 1 -6.281e-05, row 5 7.029e-06, row 10 8.191e-06, last 8.219e-06; max rows 1-10 6.281e-05, max rows 11-end 8.219e-06
  - M_tot rel: row 1 -1.636e-05, row 5 -1.636e-05, row 10 -1.636e-05, last -1.636e-05; max rows 1-10 1.636e-05, max rows 11-end 1.636e-05

## bsg without merged cells (11-24 20:24); binary athena_c4_a000d803c 3ae130e97b13e5bdfa4a7f9a14370b9f, 02-09, 4 ranks, 100 cycles per arm; nothing here goes into production
- rsts of the c1 bsg: 00001 c3575, 00002 c7300, 00003 c10883, 00004 c14605. **No rst at c5000-5500, so MG1c is skipped.** The pair NM1 / MG1 starts from bsg3d.00004.rst (c 14605, the newest, hard-linked).
- **NM0 first try did not start:** `he_star_vet: problem/he_ic_profile_seed_group needs hydro/merge_lat_rc` (he_star_vet.cpp 3525: problem/he_ic_profile_seed_group needs hydro/merge_lat_rc). The bsg deck names that seed key, so removing the merge lines alone is not enough; nothing edited, waiting for you.
- One-cycle test of dt_fixed on the restart (NMck): 2026-11-24 20:20:33 rc=0 (srun) ended: Terminating on cycle limit wall=28 s CW 0 FATAL 0 inf 0 DT_RETRY 0 stage fallback 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps ; dt of that cycle 8.0000000000000000e+01 -> accepted: yes
- NM0b: 2026-11-24 20:19:36 rc=0 (srun) ended: Terminating on cycle limit wall=116 s CW 0 FATAL 0 inf 0 DT_RETRY 0 stage fallback 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- MG0: 2026-11-24 19:27:49 rc=0 (srun) ended: Terminating on cycle limit wall=129 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- NM1: 2026-11-24 20:22:08 rc=0 (srun) ended: Terminating on cycle limit wall=94 s CW 0 FATAL 0 inf 0 DT_RETRY 0 stage fallback 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- MG1: 2026-11-24 20:24:08 rc=0 (srun) ended: Terminating on cycle limit wall=103 s CW 0 FATAL 0 inf 0 DT_RETRY 0 stage fallback 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 

/resnick/groups/carnegie_poc/jingze/b1_1124/nm.py:33: RuntimeWarning: invalid value encountered in true_divide
  r = a[m, n[c]] / np.interp(t, b[:, 0], b[:, D[y][2][c]]) - 1
| arm | rows | dt min/mean/max | s/cyc (last two elapsed prints) | its t2s2 (sum, solves) | t2s3 | CW/FATAL/inf/DT_RETRY | fallback/BE | floor_pay sum | nfb (last MLATDIAG) | wall per 1000 s model | L_top/L_in end | KE_int end |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| NM0b | 101 | 82.253/83.294/83.816 | 0.7754 (c50-100) | 17.8, 3.03 | 22.9, 3.62 | 0/0/0/0 | 0/1 | 0.000e+00 | 0 | 8.47 | 1.062339 | 7.24744e+40 |
| MG0 | 101 | 164.644/167.310/168.458 | 0.8070 (c50-100) | 20.5, 3.26 | 25.5, 3.94 | 0/0/0/0 | 0/1 | 0.000e+00 | 0 | 4.45 | 1.040158 | 9.26537e+40 |
| NM1 | 100 | 80.000/80.000/80.000 | 0.6494 (c14650-14700) | 17.1, 3.03 | 11.0, 2.71 | 0/0/0/0 | 0/0 | 0.000e+00 | - | 8.45 | 1.022834 | 1.21731e+42 |
| MG1 | 100 | 169.376/169.388/169.398 | 0.6989 (c14650-14700) | 24.0, 3.87 | 12.9, 2.80 | 0/0/0/0 | 0/0 | 0.000e+00 | - | 4.34 | 1.018221 | 1.22080e+42 |
- NM0b: MLAT lines: ['MLAT merged lateral cells: 3 shell(s): x1v < 2.43495e+12 f 1; x1v < 3.82635e+12 f 1; x1v < 5.5656e+12 f 1; mom rigid, rad fine', 'MLAT init: 0 non-uniform gas groups, 0 non-uniform radiation groups -> state kept']
  DTREP 22 lines; binders {('3', '33,130,3', '3'): 10, ('0', '5,3,3', '3'): 5, ('0', '70,3,3', '3'): 3, ('3', '130,130,3', '3'): 2}; last r 1.40233e+12 rho 5.64582e-05 v 15795.6 ceff 1.77996e+07 h 4.97371e+09 cfl*dt_hyd 83.754; min cfl*dt_hyd 82.4384 max 83.754
- MG0: MLAT lines: ['MLAT merged lateral cells: 3 shell(s): x1v < 2.43495e+12 f 2; x1v < 3.82635e+12 f 1; x1v < 5.5656e+12 f 2; mom rigid, rad fine', 'MLAT init: 0 non-uniform gas groups, 0 non-uniform radiation groups -> state kept']
  DTREP 22 lines; binders {('0', '33,3,3', '3'): 11, ('0', '53,3,3', '3'): 3, ('0', '35,3,3', '3'): 3, ('0', '42,3,3', '3'): 2}; last r 1.402326e+12 rho 5.729379e-05 v -4.191858e+04 ceff 1.768455e+07 h 4.973708e+09 cfl*dt_hyd 1.683485e+02; min cfl*dt_hyd 165.038 max 168.3485
- NM1: MLAT lines: ['MLAT merged lateral cells: 3 shell(s): x1v < 2.43495e+12 f 1; x1v < 3.82635e+12 f 1; x1v < 5.5656e+12 f 1; mom rigid, rad fine', 'MLAT init: 0 non-uniform gas groups, 0 non-uniform radiation groups -> state kept']
  DTREP 20 lines; binders {('2', '14,3,3', '3'): 20}; last r 1.40233e+12 rho 5.77758e-05 v 27769.2 ceff 1.75898e+07 h 4.97371e+09 cfl*dt_hyd 84.6946; min cfl*dt_hyd 84.6879 max 84.6946
- MG1: MLAT lines: ['MLAT merged lateral cells: 3 shell(s): x1v < 2.43495e+12 f 2; x1v < 3.82635e+12 f 1; x1v < 5.5656e+12 f 2; mom rigid, rad fine', 'MLAT init: 0 non-uniform gas groups, 0 non-uniform radiation groups -> state kept']
  DTREP 20 lines; binders {('2', '16,3,3', '3'): 8, ('2', '14,3,3', '3'): 6, ('2', '79,3,3', '3'): 4, ('3', '77,130,3', '3'): 2}; last r 1.40233e+12 rho 5.77599e-05 v -24784.1 ceff 1.7593e+07 h 4.97371e+09 cfl*dt_hyd 169.387; min cfl*dt_hyd 169.377 max 169.397
- NM0b vs MG0 at equal model time (101 rows of NM0b; interpolated):
  - L_top rel: row 1 -2.784e-01, row 5 -9.261e-01, row 10 -1.061e-01, last -3.876e-04; max rows 1-10 9.551e-01, max rows 11-end 7.765e-02
  - KE_int rel: row 1 8.683e-01, row 5 8.181e-01, row 10 5.054e-02, last -1.705e-01; max rows 1-10 1.813e+00, max rows 11-end 1.705e-01
  - Lz_tot rel: row 1 -2.044e+00, row 5 -2.041e+00, row 10 -2.090e+00, last -6.115e-01; max rows 1-10 2.090e+00, max rows 11-end 3.156e+02
  - Etot rel: row 1 -6.281e-05, row 5 7.029e-06, row 10 8.191e-06, last 8.219e-06; max rows 1-10 6.281e-05, max rows 11-end 8.219e-06
  - M_tot rel: row 1 -1.636e-05, row 5 -1.636e-05, row 10 -1.636e-05, last -1.636e-05; max rows 1-10 1.636e-05, max rows 11-end 1.636e-05
- NM1 vs MG1 at equal model time (100 rows of NM1; interpolated):
  - L_top rel: row 1 5.843e-04, row 5 2.124e-04, row 10 2.137e-04, last 4.423e-05; max rows 1-10 5.843e-04, max rows 11-end 2.322e-04
  - KE_int rel: row 1 -1.867e-04, row 5 -3.394e-04, row 10 -3.691e-04, last 1.272e-03; max rows 1-10 3.707e-04, max rows 11-end 1.272e-03
  - Lz_tot rel: row 1 4.446e-06, row 5 1.777e-06, row 10 1.629e-06, last 2.322e-05; max rows 1-10 4.446e-06, max rows 11-end 3.021e-05
  - Etot rel: row 1 5.755e-10, row 5 -4.961e-10, row 10 -7.381e-10, last -2.745e-09; max rows 1-10 7.381e-10, max rows 11-end 2.745e-09
  - M_tot rel: row 1 1.554e-15, row 5 -1.554e-15, row 10 2.220e-15, last -2.998e-15; max rows 1-10 3.109e-15, max rows 11-end 9.770e-15

### r(c) = dLz_tot + int Lz_hflx dt per cycle (zzzzf7)
| arm | rows | r RMS | r max abs | r sum | dLz_tot total | int Lz_hflx dt | Lz_tot start |
|---|---|---|---|---|---|---|---|
| NM1 | 100 | 2.0644e+33 | 5.6787e+33 | -6.0151e+34 | -1.3079e+37 | 1.3019e+37 | 2.571499e+39 |
  - NM1: CLEAN LZB lines 100; first: CLEAN LZB cycle 14605 t2s 3 ledger -3.6115369954705172e+42 | last: CLEAN LZB cycle 14704 t2s 3 ledger -3.5702695744095515e+42
  - NM1: CLEAN LZSPIN lines 0
  - NM1: CLEAN FCONE_LZ lines 0
  - NM1: CLEAN FCONE lines 200; first: CLEAN FCONE cycle 14605 t2s 2 faces 337383 momentum cut sum |dF| V/c^2 7.981077e+25 | last: CLEAN FCONE cycle 14704 t2s 3 faces 247341 momentum cut sum |dF| V/c^2 4.309105e+25
  - NM1: MLATDIAG lines 0
| MG1 | 100 | 4.7875e+34 | 1.5824e+35 | 3.6663e+35 | -2.7270e+37 | 2.7637e+37 | 2.571356e+39 |
  - MG1: CLEAN LZB lines 100; first: CLEAN LZB cycle 14605 t2s 3 ledger -3.6109082746618638e+42 | last: CLEAN LZB cycle 14704 t2s 3 ledger -3.5690425565830430e+42
  - MG1: CLEAN LZSPIN lines 0
  - MG1: CLEAN FCONE_LZ lines 0
  - MG1: CLEAN FCONE lines 200; first: CLEAN FCONE cycle 14605 t2s 2 faces 360482 momentum cut sum |dF| V/c^2 1.078604e+26 | last: CLEAN FCONE cycle 14704 t2s 3 faces 241385 momentum cut sum |dF| V/c^2 5.843416e+25
  - MG1: MLATDIAG lines 0
- NM1 vs MG1 at equal model time:
  - L_top max rel 5.843e-04
  - L_mid max rel 2.240e-05
  - L_int max rel 1.754e-03
  - KE_int max rel 1.272e-03
  - NM1 s/cyc 0.6494; dt min/max 80.000/80.000; CW 0 FATAL 0 DT_RETRY 0; e_floor sum 0; floor_pay sum 0.000e+00
  - MG1 s/cyc 0.6989; dt min/max 169.375/169.398; CW 0 FATAL 0 DT_RETRY 0; e_floor sum 0; floor_pay sum 0.000e+00
      2 MLAT init: 0 non-uniform gas groups, 0 non-uniform radiation groups -> state kept
      1 MLAT merged lateral cells: 3 shell(s): x1v < 2.43495e+12 f 1; x1v < 3.82635e+12 f 1; x1v < 5.5656e+12 f 1; mom rigid, rad fine
      1 MLAT merged lateral cells: 3 shell(s): x1v < 2.43495e+12 f 2; x1v < 3.82635e+12 f 1; x1v < 5.5656e+12 f 2; mom rigid, rad fine
- Files: caltech_results/BSG_NM/.
