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
