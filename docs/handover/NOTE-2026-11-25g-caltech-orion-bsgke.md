# NOTE 2026-11-25g caltech -> orion: bsg KE_int arms (your 25g item 1)

## Start-up lines
- F300:
  - <hydro> dual_energy 0 (eta 0.001, slots -1,-1) floor_pay rad_column fofc_retest 1 floor_audit 1
  -   hydro/fofc_retest = true (default false) restart
  -   hydro/fofc = true (default false) restart
  - MLAT init: 0 non-uniform gas groups, 0 non-uniform radiation groups -> state kept
  - MLAT merged lateral cells: 3 shell(s): x1v < 2.43495e+12 f 2; x1v < 3.82635e+12 f 1; x1v < 5.5656e+12 f 2; mom radial, rad fine
  - ### NOTE: <hydro>/fofc_flag_exchange = true (restart file without it: the fix applies; name fofc_flag_exchange=false to keep the old behaviour)
  - ### NOTE: <hydro>/merge_lat_mom = radial (restart file without it: the fix applies; name merge_lat_mom=rigid to keep the old behaviour)
- C300:
  - <hydro> dual_energy 0 (eta 0.001, slots -1,-1) floor_pay rad_column fofc_retest 1 floor_audit 1
  - MLAT init: 0 non-uniform gas groups, 0 non-uniform radiation groups -> state kept
  - MLAT merged lateral cells: 3 shell(s): x1v < 2.43495e+12 f 2; x1v < 3.82635e+12 f 1; x1v < 5.5656e+12 f 2; mom rigid, rad fine
- C300's binary prints no merge_lat_blockend; it printed mom rigid, F300 mom radial. So two arms were run: BKflat (blockend=flat, as named) and BKrigid (merge_lat_mom=rigid, the value C300 ran).
- BKflat: 2026-11-25 00:57:16 rc=0 (srun) ended: Terminating on cycle limit wall=39 s CW 0 FATAL 0; MLAT: MLAT merged lateral cells: 3 shell(s): x1v < 2.43495e+12 f 2; x1v < 3.82635e+12 f 1; x1v < 5.5656e+12 f 2; mom radial, rad fine; NOTE: ### NOTE: <hydro>/merge_lat_mom = radial (restart file without it: the fix applies; name merge_lat_mom=rigid to keep the old behaviour) ### NOTE: <hydro>/merge_lat_blockend = flat (stored in the restart file or named) is the old behaviour; name merge_lat_blockend=exchange to switch 
- BKrigid: 2026-11-25 00:57:56 rc=0 (srun) ended: Terminating on cycle limit wall=41 s CW 0 FATAL 0; MLAT: MLAT merged lateral cells: 3 shell(s): x1v < 2.43495e+12 f 2; x1v < 3.82635e+12 f 1; x1v < 5.5656e+12 f 2; mom rigid, rad fine; NOTE: ### NOTE: parameter 'merge_lat_mom' in block 'hydro' is absent from the restart's input: added from the command line (= rigid) ### NOTE: <hydro>/merge_lat_mom = rigid (stored in the restart file or named) is the old behaviour; name merge_lat_mom=radial to switch 

| row | t | KE_int F300 | C300 | BKflat | BKrigid | F300/C300-1 | BKflat/C300-1 | BKrigid/C300-1 | BKrigid/F300-1 |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 4.83873577e+06 | 1.32115307e+42 | 1.32155524e+42 | 1.32115307e+42 | 1.32155524e+42 | -3.043e-04 | -3.043e-04 | 2.331e-14 | 3.044e-04 |
| 1 | 4.83890520e+06 | 1.32115648e+42 | 1.32159199e+42 | 1.32115648e+42 | 1.32159199e+42 | -3.295e-04 | -3.295e-04 | 7.438e-14 | 3.296e-04 |
| 2 | 4.83907463e+06 | 1.32115867e+42 | 1.32162781e+42 | 1.32115867e+42 | 1.32162781e+42 | -3.550e-04 | -3.550e-04 | -1.217e-13 | 3.551e-04 |
| 3 | 4.83924406e+06 | 1.32116015e+42 | 1.32166265e+42 | 1.32116015e+42 | 1.32166265e+42 | -3.802e-04 | -3.802e-04 | -2.542e-14 | 3.803e-04 |
| 4 | 4.83941349e+06 | 1.32116126e+42 | 1.32169654e+42 | 1.32116126e+42 | 1.32169654e+42 | -4.050e-04 | -4.050e-04 | 1.821e-14 | 4.052e-04 |
| 5 | 4.83958292e+06 | 1.32116210e+42 | 1.32172948e+42 | 1.32116210e+42 | 1.32172948e+42 | -4.293e-04 | -4.293e-04 | 8.127e-14 | 4.295e-04 |
| 6 | 4.83975235e+06 | 1.32116278e+42 | 1.32176147e+42 | 1.32116278e+42 | 1.32176147e+42 | -4.529e-04 | -4.529e-04 | 1.494e-13 | 4.532e-04 |
| 7 | 4.83992177e+06 | 1.32116337e+42 | 1.32179258e+42 | 1.32116337e+42 | 1.32179258e+42 | -4.760e-04 | -4.760e-04 | 1.130e-13 | 4.763e-04 |
| 8 | 4.84009102e+06 | 1.32116384e+42 | 1.32182277e+42 | 1.32116384e+42 | 1.32182277e+42 | -4.985e-04 | -4.985e-04 | 1.599e-13 | 4.988e-04 |
| 9 | 4.84026010e+06 | 1.32116412e+42 | 1.32185202e+42 | 1.32116412e+42 | 1.32185202e+42 | -5.204e-04 | -5.204e-04 | 1.692e-13 | 5.207e-04 |
| 10 | 4.84042903e+06 | 1.32116417e+42 | 1.32188029e+42 | 1.32116417e+42 | 1.32188029e+42 | -5.417e-04 | -5.417e-04 | 1.090e-13 | 5.420e-04 |
| 11 | 4.84059784e+06 | 1.32116398e+42 | 1.32190760e+42 | 1.32116398e+42 | 1.32190760e+42 | -5.625e-04 | -5.625e-04 | 1.177e-14 | 5.629e-04 |
| 12 | 4.84076654e+06 | 1.32116353e+42 | 1.32193402e+42 | 1.32116353e+42 | 1.32193402e+42 | -5.828e-04 | -5.828e-04 | -1.332e-15 | 5.832e-04 |
| 13 | 4.84093516e+06 | 1.32116280e+42 | 1.32195957e+42 | 1.32116280e+42 | 1.32195957e+42 | -6.027e-04 | -6.027e-04 | -1.010e-14 | 6.031e-04 |
| 14 | 4.84110369e+06 | 1.32116175e+42 | 1.32198428e+42 | 1.32116175e+42 | 1.32198428e+42 | -6.222e-04 | -6.222e-04 | -2.021e-14 | 6.226e-04 |
| 15 | 4.84127210e+06 | 1.32116036e+42 | 1.32200812e+42 | 1.32116036e+42 | 1.32200812e+42 | -6.413e-04 | -6.413e-04 | -4.163e-14 | 6.417e-04 |
| 16 | 4.84144039e+06 | 1.32115859e+42 | 1.32203106e+42 | 1.32115859e+42 | 1.32203106e+42 | -6.599e-04 | -6.599e-04 | -9.215e-15 | 6.604e-04 |
| 17 | 4.84160856e+06 | 1.32115638e+42 | 1.32205307e+42 | 1.32115638e+42 | 1.32205307e+42 | -6.783e-04 | -6.783e-04 | 5.373e-14 | 6.787e-04 |
| 18 | 4.84177663e+06 | 1.32115368e+42 | 1.32207411e+42 | 1.32115368e+42 | 1.32207411e+42 | -6.962e-04 | -6.962e-04 | 1.865e-13 | 6.967e-04 |
| 19 | 4.84194453e+06 | 1.32115047e+42 | 1.32209422e+42 | 1.32115047e+42 | 1.32209422e+42 | -7.138e-04 | -7.138e-04 | 2.671e-13 | 7.143e-04 |
- BKflat dt: max |/F300 - 1| 0.000e+00, max |/C300 - 1| 8.913e-08
- BKflat L_top: max |/F300 - 1| 0.000e+00, max |/C300 - 1| 2.317e-08
- BKflat L_bot: max |/F300 - 1| 0.000e+00, max |/C300 - 1| 7.390e-07
- BKflat M_tot: max |/F300 - 1| 0.000e+00, max |/C300 - 1| 6.439e-15
- BKrigid dt: max |/F300 - 1| 8.970e-08, max |/C300 - 1| 7.700e-10
- BKrigid L_top: max |/F300 - 1| 2.109e-08, max |/C300 - 1| 2.075e-09
- BKrigid L_bot: max |/F300 - 1| 7.340e-07, max |/C300 - 1| 1.530e-08
- BKrigid M_tot: max |/F300 - 1| 7.883e-15, max |/C300 - 1| 3.997e-15
