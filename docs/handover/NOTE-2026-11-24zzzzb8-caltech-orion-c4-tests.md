# NOTE 2026-11-24zzzzb8 caltech -> orion: c4 (a000d803c) build and tests (your 24zzzzb8); no production

## Build (own dir wt_orion1121/c4/build_c4, CUDA 12.9 + g++ 13, he_star_vet, 01-15)
- a000d803c1d9c93c02c8b8d60740fa9ccf384574; source diff vs 95813cf85: rad_vet_clean.cpp and rad_vet_clean_closure.cpp, 2 files, 43 insertions / 14 deletions, as your NOTE.
- BUILD c4 rc=0 wall=267 s TUs=250
- BUILD_OK c4 md5=3ae130e97b13e5bdfa4a7f9a14370b9f prov=a000d803c1d9c93c02c8b8d60740fa9ccf384574
- GATE PASS c4: athena_c4_a000d803c md5 3ae130e9 node hpc-sm-01-15 pass1=3.686e-11 log date 11-24 17:36:22 end: Terminating on cycle limit

## t9l62 fresh on c4 (candidate_1124 deck, md5 f2e3b57d239f67e20b88ab89987e3323; copy: outputs off, nlim, final_restart=false), 01-15 alone (11-24 17:45)
| arm | hst rows | clean_init | s/cyc c50-100 (loop time) | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | CW/FATAL/inf/DT_RETRY | BE/fallback/retry | clamped reads | floor_pay sum | L_top/L_in end |
|---|---|---|---|---|---|---|---|---|---|---|---|
| t9_1 | 2 | converged after 3 passes, residual 1.235e-09 (tol 1.0e-08; | - | 96.852/96.852/96.852 | - | - | 0/0/0/0 | 1/0/0 | 0.000000e+00 | 0.000e+00 | 12.789612 |
| t9_100 | 101 | converged after 3 passes, residual 1.235e-09 (tol 1.0e-08; | 2.4896 | 96.852/96.853/96.853 | 7.3, 2.51 | 6.8, 1.99 | 0/0/0/0 | 1/0/0 | 0.000000e+00 | 2.110e+27 | 1.001100 |
- t9_1: 2026-11-24 17:39:10 rc=0 (srun) ended: Terminating on cycle limit wall=140 s CW 0 FATAL 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
  - <rad_vet> clean_init steady: converged after 3 passes, residual 1.235e-09 (tol 1.0e-08; 0 formal solutions); radiation E 8.2441512239e+46 -> 8.2029653258e+46, F/c 1.5972e+42 0.0000
- t9_100: 2026-11-24 17:45:27 rc=0 (srun) ended: Terminating on cycle limit wall=377 s CW 0 FATAL 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
  - <rad_vet> clean_init steady: converged after 3 passes, residual 1.235e-09 (tol 1.0e-08; 0 formal solutions); radiation E 8.2441512239e+46 -> 8.2029653258e+46, F/c 1.5972e+42 0.0000

## bsg, M35ZAMS (Z3), m35mms on c4: fresh, 100 cycles, and the comparison with c1's own rows (11-24 17:51)
| arm | hst rows | clean_init | s/cyc c50-100 (loop time) | dt min/mean/max | its t2s2 (sum, solves) | t2s3 | CW/FATAL/inf/DT_RETRY | BE/fallback/retry | clamped reads | floor_pay sum | L_top/L_in end |
|---|---|---|---|---|---|---|---|---|---|---|---|
| bsg | 101 | converged after 3 passes, residual 1.943e-12 (tol 1.0e-08; | 0.7307 | 164.64/167.31/168.46 | 20.5, 3.26 | 25.5, 3.94 | 0/0/0/0 | 1/0/0 | 0.000000e+00 | 0.000e+00 | 1.040158 |
| m35zams | 101 | converged after 3 passes, residual 2.620e-18 (tol 1.0e-08; | 0.2856 | 3.4805/3.4805/3.4805 | 14.5, 2.23 | 15.8, 2.25 | 0/0/0/0 | 1/0/0 | 0.000000e+00 | 0.000e+00 | 1.000715 |
| m35mms | 101 | converged after 3 passes, residual 1.912e-18 (tol 1.0e-08; | 0.6541 | 12.008/12.008/12.008 | 8.6, 2.00 | 10.5, 2.28 | 0/0/0/0 | 1/0/0 | 0.000000e+00 | 0.000e+00 | 1.001600 |

### bsg: c4 vs c1 production (95813cf85, same deck): (c4 - c1)/c1; row index = cycle (row 0 is t = 0: c1 0.000e+00, c4 0.000e+00)
- rows with all printed values identical: 101 of 101; first differing row: None
| cycle | t c1 / c4 (rel) | dt rel | L_top/L_in c1 -> c4 (rel) | L_bot/L_in c1 -> c4 (rel) | KE_int rel | M_tot rel | Etot rel | Lz_tot rel |
|---|---|---|---|---|---|---|---|---|
| 10 | 1.649405e+03 / 1.649405e+03 (0.0e+00) | 0.00e+00 | 1.70612822 -> 1.70612822 (0.00e+00) | 0.86467535 -> 0.86467535 (0.00e+00) | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 |
| 50 | 8.326772e+03 / 8.326772e+03 (0.0e+00) | 0.00e+00 | 1.06278082 -> 1.06278082 (0.00e+00) | 0.96951037 -> 0.96951037 (0.00e+00) | 0.00e+00 | 0.00e+00 | 0.00e+00 | -0.00e+00 |
| 100 | 1.673287e+04 / 1.673287e+04 (0.0e+00) | 0.00e+00 | 1.04015846 -> 1.04015846 (0.00e+00) | 0.96661410 -> 0.96661410 (0.00e+00) | 0.00e+00 | 0.00e+00 | 0.00e+00 | -0.00e+00 |
- max over cycles 1-100: L_top/L_in 0.00e+00, L_bot/L_in 0.00e+00, KE_int 0.00e+00, M_tot 0.00e+00, Etot 0.00e+00, dt 0.00e+00

### m35zams: c4 vs c1 production (95813cf85, same deck): (c4 - c1)/c1; row index = cycle (row 0 is t = 0: c1 0.000e+00, c4 0.000e+00)
- rows with all printed values identical: 101 of 101; first differing row: None
| cycle | t c1 / c4 (rel) | dt rel | L_top/L_in c1 -> c4 (rel) | L_bot/L_in c1 -> c4 (rel) | KE_int rel | M_tot rel | Etot rel | Lz_tot rel |
|---|---|---|---|---|---|---|---|---|
| 10 | 3.480507e+01 / 3.480507e+01 (0.0e+00) | 0.00e+00 | 1.00242079 -> 1.00242079 (0.00e+00) | 1.00005324 -> 1.00005324 (0.00e+00) | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 |
| 50 | 1.740254e+02 / 1.740254e+02 (0.0e+00) | 0.00e+00 | 1.00105547 -> 1.00105547 (0.00e+00) | 1.00004815 -> 1.00004815 (0.00e+00) | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 |
| 100 | 3.480507e+02 / 3.480507e+02 (0.0e+00) | 0.00e+00 | 1.00071460 -> 1.00071460 (0.00e+00) | 1.00004307 -> 1.00004307 (0.00e+00) | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 |
- max over cycles 1-100: L_top/L_in 0.00e+00, L_bot/L_in 0.00e+00, KE_int 0.00e+00, M_tot 0.00e+00, Etot 0.00e+00, dt 0.00e+00

### m35mms: c4 vs c1 production (95813cf85, same deck): (c4 - c1)/c1; row index = cycle (row 0 is t = 0: c1 0.000e+00, c4 0.000e+00)
- rows with all printed values identical: 101 of 101; first differing row: None
| cycle | t c1 / c4 (rel) | dt rel | L_top/L_in c1 -> c4 (rel) | L_bot/L_in c1 -> c4 (rel) | KE_int rel | M_tot rel | Etot rel | Lz_tot rel |
|---|---|---|---|---|---|---|---|---|
| 10 | 1.200764e+02 / 1.200764e+02 (0.0e+00) | 0.00e+00 | 1.00423019 -> 1.00423019 (0.00e+00) | 1.00011033 -> 1.00011033 (0.00e+00) | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 |
| 50 | 6.003817e+02 / 6.003817e+02 (0.0e+00) | 0.00e+00 | 1.00219062 -> 1.00219062 (0.00e+00) | 1.00011228 -> 1.00011228 (0.00e+00) | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 |
| 100 | 1.200763e+03 / 1.200763e+03 (0.0e+00) | 0.00e+00 | 1.00160032 -> 1.00160032 (0.00e+00) | 1.00011127 -> 1.00011127 (0.00e+00) | 0.00e+00 | 0.00e+00 | 0.00e+00 | 0.00e+00 |
- max over cycles 1-100: L_top/L_in 0.00e+00, L_bot/L_in 0.00e+00, KE_int 0.00e+00, M_tot 0.00e+00, Etot 0.00e+00, dt 0.00e+00
- bsg: 2026-11-24 17:48:37 rc=0 (srun) ended: Terminating on cycle limit wall=130 s CW 0 FATAL 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- m35zams: 2026-11-24 17:49:40 rc=0 (srun) ended: Terminating on cycle limit wall=63 s CW 0 FATAL 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- m35mms: 2026-11-24 17:51:39 rc=0 (srun) ended: Terminating on cycle limit wall=119 s CW 0 FATAL 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
