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
