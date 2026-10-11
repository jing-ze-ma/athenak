| task | jobs | GPU-h | CPU-h | result |
|---|---|---|---|---|
| TASK-2026-10-10-caltech-plm-gates-rerun | 4316035 (build), 4316036 (gates) | 0 | 3.7 | 12/12 rc 0; sq matches viper check; G1 cfl 1e4 NON-CONVERGED 32-33/51 |
| TASK-2026-10-11-caltech-cpu-agcarA-face | 4333342 (build), 4333462/3 (timing), 4333486/7 (cell/face 40 cyc), 4333595 (ana) | 0 | 3.4 | 2/2 rc 0, 0 NON-CONV, Picard 4.0/9, 15.3 s/cycle both; face vs cell bitwise below 0.8 R_ph and L_bot; max 9e-4 (E, f_rr) at 0.977 R_ph; L_top <= 5e-6 |
| TASK-2026-10-11-caltech-crooked-pipe | 4336052 (build), 4336114 (4 runs, 4x16 thr) | 0 | 26 | 4/4 rc 0, 0 FATAL, 0 NON-CONV, Picard 3.0-3.6/max 7; hesdirk2 NOT ADMISSIBLE 30/8/29 (c220/c2200/c22); max\|vx\| 0; t=10 T_end 0.938/0.855/0.773/0.688/0.554; c2200 vs c220 max rel 0.30/7.6e-3/4.0e-3 at t=1/3/10; c22 vs c220 6.4e-2/1.5e-2 at t=0.5/1 |
