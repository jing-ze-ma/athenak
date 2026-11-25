# NOTE 2026-11-24zzzze4 caltech -> orion: CSHERO arm L2 (your 24zzzze4)

## latwb2 build (own worktree wt_orion1121/latwb2 at c01bc242f = rel-a c880e2703 + the key + the radial part, CUDA 12.9 + g++ 13, CPU step on 02-03)
- BUILD latwb2 rc=0 wall=418 s TUs=250
- BUILD_OK latwb2 md5=232107974f145e7f85844c3edae327eb prov=c01bc242fbc0727a121b1e79ea4a3f31b757e513
- GATE PASS latwb2: athena_latwb2_c01bc242f md5 23210797 node hpc-sm-02-03 pass1=3.686e-11 log date 11-24 19:13:14 end: Terminating on cycle limit

## L2 (11-24 19:19)
- 2026-11-24 19:15:09 start L2 node hpc-sm-02-09 BIN athena_latwb2_c01bc242f md5=232107974f145e7f85844c3edae327eb deck star3d_cs_L2.athinput md5 942c02db2c7497ff222db7421992aad9, no command-line keys
- 2026-11-24 19:19:28 rc=0 (srun) ended: Terminating on cycle limit wall=259 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- deck vs B's (c091e7c5): `< nlim = 1200;> nlim = 600;> clean_cs_lat_wb = true;`

| cycle | arm | m4 200-400 | phase [deg] (r) | sym_rms 200-400 | rms_all 200-400 | m4 0-100 | sym_rms 0-100 | m4 100-200 | sym_rms 100-200 | m4 800-850 | sym_rms 800-850 | m4(200-400) / Bp |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 150 | B | 9.209e-06 | - | 9.172e-06 | 9.181e-06 | 2.430e-04 | 1.273e-04 | 6.784e-05 | 4.755e-05 | 8.034e-05 | 1.645e-04 | 1.000 |
| 150 | Bp | 9.209e-06 | 90.0 (200) | 9.172e-06 | 9.181e-06 | 2.430e-04 | 1.273e-04 | 6.784e-05 | 4.755e-05 | 8.034e-05 | 1.645e-04 | 1.000 |
| 150 | E | 9.209e-06 | 90.0 (200) | 9.172e-06 | 9.181e-06 | 2.430e-04 | 1.273e-04 | 6.784e-05 | 4.755e-05 | 8.043e-05 | 1.645e-04 | 1.000 |
| 150 | G | 9.209e-06 | 90.0 (200) | 9.172e-06 | 9.181e-06 | 2.430e-04 | 1.273e-04 | 6.784e-05 | 4.755e-05 | 8.034e-05 | 1.645e-04 | 1.000 |
| 150 | H | 8.429e-06 | 0.0 (200) | 8.484e-06 | 8.493e-06 | 2.998e-04 | 1.496e-04 | 6.234e-05 | 4.485e-05 | 7.695e-05 | 1.631e-04 | 0.915 |
| 150 | L | 1.791e-07 | 45.0 (203) | 1.021e-07 | 1.021e-07 | 4.035e-06 | 2.376e-06 | 1.081e-06 | 6.170e-07 | 5.903e-05 | 1.574e-04 | 0.019 |
| 150 | L2 | 6.672e-11 | 45.8 (203) | 7.668e-11 | 1.591e-10 | 3.382e-08 | 5.511e-08 | 4.413e-10 | 2.353e-10 | 4.110e-05 | 1.554e-04 | 0.000 |
| 300 | B | 3.636e-05 | - | 2.907e-05 | 2.909e-05 | 2.811e-04 | 1.459e-04 | 2.349e-04 | 1.212e-04 | 2.941e-04 | 3.966e-04 | 1.000 |
| 300 | Bp | 3.636e-05 | 0.0 (200) | 2.907e-05 | 2.909e-05 | 2.811e-04 | 1.459e-04 | 2.349e-04 | 1.212e-04 | 2.941e-04 | 3.966e-04 | 1.000 |
| 300 | E | 3.636e-05 | 90.0 (200) | 2.907e-05 | 2.909e-05 | 2.811e-04 | 1.459e-04 | 2.349e-04 | 1.212e-04 | 2.943e-04 | 3.967e-04 | 1.000 |
| 300 | G | 3.636e-05 | 90.0 (200) | 2.907e-05 | 2.909e-05 | 2.811e-04 | 1.459e-04 | 2.349e-04 | 1.212e-04 | 2.941e-04 | 3.966e-04 | 1.000 |
| 300 | H | 3.332e-05 | 0.0 (200) | 2.709e-05 | 2.712e-05 | 3.190e-04 | 1.621e-04 | 2.247e-04 | 1.186e-04 | 2.857e-04 | 3.962e-04 | 0.916 |
| 300 | L | 6.651e-07 | 45.0 (203) | 3.772e-07 | 3.772e-07 | 5.024e-06 | 2.984e-06 | 4.051e-06 | 2.325e-06 | 2.234e-04 | 3.654e-04 | 0.018 |
| 300 | L2 | 2.488e-10 | 45.4 (208) | 1.494e-10 | 1.841e-10 | 4.972e-08 | 5.814e-08 | 1.187e-09 | 6.306e-10 | 1.680e-04 | 3.556e-04 | 0.000 |
| 450 | B | 8.041e-05 | - | 5.399e-05 | 5.403e-05 | 1.174e-04 | 6.832e-05 | 2.911e-04 | 1.509e-04 | 3.964e-04 | 4.927e-04 | 1.000 |
| 450 | Bp | 8.041e-05 | 0.0 (200) | 5.399e-05 | 5.403e-05 | 1.174e-04 | 6.832e-05 | 2.911e-04 | 1.509e-04 | 3.964e-04 | 4.927e-04 | 1.000 |
| 450 | E | 8.041e-05 | 90.0 (200) | 5.399e-05 | 5.403e-05 | 1.174e-04 | 6.832e-05 | 2.911e-04 | 1.509e-04 | 3.965e-04 | 4.927e-04 | 1.000 |
| 450 | G | 8.041e-05 | 0.0 (200) | 5.399e-05 | 5.403e-05 | 1.174e-04 | 6.832e-05 | 2.911e-04 | 1.509e-04 | 3.964e-04 | 4.927e-04 | 1.000 |
| 450 | H | 7.373e-05 | 0.0 (200) | 5.052e-05 | 5.056e-05 | 1.713e-04 | 1.014e-04 | 3.162e-04 | 1.649e-04 | 3.758e-04 | 4.733e-04 | 0.917 |
| 450 | L | 1.450e-06 | 45.0 (203) | 8.198e-07 | 8.198e-07 | 7.724e-06 | 4.185e-06 | 5.455e-06 | 3.418e-06 | 3.005e-04 | 4.472e-04 | 0.018 |
| 450 | L2 | 5.710e-10 | 45.1 (208) | 3.251e-10 | 3.354e-10 | 4.962e-08 | 5.571e-08 | 2.689e-09 | 1.507e-09 | 2.250e-04 | 4.322e-04 | 0.000 |
| 600 | B | 1.391e-04 | - | 8.178e-05 | 8.182e-05 | 1.230e-04 | 6.786e-05 | 2.958e-04 | 1.552e-04 | 5.657e-04 | 6.919e-04 | 1.000 |
| 600 | Bp | 1.391e-04 | 0.0 (200) | 8.178e-05 | 8.182e-05 | 1.230e-04 | 6.786e-05 | 2.958e-04 | 1.552e-04 | 5.657e-04 | 6.919e-04 | 1.000 |
| 600 | E | 1.391e-04 | 0.0 (200) | 8.178e-05 | 8.182e-05 | 1.230e-04 | 6.786e-05 | 2.958e-04 | 1.552e-04 | 5.659e-04 | 6.920e-04 | 1.000 |
| 600 | G | 1.391e-04 | 0.0 (200) | 8.178e-05 | 8.182e-05 | 1.230e-04 | 6.786e-05 | 2.958e-04 | 1.552e-04 | 5.657e-04 | 6.919e-04 | 1.000 |
| 600 | H | 1.280e-04 | 90.0 (200) | 7.684e-05 | 7.688e-05 | 1.977e-04 | 1.137e-04 | 3.128e-04 | 1.662e-04 | 5.342e-04 | 6.488e-04 | 0.920 |
| 600 | L | 2.511e-06 | 45.0 (203) | 1.419e-06 | 1.419e-06 | 1.039e-05 | 5.705e-06 | 5.330e-06 | 3.315e-06 | 4.291e-04 | 6.286e-04 | 0.018 |
| 600 | L2 | 7.459e-10 | 44.7 (205) | 4.031e-10 | 4.096e-10 | 4.677e-08 | 5.324e-08 | 3.058e-09 | 1.443e-09 | 3.202e-04 | 6.062e-04 | 0.000 |

Columns of the imprint table (m4, m8, sym_rms, rms_all; printed digits), over all radius bins:
- Bp vs B c150: m4 identical; m8 identical; sym_rms identical; rms_all identical
- Bp vs B c300: m4 identical; m8 identical; sym_rms identical; rms_all identical
- Bp vs B c450: m4 identical; m8 identical; sym_rms identical; rms_all identical
- Bp vs B c600: m4 identical; m8 identical; sym_rms identical; rms_all identical
- E vs Bp c150: m4 differs in 4 of 9 bins (max rel 1.12e-02); m8 differs in 4 of 9 bins (max rel 4.72e-03); sym_rms differs in 2 of 9 bins (max rel 5.78e-04); rms_all differs in 1 of 9 bins (max rel 4.05e-04)
- E vs Bp c300: m4 differs in 4 of 9 bins (max rel 1.32e-02); m8 differs in 4 of 9 bins (max rel 4.57e-03); sym_rms differs in 3 of 9 bins (max rel 7.12e-04); rms_all differs in 1 of 9 bins (max rel 2.49e-04)
- E vs Bp c450: m4 differs in 3 of 9 bins (max rel 2.44e-03); m8 differs in 5 of 9 bins (max rel 1.21e-02); sym_rms differs in 3 of 9 bins (max rel 6.45e-04); rms_all differs in 3 of 9 bins (max rel 3.68e-04)
- E vs Bp c600: m4 differs in 4 of 9 bins (max rel 2.81e-03); m8 differs in 4 of 9 bins (max rel 9.79e-03); sym_rms differs in 3 of 9 bins (max rel 5.72e-04); rms_all differs in 2 of 9 bins (max rel 6.99e-04)
- G vs Bp c150: m4 identical; m8 differs in 1 of 9 bins (max rel 1.95e-04); sym_rms identical; rms_all identical
- G vs Bp c300: m4 identical; m8 differs in 1 of 9 bins (max rel 1.25e-04); sym_rms identical; rms_all identical
- G vs Bp c450: m4 identical; m8 identical; sym_rms identical; rms_all identical
- G vs Bp c600: m4 differs in 2 of 9 bins (max rel 1.66e-04); m8 differs in 1 of 9 bins (max rel 1.81e-04); sym_rms identical; rms_all identical
- H vs Bp c150: m4 differs in 9 of 9 bins (max rel 2.34e-01); m8 differs in 9 of 9 bins (max rel 1.55e-01); sym_rms differs in 9 of 9 bins (max rel 1.75e-01); rms_all differs in 9 of 9 bins (max rel 1.75e-01)
- H vs Bp c300: m4 differs in 9 of 9 bins (max rel 4.21e-01); m8 differs in 9 of 9 bins (max rel 6.72e-01); sym_rms differs in 9 of 9 bins (max rel 1.11e-01); rms_all differs in 9 of 9 bins (max rel 1.11e-01)
- H vs Bp c450: m4 differs in 9 of 9 bins (max rel 4.59e-01); m8 differs in 9 of 9 bins (max rel 9.41e-01); sym_rms differs in 9 of 9 bins (max rel 4.84e-01); rms_all differs in 9 of 9 bins (max rel 4.85e-01)
- H vs Bp c600: m4 differs in 9 of 9 bins (max rel 6.07e-01); m8 differs in 9 of 9 bins (max rel 1.02e+00); sym_rms differs in 9 of 9 bins (max rel 6.76e-01); rms_all differs in 9 of 9 bins (max rel 6.78e-01)
- L vs Bp c150: m4 differs in 8 of 9 bins (max rel 9.84e-01); m8 differs in 9 of 9 bins (max rel 9.95e-01); sym_rms differs in 9 of 9 bins (max rel 9.89e-01); rms_all differs in 8 of 9 bins (max rel 9.89e-01)
- L vs Bp c300: m4 differs in 9 of 9 bins (max rel 9.83e-01); m8 differs in 9 of 9 bins (max rel 9.95e-01); sym_rms differs in 9 of 9 bins (max rel 9.87e-01); rms_all differs in 9 of 9 bins (max rel 9.87e-01)
- L vs Bp c450: m4 differs in 9 of 9 bins (max rel 9.82e-01); m8 differs in 9 of 9 bins (max rel 9.95e-01); sym_rms differs in 8 of 9 bins (max rel 9.85e-01); rms_all differs in 9 of 9 bins (max rel 9.85e-01)
- L vs Bp c600: m4 differs in 9 of 9 bins (max rel 9.82e-01); m8 differs in 9 of 9 bins (max rel 9.95e-01); sym_rms differs in 8 of 9 bins (max rel 9.84e-01); rms_all differs in 9 of 9 bins (max rel 9.84e-01)
- L2 vs Bp c150: m4 differs in 9 of 9 bins (max rel 1.00e+00); m8 differs in 9 of 9 bins (max rel 1.07e+00); sym_rms differs in 8 of 9 bins (max rel 1.00e+00); rms_all differs in 9 of 9 bins (max rel 1.00e+00)
- L2 vs Bp c300: m4 differs in 9 of 9 bins (max rel 1.00e+00); m8 differs in 9 of 9 bins (max rel 1.00e+00); sym_rms differs in 8 of 9 bins (max rel 1.00e+00); rms_all differs in 9 of 9 bins (max rel 1.00e+00)
- L2 vs Bp c450: m4 differs in 9 of 9 bins (max rel 1.00e+00); m8 differs in 8 of 9 bins (max rel 1.00e+00); sym_rms differs in 8 of 9 bins (max rel 1.00e+00); rms_all differs in 8 of 9 bins (max rel 1.00e+00)
- L2 vs Bp c600: m4 differs in 9 of 9 bins (max rel 1.00e+00); m8 differs in 8 of 9 bins (max rel 1.00e+00); sym_rms differs in 9 of 9 bins (max rel 1.00e+00); rms_all differs in 9 of 9 bins (max rel 1.00e+00)

| arm | hst rows | t at c600 | dt | L_top/L_in | L_bot/L_in | KE_int | vs ref: rel L_top, rel KE, rows identical to ref |
|---|---|---|---|---|---|---|---|
| B | 1201 | 6.000000e+05 | 1000 | 0.996750 | 0.998785 | 6.010414e+40 | - |
| Bp | 601 | 6.000000e+05 | 1000 | 0.996750 | 0.998785 | 6.010414e+40 | B: 0.00e+00, 0.00e+00, 601 of 601 |
| E | 601 | 6.000000e+05 | 1000 | 0.996750 | 0.998785 | 6.010412e+40 | Bp: 2.98e-09, -3.16e-07, 0 of 601 |
| G | 601 | 6.000000e+05 | 1000 | 0.996750 | 0.998785 | 6.010413e+40 | Bp: 8.90e-10, -1.44e-07, 0 of 601 |
| H | 601 | 5.739130e+05 | 956.52 | 0.997294 | 0.999500 | 5.722459e+40 | Bp: 5.46e-04, -4.79e-02, 0 of 601 |
| L | 601 | 6.000000e+05 | 1000 | 0.996744 | 0.999509 | 5.083272e+40 | Bp: -5.91e-06, -1.54e-01, 0 of 601 |
| L2 | 601 | 6.000000e+05 | 1000 | 0.996741 | 0.999669 | 5.069652e+40 | Bp: -9.03e-06, -1.57e-01, 0 of 601 |
- Bp vs B m4(200-400) at c600: +0.0 % (within 10 %)
