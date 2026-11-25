# NOTE 2026-11-24zzzzc6 caltech -> orion: CSHERO Bp, E, G, H (your 24zzzzc6 item 1)

- B's deck (c091e7c5): `clean_stag_thinf = false`; merge_lat IS on: `merge_lat_rc = 2.0629e12:3.4307e12:5.7211e12, merge_lat_f = 8:4:2, merge_lat_rad = false, merge_lat_diag = 1000, merge_lat_blockend = exchange`. So H is run: those five lines removed (the code turns the merge off when merge_lat_rc is absent).
- Arms: Bp = B's deck (nlim 600) on 53290b5af; E = thinf true on 87af1d65c; G = thinf true on 53290b5af (same deck text as E); H = Bp without the merge. 4 ranks, 4 H200, node 02-03. B itself ran with T2 (7d7c4a71). B's c150 / 300 / 450 rows come from its old table (no phase); its c600 is also in the old table.
- Source check on arm D, carried over from my zzzzb6 correction: cs_vertex_3cell IS read (mesh.cpp) and two sites of the hydro FOFC trial state (hydro_fofc.cpp, fofc_etotgrav.cpp) do run in this deck; FOFC never fired in B or D, every other use is behind conduction / MHD, none is in rad_vet. B and D data payloads at c0 and c600 are byte-identical (double). So inert in effect; it is not literally unread by the hydro update.
- Stop rule: after E, m4(200-400) at c600 below Bp / 3 -> G and H are not run; the c600 dumps of all arms are kept until you score.

## Bp (11-24 18:13)
- 2026-11-24 18:09:23 start Bp node hpc-sm-02-03 BIN athena_sdosdo_53290b5af md5=573465bddd359381dd270fd0428e8683 deck star3d_cs_Bp.athinput md5 d3de281af1cda9971e248d1fb1efd17f, no command-line keys
- 2026-11-24 18:13:40 rc=0 (srun) ended: Terminating on cycle limit wall=257 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- deck vs B's (c091e7c5): `< nlim = 1200;> nlim = 600;`

| cycle | arm | m4 200-400 | phase [deg] (r) | sym_rms 200-400 | rms_all 200-400 | m4 0-100 | sym_rms 0-100 | m4 100-200 | sym_rms 100-200 | m4 800-850 | sym_rms 800-850 | m4(200-400) / Bp |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 150 | B | 9.209e-06 | - | 9.172e-06 | 9.181e-06 | 2.430e-04 | 1.273e-04 | 6.784e-05 | 4.755e-05 | 8.034e-05 | 1.645e-04 | 1.000 |
| 150 | Bp | 9.209e-06 | 90.0 (200) | 9.172e-06 | 9.181e-06 | 2.430e-04 | 1.273e-04 | 6.784e-05 | 4.755e-05 | 8.034e-05 | 1.645e-04 | 1.000 |
| 300 | B | 3.636e-05 | - | 2.907e-05 | 2.909e-05 | 2.811e-04 | 1.459e-04 | 2.349e-04 | 1.212e-04 | 2.941e-04 | 3.966e-04 | 1.000 |
| 300 | Bp | 3.636e-05 | 0.0 (200) | 2.907e-05 | 2.909e-05 | 2.811e-04 | 1.459e-04 | 2.349e-04 | 1.212e-04 | 2.941e-04 | 3.966e-04 | 1.000 |
| 450 | B | 8.041e-05 | - | 5.399e-05 | 5.403e-05 | 1.174e-04 | 6.832e-05 | 2.911e-04 | 1.509e-04 | 3.964e-04 | 4.927e-04 | 1.000 |
| 450 | Bp | 8.041e-05 | 0.0 (200) | 5.399e-05 | 5.403e-05 | 1.174e-04 | 6.832e-05 | 2.911e-04 | 1.509e-04 | 3.964e-04 | 4.927e-04 | 1.000 |
| 600 | B | 1.391e-04 | - | 8.178e-05 | 8.182e-05 | 1.230e-04 | 6.786e-05 | 2.958e-04 | 1.552e-04 | 5.657e-04 | 6.919e-04 | 1.000 |
| 600 | Bp | 1.391e-04 | 0.0 (200) | 8.178e-05 | 8.182e-05 | 1.230e-04 | 6.786e-05 | 2.958e-04 | 1.552e-04 | 5.657e-04 | 6.919e-04 | 1.000 |

Columns of the imprint table (m4, m8, sym_rms, rms_all; printed digits), over all radius bins:
- Bp vs B c150: m4 identical; m8 identical; sym_rms identical; rms_all identical
- Bp vs B c300: m4 identical; m8 identical; sym_rms identical; rms_all identical
- Bp vs B c450: m4 identical; m8 identical; sym_rms identical; rms_all identical
- Bp vs B c600: m4 identical; m8 identical; sym_rms identical; rms_all identical

| arm | hst rows | t at c600 | dt | L_top/L_in | L_bot/L_in | KE_int | vs ref: rel L_top, rel KE, rows identical to ref |
|---|---|---|---|---|---|---|---|
| B | 1201 | 6.000000e+05 | 1000 | 0.996750 | 0.998785 | 6.010414e+40 | - |
| Bp | 601 | 6.000000e+05 | 1000 | 0.996750 | 0.998785 | 6.010414e+40 | B: 0.00e+00, 0.00e+00, 601 of 601 |
- Bp vs B m4(200-400) at c600: +0.0 % (within 10 %)

## E (11-24 18:18)
- 2026-11-24 18:14:11 start E node hpc-sm-02-03 BIN athena_thinfcs_87af1d65c md5=fb0549d6630c36c5e7405ef6e529fec3 deck star3d_cs_E.athinput md5 93c07ab2edf575500afcd1e0556ff1e2, no command-line keys
- 2026-11-24 18:18:29 rc=0 (srun) ended: Terminating on cycle limit wall=258 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- deck vs B's (c091e7c5): `< nlim = 1200;> nlim = 600;< clean_stag_thinf = false;> clean_stag_thinf = true;`

| cycle | arm | m4 200-400 | phase [deg] (r) | sym_rms 200-400 | rms_all 200-400 | m4 0-100 | sym_rms 0-100 | m4 100-200 | sym_rms 100-200 | m4 800-850 | sym_rms 800-850 | m4(200-400) / Bp |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 150 | B | 9.209e-06 | - | 9.172e-06 | 9.181e-06 | 2.430e-04 | 1.273e-04 | 6.784e-05 | 4.755e-05 | 8.034e-05 | 1.645e-04 | 1.000 |
| 150 | Bp | 9.209e-06 | 90.0 (200) | 9.172e-06 | 9.181e-06 | 2.430e-04 | 1.273e-04 | 6.784e-05 | 4.755e-05 | 8.034e-05 | 1.645e-04 | 1.000 |
| 150 | E | 9.209e-06 | 90.0 (200) | 9.172e-06 | 9.181e-06 | 2.430e-04 | 1.273e-04 | 6.784e-05 | 4.755e-05 | 8.043e-05 | 1.645e-04 | 1.000 |
| 300 | B | 3.636e-05 | - | 2.907e-05 | 2.909e-05 | 2.811e-04 | 1.459e-04 | 2.349e-04 | 1.212e-04 | 2.941e-04 | 3.966e-04 | 1.000 |
| 300 | Bp | 3.636e-05 | 0.0 (200) | 2.907e-05 | 2.909e-05 | 2.811e-04 | 1.459e-04 | 2.349e-04 | 1.212e-04 | 2.941e-04 | 3.966e-04 | 1.000 |
| 300 | E | 3.636e-05 | 90.0 (200) | 2.907e-05 | 2.909e-05 | 2.811e-04 | 1.459e-04 | 2.349e-04 | 1.212e-04 | 2.943e-04 | 3.967e-04 | 1.000 |
| 450 | B | 8.041e-05 | - | 5.399e-05 | 5.403e-05 | 1.174e-04 | 6.832e-05 | 2.911e-04 | 1.509e-04 | 3.964e-04 | 4.927e-04 | 1.000 |
| 450 | Bp | 8.041e-05 | 0.0 (200) | 5.399e-05 | 5.403e-05 | 1.174e-04 | 6.832e-05 | 2.911e-04 | 1.509e-04 | 3.964e-04 | 4.927e-04 | 1.000 |
| 450 | E | 8.041e-05 | 90.0 (200) | 5.399e-05 | 5.403e-05 | 1.174e-04 | 6.832e-05 | 2.911e-04 | 1.509e-04 | 3.965e-04 | 4.927e-04 | 1.000 |
| 600 | B | 1.391e-04 | - | 8.178e-05 | 8.182e-05 | 1.230e-04 | 6.786e-05 | 2.958e-04 | 1.552e-04 | 5.657e-04 | 6.919e-04 | 1.000 |
| 600 | Bp | 1.391e-04 | 0.0 (200) | 8.178e-05 | 8.182e-05 | 1.230e-04 | 6.786e-05 | 2.958e-04 | 1.552e-04 | 5.657e-04 | 6.919e-04 | 1.000 |
| 600 | E | 1.391e-04 | 0.0 (200) | 8.178e-05 | 8.182e-05 | 1.230e-04 | 6.786e-05 | 2.958e-04 | 1.552e-04 | 5.659e-04 | 6.920e-04 | 1.000 |

Columns of the imprint table (m4, m8, sym_rms, rms_all; printed digits), over all radius bins:
- Bp vs B c150: m4 identical; m8 identical; sym_rms identical; rms_all identical
- Bp vs B c300: m4 identical; m8 identical; sym_rms identical; rms_all identical
- Bp vs B c450: m4 identical; m8 identical; sym_rms identical; rms_all identical
- Bp vs B c600: m4 identical; m8 identical; sym_rms identical; rms_all identical
- E vs Bp c150: m4 differs in 4 of 9 bins (max rel 1.12e-02); m8 differs in 4 of 9 bins (max rel 4.72e-03); sym_rms differs in 2 of 9 bins (max rel 5.78e-04); rms_all differs in 1 of 9 bins (max rel 4.05e-04)
- E vs Bp c300: m4 differs in 4 of 9 bins (max rel 1.32e-02); m8 differs in 4 of 9 bins (max rel 4.57e-03); sym_rms differs in 3 of 9 bins (max rel 7.12e-04); rms_all differs in 1 of 9 bins (max rel 2.49e-04)
- E vs Bp c450: m4 differs in 3 of 9 bins (max rel 2.44e-03); m8 differs in 5 of 9 bins (max rel 1.21e-02); sym_rms differs in 3 of 9 bins (max rel 6.45e-04); rms_all differs in 3 of 9 bins (max rel 3.68e-04)
- E vs Bp c600: m4 differs in 4 of 9 bins (max rel 2.81e-03); m8 differs in 4 of 9 bins (max rel 9.79e-03); sym_rms differs in 3 of 9 bins (max rel 5.72e-04); rms_all differs in 2 of 9 bins (max rel 6.99e-04)

| arm | hst rows | t at c600 | dt | L_top/L_in | L_bot/L_in | KE_int | vs ref: rel L_top, rel KE, rows identical to ref |
|---|---|---|---|---|---|---|---|
| B | 1201 | 6.000000e+05 | 1000 | 0.996750 | 0.998785 | 6.010414e+40 | - |
| Bp | 601 | 6.000000e+05 | 1000 | 0.996750 | 0.998785 | 6.010414e+40 | B: 0.00e+00, 0.00e+00, 601 of 601 |
| E | 601 | 6.000000e+05 | 1000 | 0.996750 | 0.998785 | 6.010412e+40 | Bp: 2.98e-09, -3.16e-07, 0 of 601 |
- Bp vs B m4(200-400) at c600: +0.0 % (within 10 %)

## G (11-24 18:23)
- 2026-11-24 18:18:59 start G node hpc-sm-02-03 BIN athena_sdosdo_53290b5af md5=573465bddd359381dd270fd0428e8683 deck star3d_cs_G.athinput md5 93c07ab2edf575500afcd1e0556ff1e2, no command-line keys
- 2026-11-24 18:23:22 rc=0 (srun) ended: Terminating on cycle limit wall=263 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- deck vs B's (c091e7c5): `< nlim = 1200;> nlim = 600;< clean_stag_thinf = false;> clean_stag_thinf = true;`

| cycle | arm | m4 200-400 | phase [deg] (r) | sym_rms 200-400 | rms_all 200-400 | m4 0-100 | sym_rms 0-100 | m4 100-200 | sym_rms 100-200 | m4 800-850 | sym_rms 800-850 | m4(200-400) / Bp |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 150 | B | 9.209e-06 | - | 9.172e-06 | 9.181e-06 | 2.430e-04 | 1.273e-04 | 6.784e-05 | 4.755e-05 | 8.034e-05 | 1.645e-04 | 1.000 |
| 150 | Bp | 9.209e-06 | 90.0 (200) | 9.172e-06 | 9.181e-06 | 2.430e-04 | 1.273e-04 | 6.784e-05 | 4.755e-05 | 8.034e-05 | 1.645e-04 | 1.000 |
| 150 | E | 9.209e-06 | 90.0 (200) | 9.172e-06 | 9.181e-06 | 2.430e-04 | 1.273e-04 | 6.784e-05 | 4.755e-05 | 8.043e-05 | 1.645e-04 | 1.000 |
| 150 | G | 9.209e-06 | 90.0 (200) | 9.172e-06 | 9.181e-06 | 2.430e-04 | 1.273e-04 | 6.784e-05 | 4.755e-05 | 8.034e-05 | 1.645e-04 | 1.000 |
| 300 | B | 3.636e-05 | - | 2.907e-05 | 2.909e-05 | 2.811e-04 | 1.459e-04 | 2.349e-04 | 1.212e-04 | 2.941e-04 | 3.966e-04 | 1.000 |
| 300 | Bp | 3.636e-05 | 0.0 (200) | 2.907e-05 | 2.909e-05 | 2.811e-04 | 1.459e-04 | 2.349e-04 | 1.212e-04 | 2.941e-04 | 3.966e-04 | 1.000 |
| 300 | E | 3.636e-05 | 90.0 (200) | 2.907e-05 | 2.909e-05 | 2.811e-04 | 1.459e-04 | 2.349e-04 | 1.212e-04 | 2.943e-04 | 3.967e-04 | 1.000 |
| 300 | G | 3.636e-05 | 90.0 (200) | 2.907e-05 | 2.909e-05 | 2.811e-04 | 1.459e-04 | 2.349e-04 | 1.212e-04 | 2.941e-04 | 3.966e-04 | 1.000 |
| 450 | B | 8.041e-05 | - | 5.399e-05 | 5.403e-05 | 1.174e-04 | 6.832e-05 | 2.911e-04 | 1.509e-04 | 3.964e-04 | 4.927e-04 | 1.000 |
| 450 | Bp | 8.041e-05 | 0.0 (200) | 5.399e-05 | 5.403e-05 | 1.174e-04 | 6.832e-05 | 2.911e-04 | 1.509e-04 | 3.964e-04 | 4.927e-04 | 1.000 |
| 450 | E | 8.041e-05 | 90.0 (200) | 5.399e-05 | 5.403e-05 | 1.174e-04 | 6.832e-05 | 2.911e-04 | 1.509e-04 | 3.965e-04 | 4.927e-04 | 1.000 |
| 450 | G | 8.041e-05 | 0.0 (200) | 5.399e-05 | 5.403e-05 | 1.174e-04 | 6.832e-05 | 2.911e-04 | 1.509e-04 | 3.964e-04 | 4.927e-04 | 1.000 |
| 600 | B | 1.391e-04 | - | 8.178e-05 | 8.182e-05 | 1.230e-04 | 6.786e-05 | 2.958e-04 | 1.552e-04 | 5.657e-04 | 6.919e-04 | 1.000 |
| 600 | Bp | 1.391e-04 | 0.0 (200) | 8.178e-05 | 8.182e-05 | 1.230e-04 | 6.786e-05 | 2.958e-04 | 1.552e-04 | 5.657e-04 | 6.919e-04 | 1.000 |
| 600 | E | 1.391e-04 | 0.0 (200) | 8.178e-05 | 8.182e-05 | 1.230e-04 | 6.786e-05 | 2.958e-04 | 1.552e-04 | 5.659e-04 | 6.920e-04 | 1.000 |
| 600 | G | 1.391e-04 | 0.0 (200) | 8.178e-05 | 8.182e-05 | 1.230e-04 | 6.786e-05 | 2.958e-04 | 1.552e-04 | 5.657e-04 | 6.919e-04 | 1.000 |

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

| arm | hst rows | t at c600 | dt | L_top/L_in | L_bot/L_in | KE_int | vs ref: rel L_top, rel KE, rows identical to ref |
|---|---|---|---|---|---|---|---|
| B | 1201 | 6.000000e+05 | 1000 | 0.996750 | 0.998785 | 6.010414e+40 | - |
| Bp | 601 | 6.000000e+05 | 1000 | 0.996750 | 0.998785 | 6.010414e+40 | B: 0.00e+00, 0.00e+00, 601 of 601 |
| E | 601 | 6.000000e+05 | 1000 | 0.996750 | 0.998785 | 6.010412e+40 | Bp: 2.98e-09, -3.16e-07, 0 of 601 |
| G | 601 | 6.000000e+05 | 1000 | 0.996750 | 0.998785 | 6.010413e+40 | Bp: 8.90e-10, -1.44e-07, 0 of 601 |
- Bp vs B m4(200-400) at c600: +0.0 % (within 10 %)

## H (11-24 18:28)
- 2026-11-24 18:23:46 start H node hpc-sm-02-03 BIN athena_sdosdo_53290b5af md5=573465bddd359381dd270fd0428e8683 deck star3d_cs_H.athinput md5 7f4a6ab7ce74556dd603ac7e6eba484d, no command-line keys
- 2026-11-24 18:27:56 rc=0 (srun) ended: Terminating on cycle limit wall=250 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- deck vs B's (c091e7c5): `< nlim = 1200;> nlim = 600;< merge_lat_rc = 2.0629e12:3.4307e12:5.7211e12;< merge_lat_f = 8:4:2;< merge_lat_rad = false;< merge_lat_diag = 1000;< merge_lat_blockend = exchange;`

| cycle | arm | m4 200-400 | phase [deg] (r) | sym_rms 200-400 | rms_all 200-400 | m4 0-100 | sym_rms 0-100 | m4 100-200 | sym_rms 100-200 | m4 800-850 | sym_rms 800-850 | m4(200-400) / Bp |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 150 | B | 9.209e-06 | - | 9.172e-06 | 9.181e-06 | 2.430e-04 | 1.273e-04 | 6.784e-05 | 4.755e-05 | 8.034e-05 | 1.645e-04 | 1.000 |
| 150 | Bp | 9.209e-06 | 90.0 (200) | 9.172e-06 | 9.181e-06 | 2.430e-04 | 1.273e-04 | 6.784e-05 | 4.755e-05 | 8.034e-05 | 1.645e-04 | 1.000 |
| 150 | E | 9.209e-06 | 90.0 (200) | 9.172e-06 | 9.181e-06 | 2.430e-04 | 1.273e-04 | 6.784e-05 | 4.755e-05 | 8.043e-05 | 1.645e-04 | 1.000 |
| 150 | G | 9.209e-06 | 90.0 (200) | 9.172e-06 | 9.181e-06 | 2.430e-04 | 1.273e-04 | 6.784e-05 | 4.755e-05 | 8.034e-05 | 1.645e-04 | 1.000 |
| 150 | H | 8.429e-06 | 0.0 (200) | 8.484e-06 | 8.493e-06 | 2.998e-04 | 1.496e-04 | 6.234e-05 | 4.485e-05 | 7.695e-05 | 1.631e-04 | 0.915 |
| 300 | B | 3.636e-05 | - | 2.907e-05 | 2.909e-05 | 2.811e-04 | 1.459e-04 | 2.349e-04 | 1.212e-04 | 2.941e-04 | 3.966e-04 | 1.000 |
| 300 | Bp | 3.636e-05 | 0.0 (200) | 2.907e-05 | 2.909e-05 | 2.811e-04 | 1.459e-04 | 2.349e-04 | 1.212e-04 | 2.941e-04 | 3.966e-04 | 1.000 |
| 300 | E | 3.636e-05 | 90.0 (200) | 2.907e-05 | 2.909e-05 | 2.811e-04 | 1.459e-04 | 2.349e-04 | 1.212e-04 | 2.943e-04 | 3.967e-04 | 1.000 |
| 300 | G | 3.636e-05 | 90.0 (200) | 2.907e-05 | 2.909e-05 | 2.811e-04 | 1.459e-04 | 2.349e-04 | 1.212e-04 | 2.941e-04 | 3.966e-04 | 1.000 |
| 300 | H | 3.332e-05 | 0.0 (200) | 2.709e-05 | 2.712e-05 | 3.190e-04 | 1.621e-04 | 2.247e-04 | 1.186e-04 | 2.857e-04 | 3.962e-04 | 0.916 |
| 450 | B | 8.041e-05 | - | 5.399e-05 | 5.403e-05 | 1.174e-04 | 6.832e-05 | 2.911e-04 | 1.509e-04 | 3.964e-04 | 4.927e-04 | 1.000 |
| 450 | Bp | 8.041e-05 | 0.0 (200) | 5.399e-05 | 5.403e-05 | 1.174e-04 | 6.832e-05 | 2.911e-04 | 1.509e-04 | 3.964e-04 | 4.927e-04 | 1.000 |
| 450 | E | 8.041e-05 | 90.0 (200) | 5.399e-05 | 5.403e-05 | 1.174e-04 | 6.832e-05 | 2.911e-04 | 1.509e-04 | 3.965e-04 | 4.927e-04 | 1.000 |
| 450 | G | 8.041e-05 | 0.0 (200) | 5.399e-05 | 5.403e-05 | 1.174e-04 | 6.832e-05 | 2.911e-04 | 1.509e-04 | 3.964e-04 | 4.927e-04 | 1.000 |
| 450 | H | 7.373e-05 | 0.0 (200) | 5.052e-05 | 5.056e-05 | 1.713e-04 | 1.014e-04 | 3.162e-04 | 1.649e-04 | 3.758e-04 | 4.733e-04 | 0.917 |
| 600 | B | 1.391e-04 | - | 8.178e-05 | 8.182e-05 | 1.230e-04 | 6.786e-05 | 2.958e-04 | 1.552e-04 | 5.657e-04 | 6.919e-04 | 1.000 |
| 600 | Bp | 1.391e-04 | 0.0 (200) | 8.178e-05 | 8.182e-05 | 1.230e-04 | 6.786e-05 | 2.958e-04 | 1.552e-04 | 5.657e-04 | 6.919e-04 | 1.000 |
| 600 | E | 1.391e-04 | 0.0 (200) | 8.178e-05 | 8.182e-05 | 1.230e-04 | 6.786e-05 | 2.958e-04 | 1.552e-04 | 5.659e-04 | 6.920e-04 | 1.000 |
| 600 | G | 1.391e-04 | 0.0 (200) | 8.178e-05 | 8.182e-05 | 1.230e-04 | 6.786e-05 | 2.958e-04 | 1.552e-04 | 5.657e-04 | 6.919e-04 | 1.000 |
| 600 | H | 1.280e-04 | 90.0 (200) | 7.684e-05 | 7.688e-05 | 1.977e-04 | 1.137e-04 | 3.128e-04 | 1.662e-04 | 5.342e-04 | 6.488e-04 | 0.920 |

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

| arm | hst rows | t at c600 | dt | L_top/L_in | L_bot/L_in | KE_int | vs ref: rel L_top, rel KE, rows identical to ref |
|---|---|---|---|---|---|---|---|
| B | 1201 | 6.000000e+05 | 1000 | 0.996750 | 0.998785 | 6.010414e+40 | - |
| Bp | 601 | 6.000000e+05 | 1000 | 0.996750 | 0.998785 | 6.010414e+40 | B: 0.00e+00, 0.00e+00, 601 of 601 |
| E | 601 | 6.000000e+05 | 1000 | 0.996750 | 0.998785 | 6.010412e+40 | Bp: 2.98e-09, -3.16e-07, 0 of 601 |
| G | 601 | 6.000000e+05 | 1000 | 0.996750 | 0.998785 | 6.010413e+40 | Bp: 8.90e-10, -1.44e-07, 0 of 601 |
| H | 601 | 5.739130e+05 | 956.52 | 0.997294 | 0.999500 | 5.722459e+40 | Bp: 5.46e-04, -4.79e-02, 0 of 601 |
- Bp vs B m4(200-400) at c600: +0.0 % (within 10 %)

## latwb build (own worktree wt_orion1121/latwb at 83663792f = 53290b5af + the key, CUDA 12.9 + g++ 13, CPU step on 01-17)
- BUILD latwb rc=0 wall=273 s TUs=250
- BUILD_OK latwb md5=a4c12aa1edc0718af2cb9f6bd10de9b2 prov=83663792f0079aa904d51053b0da1b01bc04741a
- GATE PASS latwb: athena_latwb_83663792f md5 a4c12aa1 node hpc-sm-01-17 pass1=3.686e-11 log date 11-24 18:29:33 end: Terminating on cycle limit

## L (11-24 18:34)
- 2026-11-24 18:30:20 start L node hpc-sm-01-17 BIN athena_latwb_83663792f md5=a4c12aa1edc0718af2cb9f6bd10de9b2 deck star3d_cs_L.athinput md5 942c02db2c7497ff222db7421992aad9, no command-line keys
- 2026-11-24 18:34:31 rc=0 (srun) ended: Terminating on cycle limit wall=251 s CW 0 FATAL 0 inf 0 DT_RETRY 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- deck vs B's (c091e7c5): `< nlim = 1200;> nlim = 600;> clean_cs_lat_wb = true;`

| cycle | arm | m4 200-400 | phase [deg] (r) | sym_rms 200-400 | rms_all 200-400 | m4 0-100 | sym_rms 0-100 | m4 100-200 | sym_rms 100-200 | m4 800-850 | sym_rms 800-850 | m4(200-400) / Bp |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 150 | B | 9.209e-06 | - | 9.172e-06 | 9.181e-06 | 2.430e-04 | 1.273e-04 | 6.784e-05 | 4.755e-05 | 8.034e-05 | 1.645e-04 | 1.000 |
| 150 | Bp | 9.209e-06 | 90.0 (200) | 9.172e-06 | 9.181e-06 | 2.430e-04 | 1.273e-04 | 6.784e-05 | 4.755e-05 | 8.034e-05 | 1.645e-04 | 1.000 |
| 150 | E | 9.209e-06 | 90.0 (200) | 9.172e-06 | 9.181e-06 | 2.430e-04 | 1.273e-04 | 6.784e-05 | 4.755e-05 | 8.043e-05 | 1.645e-04 | 1.000 |
| 150 | G | 9.209e-06 | 90.0 (200) | 9.172e-06 | 9.181e-06 | 2.430e-04 | 1.273e-04 | 6.784e-05 | 4.755e-05 | 8.034e-05 | 1.645e-04 | 1.000 |
| 150 | H | 8.429e-06 | 0.0 (200) | 8.484e-06 | 8.493e-06 | 2.998e-04 | 1.496e-04 | 6.234e-05 | 4.485e-05 | 7.695e-05 | 1.631e-04 | 0.915 |
| 150 | L | 1.791e-07 | 45.0 (203) | 1.021e-07 | 1.021e-07 | 4.035e-06 | 2.376e-06 | 1.081e-06 | 6.170e-07 | 5.903e-05 | 1.574e-04 | 0.019 |
| 300 | B | 3.636e-05 | - | 2.907e-05 | 2.909e-05 | 2.811e-04 | 1.459e-04 | 2.349e-04 | 1.212e-04 | 2.941e-04 | 3.966e-04 | 1.000 |
| 300 | Bp | 3.636e-05 | 0.0 (200) | 2.907e-05 | 2.909e-05 | 2.811e-04 | 1.459e-04 | 2.349e-04 | 1.212e-04 | 2.941e-04 | 3.966e-04 | 1.000 |
| 300 | E | 3.636e-05 | 90.0 (200) | 2.907e-05 | 2.909e-05 | 2.811e-04 | 1.459e-04 | 2.349e-04 | 1.212e-04 | 2.943e-04 | 3.967e-04 | 1.000 |
| 300 | G | 3.636e-05 | 90.0 (200) | 2.907e-05 | 2.909e-05 | 2.811e-04 | 1.459e-04 | 2.349e-04 | 1.212e-04 | 2.941e-04 | 3.966e-04 | 1.000 |
| 300 | H | 3.332e-05 | 0.0 (200) | 2.709e-05 | 2.712e-05 | 3.190e-04 | 1.621e-04 | 2.247e-04 | 1.186e-04 | 2.857e-04 | 3.962e-04 | 0.916 |
| 300 | L | 6.651e-07 | 45.0 (203) | 3.772e-07 | 3.772e-07 | 5.024e-06 | 2.984e-06 | 4.051e-06 | 2.325e-06 | 2.234e-04 | 3.654e-04 | 0.018 |
| 450 | B | 8.041e-05 | - | 5.399e-05 | 5.403e-05 | 1.174e-04 | 6.832e-05 | 2.911e-04 | 1.509e-04 | 3.964e-04 | 4.927e-04 | 1.000 |
| 450 | Bp | 8.041e-05 | 0.0 (200) | 5.399e-05 | 5.403e-05 | 1.174e-04 | 6.832e-05 | 2.911e-04 | 1.509e-04 | 3.964e-04 | 4.927e-04 | 1.000 |
| 450 | E | 8.041e-05 | 90.0 (200) | 5.399e-05 | 5.403e-05 | 1.174e-04 | 6.832e-05 | 2.911e-04 | 1.509e-04 | 3.965e-04 | 4.927e-04 | 1.000 |
| 450 | G | 8.041e-05 | 0.0 (200) | 5.399e-05 | 5.403e-05 | 1.174e-04 | 6.832e-05 | 2.911e-04 | 1.509e-04 | 3.964e-04 | 4.927e-04 | 1.000 |
| 450 | H | 7.373e-05 | 0.0 (200) | 5.052e-05 | 5.056e-05 | 1.713e-04 | 1.014e-04 | 3.162e-04 | 1.649e-04 | 3.758e-04 | 4.733e-04 | 0.917 |
| 450 | L | 1.450e-06 | 45.0 (203) | 8.198e-07 | 8.198e-07 | 7.724e-06 | 4.185e-06 | 5.455e-06 | 3.418e-06 | 3.005e-04 | 4.472e-04 | 0.018 |
| 600 | B | 1.391e-04 | - | 8.178e-05 | 8.182e-05 | 1.230e-04 | 6.786e-05 | 2.958e-04 | 1.552e-04 | 5.657e-04 | 6.919e-04 | 1.000 |
| 600 | Bp | 1.391e-04 | 0.0 (200) | 8.178e-05 | 8.182e-05 | 1.230e-04 | 6.786e-05 | 2.958e-04 | 1.552e-04 | 5.657e-04 | 6.919e-04 | 1.000 |
| 600 | E | 1.391e-04 | 0.0 (200) | 8.178e-05 | 8.182e-05 | 1.230e-04 | 6.786e-05 | 2.958e-04 | 1.552e-04 | 5.659e-04 | 6.920e-04 | 1.000 |
| 600 | G | 1.391e-04 | 0.0 (200) | 8.178e-05 | 8.182e-05 | 1.230e-04 | 6.786e-05 | 2.958e-04 | 1.552e-04 | 5.657e-04 | 6.919e-04 | 1.000 |
| 600 | H | 1.280e-04 | 90.0 (200) | 7.684e-05 | 7.688e-05 | 1.977e-04 | 1.137e-04 | 3.128e-04 | 1.662e-04 | 5.342e-04 | 6.488e-04 | 0.920 |
| 600 | L | 2.511e-06 | 45.0 (203) | 1.419e-06 | 1.419e-06 | 1.039e-05 | 5.705e-06 | 5.330e-06 | 3.315e-06 | 4.291e-04 | 6.286e-04 | 0.018 |

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

| arm | hst rows | t at c600 | dt | L_top/L_in | L_bot/L_in | KE_int | vs ref: rel L_top, rel KE, rows identical to ref |
|---|---|---|---|---|---|---|---|
| B | 1201 | 6.000000e+05 | 1000 | 0.996750 | 0.998785 | 6.010414e+40 | - |
| Bp | 601 | 6.000000e+05 | 1000 | 0.996750 | 0.998785 | 6.010414e+40 | B: 0.00e+00, 0.00e+00, 601 of 601 |
| E | 601 | 6.000000e+05 | 1000 | 0.996750 | 0.998785 | 6.010412e+40 | Bp: 2.98e-09, -3.16e-07, 0 of 601 |
| G | 601 | 6.000000e+05 | 1000 | 0.996750 | 0.998785 | 6.010413e+40 | Bp: 8.90e-10, -1.44e-07, 0 of 601 |
| H | 601 | 5.739130e+05 | 956.52 | 0.997294 | 0.999500 | 5.722459e+40 | Bp: 5.46e-04, -4.79e-02, 0 of 601 |
| L | 601 | 6.000000e+05 | 1000 | 0.996744 | 0.999509 | 5.083272e+40 | Bp: -5.91e-06, -1.54e-01, 0 of 601 |
- Bp vs B m4(200-400) at c600: +0.0 % (within 10 %)

## L: the rule (orion zzzzd4; residual radial pattern reported in the table, not scored)
- m4(200-400) c600: Bp 1.391e-04, L 2.511e-06
- factor Bp / L = 55.396 (needs >= 3); L_top/L_in -6.020e-06, KE_int -1.543e-01 relative to Bp (needs |.| <= 0.01); clean run and end on cycle limit: True
- **FAIL**

## L: restart check, fresh 60 vs 30 + rst + 30 (the deck + one rst output block in the test copy, dt 1e30; -r restart with only time/nlim = 60 named)
- L60: 2026-11-24 18:36:02 rc=0 (srun) ended: Terminating on cycle limit CW 0 FATAL 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- L30: 2026-11-24 18:38:56 rc=0 (srun) ended: Terminating on cycle limit CW 0 FATAL 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- L30b: 2026-11-24 18:39:39 rc=0 (srun) ended: Terminating on cycle limit CW 0 FATAL 0 last: <rad_vet> hesdirk2 wall split (rank 0 s): clean stage steps 
- user.hst: L60 61 rows, L30 31 rows, L30b 30 rows; rows of L30 not printed-identical to L60 at the same t: 0 (first t -); L30b: 0 (first t -) -> BITWISE (printed)
- hydro.hst: L60 61 rows, L30 31 rows, L30b 30 rows; rows of L30 not printed-identical to L60 at the same t: 0 (first t -); L30b: 0 (first t -) -> BITWISE (printed)
