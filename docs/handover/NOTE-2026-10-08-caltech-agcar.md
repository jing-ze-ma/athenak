# NOTE 2026-10-08 Caltech: AG Car A + B smokes queued (TASK-2026-10-08-caltech-agcar)

From: Caltech session. To: viper, Raven, DeltaAI. Docs only.

## Status: build + smokes queued, NO runs yet
- Files: agcar-files-1008 copied to /resnick/groups/carnegie_poc/jingze/agcar_1008/files; SETUP.sh done (inputs with the
  10-08 A change: tfloor_kelvin 3.0e3, eos_logt_min 3.4).
- Build 4235599 (running, CPU node): incremental GPU build at 2fbc6aa1, PROBLEM=he_star_m1, CUDA 12.9 / H200
  -> builds/bin/athena_gpu_he_star_m1_2fbc6aa1 (md5 to follow).
- Smokes, afterok on the build: 4235600 (A) and 4235601 (B); `time/nlim=10`, 1 H200, 1 rank, 4 MeshBlocks on that GPU.
  To be compared with the DeltaAI TASK section 5 table (viper 12131991 A / 12130898 B).
- The B and A runs (step 5 of the TASK) are NOT queued: the user decides after the smokes.

## Other copies
Raven (31004544 B, 31004545 A), viper (12131308 B, 12131309 A) and DeltaAI keep their copies. Do not cancel anything
because of this note; a later Caltech NOTE will say when a Caltech run actually starts.

## UPDATE 10-08 08:35 PDT: both smokes PASS
Binary builds/bin/athena_gpu_he_star_m1_2fbc6aa1 md5 5631b59b6a48f3a42d4320fc9a89433f (CUDA 12.9, H200; build 4235599).
Smokes 4235600 (A) / 4235601 (B), 1 H200, rc 0, 0 FATAL, 0 NON-CONVERGED.

| | A Caltech | A viper 12131991 | B Caltech | B viper 12130898 |
|---|---|---|---|---|
| t after 10 cycles | 1.2525591814773266e+04 | ...279e+04 | 2.4186988199012371e+03 | ...844e+03 |
| dt | 1.2525587305197887e+03 | ...751e+03 | 2.4186976857008051e+02 | ...841e+02 |
| mass (col 3) | 6.1192879615024131e+32 | ...7439e+32 | 1.7891809043881911e+31 | ...3109e+31 |
| tot-E (col 7) | 1.6552578844678480e+47 | ...9151e+47 | 1.0633773260127060e+46 | ...7629e+46 |
| IC T check | 7.28e-14 | 7.1e-14 | 5.68e-14 | 6.4e-14 |
| he_ic_balance max | 3.24e-7 | 3.24e-7 | 1.08e-7 | 1.08e-7 |
| Picard mean / max | 14.2 / 16 | 13.5 / 16 | 6.9 / 24 | 6.9 / 24 |
| zone-cycles/s | 3.33e6 (1 H200) | 3.06e6 (2 MI300A) | 3.80e6 (1 H200) | 3.18e6 (2 MI300A) |

All hst differences <= 5e-15 relative (FMA level). Speed: 1 H200 ~ 1 viper node; A 2.36 s/cycle, B 2.07 s/cycle.
Wall estimate at the smoke dt (dt will change): B ~1900 cycles ~1.1 h, A ~6400 cycles ~4.2 h on 1 H200.
Runs still NOT queued on Caltech: waiting for the user.

## UPDATE 10-08 08:45 PDT: runs QUEUED (user OK), not started
1 H200 each, script agcar_1008/link_agcar.sh, run dirs agcar_1008/run{B,A} (DONE/STOP/CANCEL files):
B 4235665 (12 h, tlim 4.624e5) -> A 4235666 (24 h, afterany B, tlim 8.064e6) -> A 4235667 (24 h, afterany).
Start estimate: none yet (Priority). A later update will say when B actually starts; keep other copies until then.
