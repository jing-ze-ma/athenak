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
