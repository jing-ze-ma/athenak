# NOTE 2026-10-13 from Caltech: open-ot CVALL (re TASK-2026-10-13-caltech-open-ot-cvall.md)

## QUEUED (10-12 ~22:45 PDT)
- OLD = athena_cpu_omp_281dd686, md5 9ae9f09aeb81841a275e8d046a0eebe4; NEW = 281dd686 + patch/0001 (CV binary), md5 dc58375a9206e2c034fc5d57cfd8ac9d.
- Debug QoS chain (1 node x 32 cores, 28-min cutoff each, finished steps skipped): CV rounds 4429421 -> 4429434 -> 4429435, then
  CVALL rounds 4429436 -> 4429437 -> 4429438 -> 4429439 -> 4429440 -> 4429441.
- run_cvall.sh unchanged except /usr/bin/time -> stand-in; P = 6 (each dir runs OLD + NEW side by side at 2-4 threads; 16 would oversubscribe 32 cores ~4x).
  m1_shock_ref in the 5 rs decks -> our copies of caltech_od/shockref/ref_m*.txt.
- FULL steps longer than ~28 min cannot finish in debug (restart from scratch); they will be listed. Interim bit_all counts here every few hours.
