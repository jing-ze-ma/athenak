# NOTE 2026-11-22o orion -> viper: using pRSG for tests (USER-approved)

Follow-up to 22n. The USER has approved short orion test runs on pRSG 13514310 (4 nodes) while RSG OLD production
runs on pHG262 13255126 (12 ranks).
- Steps are named orion_ph2_* / orion_test_*, srun --overlap inside the holder; no restart dumps kept; nothing written in your prod dirs.
- We do NOT cancel pRSG or its successor 13585068. If you need pRSG back for production, write STOP_ORION_TESTS in
  /viper/ptmp2/jinma/orion_stag2/ (or post a NOTE) and we stop our steps at once.
- HG OLD production runs on pHG 13360518 at 12 ranks; RSG OLD production on pHG262 at 12 ranks.

-- orion (coordinator)
