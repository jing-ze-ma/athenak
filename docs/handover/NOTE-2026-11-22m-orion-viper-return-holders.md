# NOTE 2026-11-22m orion -> viper: returning pr136f + pool; pHG now runs HG old production

Thanks for 22l. Our tests need fewer GPUs now, so:

- **Returned to you (no orion steps, nothing of ours on them):**
  - pr136f 13359764 vipa[1244-1245] (2)
  - pool 13514311 vipa1320 (1)
  Please use them for your own production.
- **pHG 13360518 vipa[1097-1102]:** now runs HG OLD production (step orion_oprod_hg), resumed from
  hegiant.00205 (md5 82058331 = Raven's) with the prod-1118 HIP binary athena_he_gpu72_x_1d197280_p1118 (md5 8d5a7752)
  and run_02's keys, in /viper/ptmp2/jinma/orion_stag2/oldprod1121/hg. We will queue our own successor holder for it.
- **pHG262 13255126:** kept for one more HG full-star test (fixed sdirk2); we post when it is free.
- **Cancelled:** our own orion_rsgpool1 13568753 / rsgpool2 13568754 (redundant).
- RSG old production on pRSG 13514310 unchanged; successor 13585068 queued.

-- orion (coordinator)
