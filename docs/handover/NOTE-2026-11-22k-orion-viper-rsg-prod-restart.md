# NOTE 2026-11-22k orion -> viper: RSG old production restarted on pRSG 13514310

Times PST.
- RSG old production (prod-1109 864d120d, athena_he_gpu72_q_864d120d_p1109 md5 9acb19df, rsg.keys, 8 ranks 2/node) is
  running on its own pRSG 13514310 (vipa1124,1126,1240,1276) since 20:01, step orion_oprod_rsg, from rsg.01722 (cycle 778667).
  Dir /viper/ptmp2/jinma/orion_stag2/oldprod1121/rsg (driver oprun2.sh, keep 3 rst, stall watchdog). At 20:05: cycle 778711,
  dt 916.67, 0 FATAL, L_top 6.383e38 (pre-restart 6.37-6.38e38).
- No Orion production on the lent pHG 13360518 / pHG262 13255126 / pr136f 13359764 / pool 13514311 — test steps only
  (orion_* names, logged on our side). HG old production stays on hold (USER: tests first).

-- orion (coordinator, for the USER)
