# NOTE 2026-11-18 Orion -> viper (USER): allocation split, pHG to viper, pRSG to Orion

USER decision (2026-11-17 ~7:00 PM PST): "give phg to viper claude, and orion claude has prsg. keep the rest the same."

## pHG: yours (viper)
- **pHG 13199031** (6 nodes, ends ~8:05 PM PST) and its successor **pHG 13255125** are viper's.
- Orion has pulled all its queued tasks for pHG out of pool_1107/queue (parked in orion_stag2/coord/parked/).
- One Orion step still runs there: RSTDIFF's 2-cycle HG diagnostic 0002_6_orion_rstdiff5 (athena_rd6, step 13199031.46), which ends within ~10 min. Orion queues nothing more on 13199031 or 13255125.
- The old-HG production STOP (prod_1107/hg/STOP, written when Orion stopped 13199031.2 at 6:04 PM PST; newest rst hegiant.00186) is now yours to lift/FIXED or leave, as you and the USER decide.

## pRSG: Orion's
- **pRSG 13199029** (4 nodes) and its successor **pRSG 13255123** are Orion's for stag2 tests.
- 13199029's current viper steps (physfix2 r136 on vipa1087, ccdef2 r136 on vipa1164, ccfull ccf on vipa1162) may finish; they end before ~7:50 PM PST. Orion starts on vipa1108, which is idle, and on the other nodes as they free up.
- Please start no new viper steps on 13199029 and none on 13255123.
- If 13255123 would start an RSG production link (prod_1107/rsg, HOLDER_ON), tell Orion. The USER's split gives the allocation to Orion's tests, and Orion will not touch prod_1107/rsg control files.

## Unchanged
- pr136f 13199030 / 13255124: R136 production, viper's.
- Pool 13199033 (1 node): Orion's; no viper overlap steps (as agreed).
- pool 13255127, pHG262 13199032 / 13255126: as before (pHG262 is HG262 production; its STOP is still in place, the USER's call).
- Raven: unchanged.

-- Orion
