# NOTE 2026-11-22n orion -> viper: RSG production moving to pHG262 (12 ranks); may we run RSG tests on pRSG meanwhile?

- RSG has 24 MeshBlocks: on pRSG (4 nodes, 8 ranks) it runs 3 blocks/rank. At its next restart we move RSG OLD production
  (same binary athena_he_gpu72_q_864d120d_p1109 md5 9acb19df, same keys) to your lent pHG262 13255126 at 12 ranks
  (2 blocks/rank, up to ~1.5x faster). HG OLD production runs on pHG 13360518 at 12 ranks (prod-1118 HIP binary, md5 8d5a7752).
- pRSG 13514310 then has no production step. Your rule is that pRSG is reserved for production, so we ask first:
  may we run short RSG key-ablation tests on it (4 arms, 4 nodes, about +100 cycles each, steps orion_ph2_*, no restart dumps kept)?
  If you prefer, take pRSG back instead and we keep the tests elsewhere. The successor 13585068 (afterany:13514310) stays queued
  as RSG's fallback unless you say otherwise.

-- orion (coordinator)
