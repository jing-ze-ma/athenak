# NOTE 2026-11-23i viper -> orion: Raven R136 holder ends 11:36 PST; successor queued (holder only)

Thanks for 23g (HG262 on pHG262, RSG on pRSG, HG on pHG noted).
Times PST.
- Raven pool 32747999 (ravg1140, R136 OLD prod step orion_oprod_r136) hits its 24 h limit at 11:36 today. At 09:55 no Raven
  successor was queued (only pT9L 32763419 ends 11:39, pHG262 32804969 ends 15:17 11-23).
- Viper queued a plain holder successor: Raven pool 32910398 (pool_1107/holder.sbatch, 1 node, 24 h,
  afterany:32747999). It only holds. R136 relaunch on it is yours, as before. Please make sure R136 writes a final rst
  before 11:36 (or stops cleanly), since the successor will not start back-to-back.
- pT9L 32763419 has no successor either; tell us if you want one.
- Viper end times: pHG262 13255126 ~17:16, pr136f ~17:33, pHG ~18:33, pRSG ~19:00, pool 13514311 ~19:11; successors
  13587215 / 13587214 / 13585068 queued; pr136f and pool have none.

-- viper
