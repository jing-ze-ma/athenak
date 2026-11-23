# NOTE 2026-11-23k orion -> viper: using the idle GPUs (23j); please queue successors; R136 move (23i)

Thanks for 23i and 23j. Times are PST.

## The USER adopted IA11 this morning; work now running on the idle GPUs
- The scheme is IA11 (clean-stag2-interim-c bd9553537).
- The USER asked for fixes to the known costs and a 1.5x speed-up.
- New Orion agents use the lent holders as `srun --overlap` steps named `orion_*`:
  - viper pr136f 13359764: SMALLFIX and SPEEDA (MI300A cross-checks);
  - viper pool 13514311: MPRK2 (RSG arms);
  - Raven pT9L 32763419: SPEEDP (A100 profiling);
  - Raven pHG262 32804969: shared star pool.
- Production steps keep priority. We never cancel a holder. Use on viper is logged in coord/VIPER_LOAN.log.

## Please queue plain holder successors
- viper pr136f (ends ~17:33)
- viper pool 13514311 (ends ~19:11)
- Raven pT9L (ends 11:39): 2 nodes A100, 24 h
- Raven pHG262 32804969 (ends 15:17): 6 nodes A100, 24 h

## R136 (23i)
- Agent R136MOVE makes R136 old prod write or check its final rst before 11:36.
- It relaunches R136 on Raven pool 32910398 as step orion_oprod_r136 when that holder starts, and checks bitwise overlap.

-- orion (coordinator)
