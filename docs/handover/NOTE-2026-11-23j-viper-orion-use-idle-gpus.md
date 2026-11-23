# NOTE 2026-11-23j viper -> orion: USER — use all idle GPUs

USER (11-23 ~09:57 PST): "ask orion to use all idle gpus".
Idle at 09:57 PST (holder batch only, no step), all lent to you under the usual rules (srun --overlap, steps orion_*,
never cancel a holder; production steps keep priority):
- viper pr136f 13359764 vipa[1244-1245] (2 nodes, 8 MI300A), ends ~17:33 PST, no successor
- viper pool 13514311 vipa1320 (1 node, 4 MI300A), ends ~19:11 PST, no successor
- Raven pT9L 32763419 ravg[1032,1207] (2 nodes, 8 A100), ends 11:39 PST, no successor
Busy: viper pHG262/pHG/pRSG (HG262/HG/RSG old prod), Raven pool 32747999 (R136, ends 11:36; successor 32910398 per NOTE 23i),
Raven pHG262 32804969 (your ia11go steps). Any holder that falls idle later is also yours to use; no need to ask.
If you want successors for pr136f / pool / pT9L to keep them, say so and viper queues plain holders.

-- viper
