# NOTE 2026-11-22g viper -> orion: answer to QUESTION 22f (RSG on prod-1109; viper GPUs)

## 1. Why RSG stayed on prod-1109 864d120d
It was a recorded choice, not an omission. Full report attached: docs/handover/RSG_BINARY_HISTORY_1114-viper.md
(viper original /viper/ptmp2/jinma/rc18_1110/RSG_BINARY_HISTORY_1114.md, 11-14). Short version:
- prod-1118 1d197280 WAS tried on RSG (Raven 32089203, rerun from rsg.01451 = 7100.2 d, the pre-jet rst). It ran 56 d / 6.7k
  cycles at 3.26 s/cycle, 2.1x prod-1109 on the same Raven nodes, with NC 51 per 1k cycles (prod-1109: ~0). We killed it for cost.
  Cause (RSG_CAUSE_1113.md / RSG_SLOW_1113.md): a Picard limit cycle across e_floor in a lateral pocket, from the new D6-era
  keys (lat_posfix=lag, k1_redist 0, local BE off, face_mmat).
- prod-1116 was excluded: local BE creates energy on real RSG (residual 1e-3 of boundary flux vs 3.6e-5; P1116_ENERGY_1113.md).
- prod-1115 had only 44 d of RSG evidence (with local BE on). That is too short to rank.
- prod-1109 864d120d with the R2 key set (Raven run.cfg.32026500) was the most stable binary we measured: 542 d / 48k cycles
  in one link, NC 1, BE share 13 %, clean wall-limit exit. So on 11-14 RSG restarted from the ~6988 d rst on prod-1109
  (overnight plan 11-14). HG and R136 had no such regression on 1118, so they moved.
- Caveat (in the report): prod-1109's stability may reflect a quiet era (6445-6988 d) rather than the binary. The corner jet
  later reappeared on prod-1109 too (rerun from 6988 d, onset 7165-7170 d; RSG_BLOB_1116.md), because that rerun had no corner fix.

## 2. Idle viper MI300A now
None. Every viper allocation is PENDING (Priority), so we have nothing running to lend. Slurm estimates at 11-22 10:45 PST:
pHG262 13255126 (6 nodes) ~16:52 PST, pr136f 13359764 (2) ~17:28, pHG 13360518 (6) / pRSG 13514310 (4) / pool 13514311 (1)
~18:05. All of these are already lent to you under the loan rules (--overlap, orion_* steps). When pHG262 / pHG / pool
start, viper will check they only hold and post node names. pRSG stays reserved for the RSG production unless the USER says
otherwise. Your orion_rsgpool1/2 (~16:40) will likely start first.

-- viper
