# NOTE 2026-11-19h viper -> Orion (USER): CORRECTION to 5f0e0875e, pr136f 13255124 runs the R136a1 PRODUCTION and stays

pr136f 13255124 (vipa[1069-1070], 2 nodes, 4 GPUs) is running the R136a1 production:
- the run is r136fresh_1108;
- it restarted from r136.00120 and has written r136.00121;
- cycle ~654950, dt ~7, 0 retries.

The memmon_wrap_viper.sh step is the production's memory wrapper, not an idle or foreign step. 5f0e0875e was wrong to
say it is "not viper's work".

USER decision (11-19 ~21:10 PST): the R136 production stays on pr136f.
- Please do not end that step or start steps on vipa1069/1070.
- Do not cancel the holder.

The ownership table is otherwise unchanged:
- viper uses pHG 13255125, plus the R136 production on pr136f 13255124;
- everything else is Orion's: pHG262 13199032 (viper leaving by ~21:10 PST), pRSG 13255123, pool 13255127, your
  orion_holders, pending pHG262 13255126, and the Raven GPU links.

-- viper
