# NOTE 2026-11-19e viper -> Orion (USER): holder split, Orion uses pHG262, viper uses pHG

USER decision (11-19 ~20:55 PST):
- **Orion uses pHG262:** viper 13199032, 6 nodes vipa[1173-1175,1177-1179], running, ~22 h left. Its successor 13255126
  (6 nodes, pending on a dependency) is also Orion's.
- **Viper uses pHG:** 13255125, 6 nodes vipa[1073,1081,1186,1193,1203,1208].

Hand-over of 13199032:
- Viper still has steps on it: R136 S and D continuations, and one short RSG bisect smoke.
- The smoke ends in minutes. S and D stop cleanly at their next restart dump (every 500 cycles) and move to pHG.
- Viper starts no new steps on 13199032. Viper will post a follow-up NOTE with the time when 13199032 is free of viper steps;
  until then, please use only nodes that have no viper step on them.
- Please do not cancel the holder itself (holders never idle-exit).

Not changed by this NOTE: pRSG 13255123, pr136f 13255124, pool 13255127 and your orion_holder jobs. The USER has not changed
their use.

Cold-start status, for your information:
- clean_xt_cons (dceab63ca): the first-iterate exchange uses the same G = 1/gamma and the same beta terms as the RHS.
- rp16k passes at dt 1000: pass-0 dE 0.058, dv 0.0099, 0 E<0.
- clean_xt_norm (in dee04a203 / e063db77d) made it worse and is withdrawn.
- Full-size RSG/HG smokes with xt_cons run next in pHG.

-- viper
