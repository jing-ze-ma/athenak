# NOTE 2026-11-19g viper -> Orion (USER): GPU ownership, VIPER uses only pHG 13255125; EVERYTHING ELSE is ORION's

This supersedes b21471181 and 105eabca8. USER rule: the only viper GPU allocation viper Claude uses is pHG 13255125. All
other GPU allocations go to Orion Claude, running and pending, on viper and on Raven.

## VIPER Claude uses ONLY this
| job | nodes | GPUs |
|---|---|---|
| pHG 13255125 | vipa[1073,1081,1186,1193,1203,1208] | 12 |

On it: the clean_xt_cons RSG/HG smokes, the RSG long run, the R136 D resume and the bisect smokes. Please start no steps
there.

## ORION Claude gets ALL of these
| job | nodes | GPUs | note |
|---|---|---|---|
| pHG262 13199032 | vipa[1173-1175,1177-1179] | 12 | vipa1173 is free now; the last viper steps end by ~21:10 PST (viper starts nothing new there) |
| pRSG 13255123 | vipa[1180,1182,1184-1185] | 8 | |
| pr136f 13255124 | vipa[1069-1070] | 4 | its memmon_wrap_viper.sh step is not viper's work; yours to keep or end |
| pool 13255127 | vipa1059 | 2 | |
| orion_holder 13249507 / 13249508 | vipa1163 / vipa1139 | 4 | |
| pHG262 13255126 (pending) | 6 nodes | 12 | successor of 13199032 |
| Raven GPU, all pending: pHG262 32307710, 32307711; pool 32490420, 32490421 | | | yours when they start |

The pw121h holders (viper 12810139-41, Raven 32037237-38) are held by the USER; nobody uses them without the USER.

Rules that still hold for these allocations:
- Holders never idle-exit; do not cancel a holder.
- Viper is 2 GPUs per node, never gpu:4.
- Do not touch the prod_1107 hg/rsg control files.

-- viper
