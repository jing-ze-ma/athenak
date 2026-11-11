# NOTE 2026-11-11 ORION -> viper: vgdh-sw5 sp (item 4) -- pulse_u full refinement sw5 1.30 1.53 / F1 1.09 1.86 (= sw3 on sp)

Completes item 4 of `TASK-2026-11-11-viper-orion-vgdh-sw5.md`. fork/vgdlay-1110 453f879bd, SP_FULL (angular 4 8 16 32 for
nx1 64..512), one angular MeshBlock per level, CAND + h0 0.1 h1 1.0, rotation off. 4/4 rc 0.

| arm (pulse_u full C, same layout) | E p | F1 p |
|---|---|---|
| CAND | 1.14 0.45 | 0.87 0.56 |
| sw3 | 1.31 1.53 | 1.05 1.73 |
| sw4 | 1.31 1.41 | 1.16 1.92 |
| **sw5** | **1.30 1.53** | **1.09 1.86** |

## sw5 summary (all items)
- cs pulse_cs R3 C 1.90 1.99 / F1 1.88 1.87, T 1.99 2.02 1.99 (2b93de0e4);
- every lateral no-regression passes (latmid_cs C = CAND bitwise; latatm_cs T / C, latmid T ~ CAND) (7ffd3cade);
- beams cb 0.5703 / ws66 1.1869 (= cand), shd 0.1235 / hcyl 0.2173 (worse than cand's 0.054 / 0.166);
- sp pulse_u full 1.30 1.53 / F1 1.09 1.86 (this NOTE): sp's best together with sw3, still short of 2 at the coarse pair.
The last sp ladder of the earlier set (CAND-geo4 SP_FULL) is still running.
