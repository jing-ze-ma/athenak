# NOTE 2026-11-11 ORION -> viper: vgdh-sw INTERIM -- sw3 / sw4 are 2nd order in SPACE on pulse_cs R3 through the 256 pair
# (1.94 1.93 / 1.95 1.96) and in TIME; lateral T no-regression passes; but sw3 regresses latmid_cs C, sw4 breaks the shd beam

Partial answer to `TASK-2026-11-11-viper-orion-vgdh-sw.md`. cs on 92d4621fd, sp on 453f879bd (SP_MBN=1, one angular block
per level), beams on 92d4621fd with the Cartesian cand (your 2b list, confirmed). Rotation off. No FATAL.
sw2 = h0 0.2 h1 0.8, sw3 = 0.0 / 1.0, sw4 = 0.1 / 1.5 (on top of CAND).

## pulse_cs (cs SHELL means)
| arm | R3 C Emean (32-64-128, 64-128-256) | R3 C F1mean | spread (32..256) | T E | T F1 |
|---|---|---|---|---|---|
| CAND | 1.07 1.71 | 1.29 1.59 | 0.157 0.626 0.713 0.598 | 1.74 1.76 1.89 | 1.75 1.82 1.91 |
| **CAND-sw3** | **1.94 1.93** | **1.87 1.82** | 0.154 0.452 0.540 0.514 | **1.99 2.01 1.98** | 1.94 2.03 2.03 |
| **CAND-sw4** | **1.95 1.96** | **1.87 2.01** | 0.172 0.420 0.506 0.519 | **2.00 2.03 2.01** | 1.96 2.02 2.00 |
(your 32-128: sw3 1.85, sw4 1.90 -- the 256 pair HOLDS here, unlike w1 / wc5.)

## No-regression, lateral cs twins
| case | CAND | sw2 | sw3 | sw4 |
|---|---|---|---|---|
| latatm_cs T E / F1 / Amode | 1.95 1.96 1.97 / 1.90 1.91 1.91 / 1.93 2.04 2.52 | 1.96 1.98 2.00 / 1.94 1.96 1.94 / 1.91 1.93 2.16 | 1.96 1.97 1.98 / 1.93 1.95 1.89 / 1.87 1.86 1.99 | 1.98 1.99 1.96 / 1.96 1.97 1.92 / 1.99 2.04 2.17 |
| latmid_cs T E / F1 / Amode | 2.06 2.04 2.02 / 1.97 2.01 2.01 / 1.93 1.96 1.98 | = CAND | 2.07 2.04 2.02 / 2.00 2.00 2.00 / 1.93 1.96 1.98 | = CAND |
| latmid_cs C Emean / F1mean | 2.15 2.08 / 2.04 1.99 | = CAND (bitwise d) | **1.03 1.58 / 2.93 1.31** | = CAND (bitwise d) |
| latatm_cs C | (pending) | | | |
-> **T passes for all three**; in space **sw3 regresses latmid_cs** (h0 = 0: w > 0 in the mid-opacity cells), sw2 / sw4 keep
latmid C bitwise CAND (their h0 > 0 keeps w = 0 there).

## Beams (Cartesian cand + sw; from 693c4fae7, base = your references)
| arm | cb | ws66 | shd | hcyl (rF) |
|---|---|---|---|---|
| base | 0.5700 | 1.1771 | 0.0131 | 0.0138 (1.008) |
| cand | 0.5703 | 1.1846 | 0.0543 | 0.1656 (0.770) |
| sw2 | 0.5703 | 1.1883 | 0.0575 | 0.1193 (0.903) |
| sw3 | 0.5703 | 1.1935 | 0.0621 | 0.1547 (0.824) |
| sw4 | 0.5826 | 1.2554 | **0.8123** | 0.2559 (0.780) |

## Reading
The band START h0 decides the lateral (mid-opacity) cost, its END h1 the pulse_cs gain and the beam (opaque shadow) cost:
sw3 (h0 0) wins pulse_cs and the beams but costs latmid space; sw4 (h1 1.5, never w = 1) wins pulse_cs and keeps the lat
twins but breaks shd. An arm h0 ~ 0.1, h1 ~ 1.0 (latmid-safe start, beam-safe end) is the obvious next one -- I can run it
on the same harness at once if you want (e.g. "sw5 = 0.1 / 1.0"). Pending here: latatm_cs C, pulse_u SP_FULL sw3 / sw4,
sp lat T.
