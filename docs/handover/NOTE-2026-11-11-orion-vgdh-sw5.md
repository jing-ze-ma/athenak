# NOTE 2026-11-11 ORION -> viper: vgdh-sw5 cs COMPLETE -- **sw5 is 2nd order in space AND time on pulse_cs through the 256
# pair (R3 C E 1.90 1.99 / F1 1.88 1.87) and passes every lateral no-regression** -> the best cs setting so far

Completes item 1 of `TASK-2026-11-11-viper-orion-vgdh-sw5.md` (interim 7ffd3cade had 2, 3 and the T part of 1).
pulse_cs C R3 (N per face nx1/4, nx1 32..256, 56 ranks), CAND + implicit_thin_h0 0.1 implicit_thin_h1 1.0, 92d4621fd,
rotation off. 4/4 levels rc 0, no FATAL.

| arm | R3 C Emean (32-64-128, 64-128-256) | R3 C F1mean | spread (32..256) | pulse_cs T E / F1 | latmid_cs C | lat twins T |
|---|---|---|---|---|---|---|
| CAND | 1.07 1.71 | 1.29 1.59 | 0.157 0.626 0.713 0.598 | 1.74 1.76 1.89 / 1.75 1.82 1.91 | 2.15 2.08 | ~2 |
| sw3 (0.0 / 1.0) | 1.94 1.93 | 1.87 1.82 | 0.154 0.452 0.540 0.514 | 1.99 2.01 1.98 | **1.03 1.58 (FAIL)** | ~2 |
| sw4 (0.1 / 1.5) | 1.95 1.96 | 1.87 2.01 | 0.172 0.420 0.506 0.519 | 2.00 2.03 2.01 | 2.15 2.08 | ~2 |
| **sw5 (0.1 / 1.0)** | **1.90 1.99** | **1.88 1.87** | 0.163 0.462 0.557 0.536 | **1.99 2.02 1.99 / 1.98 2.00 2.10** | **2.15 2.08** | **~2** |

Beams (Cartesian cand + arm): cand 0.5703 / 1.1846 / 0.0543 / 0.1656; sw3 0.5703 / 1.1935 / 0.0621 / 0.1547; sw4 0.5826 /
1.2554 / **0.8123** / 0.2559; **sw5 0.5703 / 1.1869 / 0.1235 / 0.2173 (rF 0.76)**.

## Reading
sw5 = sw4's lateral-safe start (h0 0.1: w = 0 in the mid-opacity cells, latmid C bitwise CAND) + sw3's beam-friendlier end
(h1 1.0: w reaches 1 at g = 1) and keeps the 256 pair at ~1.9-2.0 in both E and F1. Its one cost vs CAND is on the opaque
beam decks (shd 0.054 -> 0.12, hcyl 0.17 -> 0.22), i.e. on top of plm-on-h's known opaque-body cost; far from sw4's shd
failure. pulse_u SP_FULL sw5 (453f879bd) still running -> follow-up.
