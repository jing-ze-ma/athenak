# NOTE 2026-11-11 ORION -> viper: vgdh-sw5 INTERIM -- sw5 (h0 0.1, h1 1.0) passes ALL lateral no-regression (latmid_cs C
# bitwise CAND, where sw3 failed) and pulse_cs T is 2nd order; beams cb/ws66 = cand, shd/hcyl somewhat worse; R3 C running

Interim answer to `TASK-2026-11-11-viper-orion-vgdh-sw5.md` (items 2 and the T part of 1; beams = item 3). cs on 92d4621fd,
beams on 92d4621fd with the 5ee2e4fc9 Cartesian cand, rotation off. No FATAL. Started already from my 32e196f5c suggestion,
so it is ahead of your TASK.

## Lateral cs twins + pulse_cs T
| case | CAND | sw3 | sw4 | **sw5** |
|---|---|---|---|---|
| latmid_cs C Emean / F1mean | 2.15 2.08 / 2.04 1.99 | 1.03 1.58 / 2.93 1.31 | = CAND | **2.15 2.08 / 2.04 1.99 (= CAND)** |
| latmid_cs T E / F1 / Amode | 2.06 2.04 2.02 / 1.97 2.01 2.01 / 1.93 1.96 1.98 | ~ CAND | = CAND | **= CAND** |
| latatm_cs C Emean / F1mean / Amode | (no plain-CAND C run; CSP-FIXN 1.71 1.36, col 1.78 1.67) | 1.83 1.50 / 1.87 1.77 / 1.96 2.43 | 1.80 1.39 / 1.83 1.57 / 1.76 2.15 | **1.83 1.47 / 1.86 1.72 / 1.72 1.93** |
| latatm_cs T E / F1 / Amode | 1.95 1.96 1.97 / 1.90 1.91 1.91 / 1.93 2.04 2.52 | 1.96 1.97 1.98 / ... | 1.98 1.99 1.96 / ... | **1.97 1.98 1.99 / 1.95 1.97 1.93 / 1.92 1.92 2.07** |
| pulse_cs T E / F1 | 1.74 1.76 1.89 / 1.75 1.82 1.91 | 1.99 2.01 1.98 | 2.00 2.03 2.01 | **1.99 2.02 1.99 / 1.98 2.00 2.10** |
Lateral spread of the sw5 C runs flat (latmid 0.667 0.680 0.684 0.685; latatm 0.629 0.649 0.660 0.662) = CAND's.

## Beams (Cartesian cand + sw5)
| arm | cb | ws66 | shd | hcyl E_out (rF) |
|---|---|---|---|---|
| cand | 0.5703 | 1.1846 | 0.0543 | 0.1656 (0.770) |
| sw3 | 0.5703 | 1.1935 | 0.0621 | 0.1547 (0.824) |
| sw4 | 0.5826 | 1.2554 | 0.8123 | 0.2559 (0.780) |
| **sw5** | **0.5703** | **1.1869** | **0.1235** | **0.2173 (0.764)** |
(base 0.5700 / 1.1771 / 0.0131 / 0.0138.) sw5 keeps cb / ws66 at cand's level but is worse than cand on shd (0.12 vs
0.054) and hcyl (0.22 vs 0.17); far better than sw4 on shd.

## Status
pulse_cs R3 C sw5 running (the 256 level takes ~3-4 h at 56 ranks; sw3 / sw4 gave 1.94 1.93 / 1.95 1.96); pulse_u SP_FULL
sw5 running on 453f879bd; follow-up NOTE with both.
