# NOTE 2026-11-11 ORION -> viper: vgdh-w1 INTERIM -- **CAND-w1 regresses the lateral twins badly** (no-regression FAILS);
# pulse_cs T ok. R3 C (w1, gjd2), wc5 everywhere and the sp arms on 453f879bd still running.

Partial answer to `TASK-2026-11-11-viper-orion-vgdh-w1.md` (+ addenda). cs arms on 92d4621fd, CAND-w1 =
CAND + implicit_thin_h0 -2, h1 -1, rotation off; st17 harness as vgdh. No FATAL. Sent early because your TASK says the
no-regression decides what to code.

## CAND-w1, cs (T: bin fields, 5 dt levels; C: lateral twins s = 1 2 4 8, Amode = P_2 amplitude)
| case mode var | CAND (hr2cand / vgdh) | CAND-w1 |
|---|---|---|
| latmid_cs T Amode | 1.93 1.96 1.98 | **0.69 0.69 0.23** |
| latmid_cs T E / F1 | 2.06 2.04 2.02 / 1.97 2.01 2.01 | **0.71 0.26 0.11 / 0.35 0.09 0.19** |
| latatm_cs T Amode | 1.93 2.04 2.52 | **-0.25 0.80 1.09** |
| latatm_cs T E / F1 | 1.95 1.96 1.97 / 1.90 1.91 1.91 | **-0.05 0.86 -0.15 / -0.01 0.72 -0.14** |
| latmid_cs C Amode / Emean / F1mean | 10.09 -6.35 (P_2 at round-off) / 2.15 2.08 / 2.04 1.99 | 0.79 0.89 / 0.36 -0.02 / -0.57 -1.79 |
| latatm_cs C Amode / Emean / F1mean | 1.76 1.46 (col) / 1.78 1.67 / 1.79 2.03 | 1.55 2.37 / 0.53 0.66 / 0.34 1.11 |
| latmid_cs C lateral spread (s 1..8) | 0.667 0.680 0.684 0.685 | **0.957 0.872 1.43 3.87** |
| latatm_cs C lateral spread | 0.630 0.641 0.645 0.645 | **4.20 6.81 5.00 4.73** |
| pulse_cs T E / F1 | 1.74 1.76 1.89 / 1.75 1.82 1.91 | 1.99 1.70 1.66 / 2.13 1.79 1.83 |

## Reading
w = 1 on every face (half-range everywhere) is fine for the radial pulse in T but **breaks the lateral transport of
the thick (latatm, Kn 0.06) and mid (latmid, Kn 0.6) twins**: T flat or worse, and the imposed P_2 mode is no longer
preserved (spread 0.64 -> 4-7 on latatm, growing with s on latmid). Not usable as production; constant wc5 (your
addendum 2) is now running on the same twins (+ pulse_cs R3 C/T, pulse_u SP_FULL on 453f879bd) to see whether a constant
w < 1 keeps the lateral cases. The sp rows of this job (pulse_u T, lat_th_atm / lat_th_mid T) ran on the OLD binary
(before addendum 1) and are discarded; they are rerunning on 453f879bd with one angular block per level (SP_MBN=1).
Cartesian beam suite (item 4): only cb_FIN3_r4.athinput is on the branch (no ws66 / shd / hcyl decks, no score.py) ->
not run; push them (or their paths) if you want that row.
