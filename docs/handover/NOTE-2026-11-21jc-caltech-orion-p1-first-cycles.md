# NOTE 2026-11-21jc Caltech -> orion: first cycles on B_p1 (16:41 PST); t9l62 -> cfl 0.6

| run | cfl | status at 16:40 |
|---|---|---|
| m35mms | 0.9 | 42.7 % of tlim, c45150, dt 16.4, CW 0, FATAL 0, retry 0 (3 short hydro dt dips). Not yet at c45706 (the old FATAL cycle). |
| i1_m35mms_f1k | 0.9 | 45.1 %, c43250, dt 23, CW 0, FATAL 0, retry 0 (3 short hydro dt dips) |
| t9r | 0.6 + tfloor 1000 K | 54.9 %, c26850, dt 48, CW 0, FATAL 0, retry 0 |
| **t9l62** | 0.9 -> **0.6** | at cfl 0.9: hydro dt COLLAPSE 74 -> 15 at c29438, then CW 1 at c29469 (t2s 2 stage fallback 1) within the first 100 cycles -> per the fallback rule restarted at cfl 0.6 from t9l62.00035 at 16:41 (B_p1, same keys) |

Holders end 17:46 PST; the runs pause at their last periodic rst then (successors still PENDING, est. 12-01..12-03). Tell us if the USER wants them moved to viper / Raven.

-- Caltech
