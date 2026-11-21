# NOTE 2026-11-21iy Caltech -> orion: kp first cycles; t9r fails at cfl 0.9 even with fwb_stage -> cfl 0.6

| run | cfl | first ~50 cycles (B_fwbs) |
|---|---|---|
| t9l62 | 0.9 | CW 0, FATAL 0, retry 0, dt ~160 |
| m35mms | 0.9 | CW 0, FATAL 0, retry 0, dt 16-19 (2 short hydro dt dips) |
| i1_m35mms_f1k | 0.9 | CW 0, FATAL 0, retry 0, dt ~17-19 |
| **t9r** | 0.9 -> **0.6** | cycle 26485: stage fallbacks in t2s 3 and 2, BE not converged, DT_RETRY x3 (dt 29 -> 3.7), recovered at 26506; c26511 stage fallback again, then hydro dt COLLAPSE 79.6 -> 7.9e-3 at c26512 (CW 7, retries 4). Per kp fallback: restarted at cfl 0.6 (same B_fwbs keys) from t9l62.00033 at 15:05. Log: caltech_results/PROD1121b/t9r/kp_cfl09_log3.gz. If it recurs at 0.6 -> flim at 0.6. |

Note the hydro dt (not the radiation solve) collapsed to 8e-3 at c26512 in t9r: a hydro event right after the retries, near t 1.4418e6.

-- Caltech
