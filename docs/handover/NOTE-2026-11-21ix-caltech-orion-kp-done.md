# NOTE 2026-11-21ix Caltech -> orion: kp restart DONE (USER-confirmed "yes switch all four"), all from the PRE-cfl-0.9 rsts (USER)

B_fwbs = athena_fwbs_f496d606b (caltech-test-fwbfast f496d606b, CUDA 12.9 + g++ 13, in place, 53 TUs ~4 min), md5 b2ae1be413b527fab02023f159832612.
Smoke gate: pass 1 5.875e-06, NOT converged after 50, FATAL 0 (healthy; differs from 4.020e-06 by the fwb/fcone fresh defaults).
Keys: rst keys + conv_max 1e-4 + face_opac harmonic + stag_opac_update true + stag_flim false + stag_fwb true + stag_fcone true + stag_fwb_stage true + stag_thinf_held false + cfl 0.9. No damp, no gate.
USER: restart from the rst from BEFORE the cfl-0.9 switch, so the cfl-0.9 (no fwb_stage) segment is discarded; any rst written in that segment moved to <run>/rst_cfl09/ (kept).

| run | node | restarted from | relaunched (PST) |
|---|---|---|---|
| m35mms (first, per kq) | 02-10 | m35mms.00015 | 14:57 |
| t9r | 02-03 | t9l62.00033 | 14:58 |
| i1_m35mms_f1k | 01-17 | m35mms.00014 (00015 from the 0.9 segment moved aside) | 14:58 |
| t9l62 | 02-09 | t9l62.00031 (00032, 00033 from the 0.9 segment moved aside) | 14:59 |

First 50+ cycles per run follow. Wall END 17:25 (holders end 17:46).

-- Caltech
