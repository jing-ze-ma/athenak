# NOTE DeltaAI: BSG true repro R / RC production (TASK-2026-10-11-deltaai-bsgrec-prod) -- NOT RUN, superseded

DeltaAI did **not** start this production (0 GPU-h). It was held for a user decision (DeltaAI's role is short tests), and
meanwhile Raven finished arm R to t = 6.48 d (job 31078702, 8c206e54) and arm RC moved to Raven with conservative2 + drag.
DeltaAI's contributions to this line: the smoke (NOTE-2026-10-11-deltaai-bsgrec-smoke: R clean, old conservative unhealthy, IC bitwise),
the full-size bsgaccel timing + conservative2 check (NOTE-2026-10-12-deltaai-bsgaccel-fullbsg), and the conservative2 + drag gate
(NOTE-2026-10-12-deltaai-consfix-c2drag: OK on the full BSG and AG Car A). Closed on the DeltaAI side 2026-10-11 22:37 CDT.
