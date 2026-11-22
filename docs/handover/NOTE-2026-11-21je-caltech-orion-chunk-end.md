# NOTE 2026-11-21je Caltech -> orion: all 4 production runs paused cleanly at wall END (17:15 PST)

B_p1 graceful wall-time stop worked on every run ("stop requested (wall time) ... restart written", rc 0):
| run | final rst | cycle | CW / FATAL / retry (B_p1 chunk) |
|---|---|---|---|
| t9l62 (cfl 0.6) | t9l62.00037 | 30085 | 1 / 0 / 0 |
| t9r (cfl 0.6, 1000 K) | t9l62.00035 | 27579 | 0 / 0 / 0 |
| m35mms (cfl 0.9) | m35mms.00018 | 46957 | 3 / 0 / 0 (passed the old FATAL c45706) |
| i1_m35mms_f1k (cfl 0.9) | m35mms.00020 | 44843 | 2 / 0 / 0 |
Holders end 17:46; successors still PENDING (est. 12-01..12-03). Continue elsewhere only on the USER decision.

-- Caltech
