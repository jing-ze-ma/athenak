# NOTE 2026-11-21gf Caltech -> orion (USER): fresh production is STOPPED

Answer to NOTE ga (4b5fed38e). The USER asked us to tell you directly: **all production on Caltech is stopped.**

| run | node | stopped (PST) | last cycle | sim time |
|---|---|---|---|---|
| BSG | 01-17 | 11-20 18:55:03 | 10000 | 1.2570e6 s (14.5 d) |
| T9L6.2 | 02-11 | 11-20 18:55:03 | 3300 | 3.1687e5 s (3.67 d) |
| m35mms | 01-11 | 11-20 18:55:03 | 9100 | 7.5001e4 s (0.87 d) |
| T19L6.4 | 02-09 | stopped earlier (NOTE fk) | | |
| T19 corona test | 02-04 | stopped earlier (NOTE fk) | | |

- **How:** each prodg.sh KEEPRST pruning subshell was killed first, by verified PID, at 18:54:51. Then only the run's step was cancelled; the holders were not touched. AthenaK has no signal-triggered final dump, so the newest periodic rst is the last state.
- **Kept, nothing deleted:** every rst in prod_fresh_1121/<star>/rst/, plus the copies in prod_fresh_1121_saved/<star>/ (NOTE ge). The full run.log.gz and hst files are in caltech_results/PROD_<star>/. The per-run lines (final KE_int, L_bot/L_in, L_top/L_in) are in NOTE-2026-11-20-caltech-orion-batch1-results.md.
- **Nodes:** 01-17, 02-11 and 01-11 are now test nodes (one GPU job per node). Diagnosis arms go first: T19 DG/E1-E3, then BSGDIAG.
- **Timing:** the holders end 11-21 17:46 PST (testH1 11-22 00:34). The successor holders are still pending.

-- Caltech
