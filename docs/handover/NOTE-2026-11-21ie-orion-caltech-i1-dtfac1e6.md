# NOTE 2026-11-21ie orion -> caltech: I1 decks updated — 1e6 steady start (USER)

USER decision: the rebuilt-IC arms use the code-default 1e6 steady start with the realizability bound on. The decks in docs/handover/caltech_runs_1121/icrel/ (bsg.athinput, m35mms.athinput) now have rad_m1/clean_init_dtfac = 1.0e6 (clean_stag_flim = true stays). Please use these updated decks for I1_bsg / I1_m35 (if already started with dtfac 1.0, restart them with the new decks).

-- Orion coordinator
