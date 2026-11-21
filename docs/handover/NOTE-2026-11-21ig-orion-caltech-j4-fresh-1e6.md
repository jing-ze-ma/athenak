# NOTE 2026-11-21ig orion -> caltech: restart the 4 star tests FRESH with the 1e6 steady start (USER)

USER decision: run the four stars fresh with the default clean_init_dtfac = 1.0e6 plus the realizability bound.
1. Stop the four J3 continuations (J3_t19, J3_t9, J3_m35, J3_bsg; they use dtfac 1.0) — step cancel, commit partial hst, delete their own rst/bins.
2. On the same nodes start J4_t19, J4_t9, J4_m35, J4_bsg FRESH (t = 0) on B_j3 (joint3 a5257b7ee): the same *_wngl_smoke decks and keys as J3
   (rad_m1/clean_stag_flim = true, rad_m1/clean_stag_thinf_held = true, rad_m1/floorcap_rho = 1e-14, floorcap_v = 3e7, floorcap_report = 25)
   but rad_m1/clean_init_dtfac = 1.0e6 (replace the 1.0 line in place).
   tlim: J4_t19 1.0e5 s, J4_t9 5.0e4 s, J4_m35 3.0e4 s, J4_bsg 1.12e6 s (13 d). rst every ~1 h wall (keep 2), m1prof/vetprof at the end and at
   BSG 6 d. Progress every ~20 min (dt, KE, L_bot/L_in, L_top/L_in, CW, FLIM/floorcap lines, clean_init steady line).
These are TESTS on test nodes (not production). I1_bsg / I1_m35 and the C2F arms keep running.

-- Orion coordinator
