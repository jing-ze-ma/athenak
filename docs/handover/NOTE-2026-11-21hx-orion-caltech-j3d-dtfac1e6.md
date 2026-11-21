# NOTE 2026-11-21hx orion -> caltech: 2 short arms — does the 1e6 steady start work WITH the realizability bound? (USER: run it)

On B_j3 (joint3 a5257b7ee), production layout, free test nodes, alongside the continuing J3 arms:
  J3D_t19 : as J3_t19 (t19l64_wngl_smoke + joint3 keys) but rad_m1/clean_init_dtfac = 1.0e6, tlim 6.0e3 s
  J3D_t9  : as J3_t9 but rad_m1/clean_init_dtfac = 1.0e6, tlim 6.0e3 s
Report: the clean_init steady lines (passes, converged or not), dt every 20 cycles, CW/FATAL, hot cells / vceil, FLIM counts, KE, L_bot/L_in, L_top/L_in. Results to caltech_results/J3D/.
Question: if these match J3 (no floor runaway), the code default can stay 1e6 and the realizability bound alone protects the stars.

-- Orion coordinator
