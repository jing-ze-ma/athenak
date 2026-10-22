# NOTE 2026-10-22 caltech -> viper: CORRECTION to 7afcaf1d: opac-lag / convergence keys do NOT cure the photospheric runaway

Restart arms from a GROWING state (HG N=32, seed rst t=0.75 d -> 1.25 d; base photospheric E sd 4.6e-3 -> 1.2e-2):
implicit_opac_lag_passes = 4, = 99, implicit_tol_rel = 0 with opac_lag none (full convergence + always update), and tol_rel = 0 with
passes 2 ALL match base to 3 digits through 1.126 d (the keys did act: e.g. p4 "lagged solves 509"). The flat fresh-start arms
(passes 4, face_weight equal, opac_lag none) were smaller seeds from a different start-up transient, not cures.
=> The runaway lives in the CONVERGED discrete system: one-cell MLT deposit + under-resolved kappa cliff (your criterion w n dlnT < -1).
Not solver iteration. mlt_fsub_smooth n=2 and 2x nx1 slow it (e-fold ~0.7 d) but do not stop it.
Next at Caltech: a discrete face-flux fix that keeps the cliff cell's outflow monotone in its own T (Kirchhoff / T-integral-mean face
diffusivity across the cliff), tested on the restart setup. Suggest viper uses a restart-from-growing-state setup on RSG E too
(fresh-start comparisons are confounded by the seed).
