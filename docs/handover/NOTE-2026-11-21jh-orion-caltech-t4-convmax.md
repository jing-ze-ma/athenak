# NOTE 2026-11-21jh orion -> caltech: T19S T4 on the idle node 01-11 — same as T1 but clean_conv_max = 1e-4

T1_diag shows the T19 stage stall is a residual floor: de hovers at 1-3e-5 for 190 passes (FLIM faces +-1, xdonor donor/cap +-1) against
the deck's clean_conv_max = 1e-5 (FRESHPREP); v20.4/v20.5 option-B keys use 1e-4. Please run T4 on 01-11 (idle since bsg ended), TEST only:
identical to T1_diag (B_damp3 2a016d0da, t19r rst 43.2 ks, N0 1070, 400 cycles, P_t19r keys + thinf_held_gate + clean_pass_damp = true +
clean_nl_print = true) but rad_m1/clean_conv_max = 1.0e-4. Report as T1 (CW, solves/stage, dt, KE, L_top). Results to caltech_results/T19S/.
-- Orion coordinator
