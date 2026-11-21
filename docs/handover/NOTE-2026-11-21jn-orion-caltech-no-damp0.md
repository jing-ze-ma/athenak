# NOTE 2026-11-21jn orion -> caltech: stop the damp-0 V7 arms, rerun them with clean_pass_damp = 2 (USER)

USER: clean_pass_damp = 0 does not work (viper RSG rsg.01510 at conv_max 1e-4: 9 CW, 200-solve stalls), so no more damp-0 tests.

1. Stop V7a (t19r, pass_damp 0) and V7c (t19l64, pass_damp 0) now. These are our own test jobs; touch nothing else.
2. On the same nodes, start the same arms with `rad_m1/clean_pass_damp = 2` and every other key unchanged (B_v207 md5 86d31021,
   600 cycles, rst off, clean_conv_max 1e-4, clean_stag_thinf_held_gate): V7a2 = t19r damp 2 on 02-04, V7c2 = t19l64 damp 2 on 01-11.
3. Keep V7b and V7d (damp 1) running. Results go to caltech_results/V207/ as before. V8 stays on hold.
The question is now only damp N = 1 vs N = 2.

-- orion (for the USER)
