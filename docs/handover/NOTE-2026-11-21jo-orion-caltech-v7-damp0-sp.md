# NOTE 2026-11-21jo orion -> caltech: rerun the sp damp-0 arms V7a and V7c (USER), keep V7a2/V7c2

USER: still run the damp-0 tests for the spherical-polar T19 runs that were stopped by NOTE jn. (The no-damp-0 rule is for the
cubed-sphere RSG; sp T19 damp 0 is still wanted.)

1. Rerun V7a (t19r restart t19l64.00001, clean_pass_damp 0) and V7c (t19l64 restart t19l64.00001, J4 keys from the rst, clean_pass_damp 0)
   exactly as in NOTE ji: B_v207 md5 86d31021, 600 cycles, rst off, + clean_conv_max 1e-4 + clean_stag_thinf_held_gate.
2. Nodes: any idle NON-production node now (e.g. those freed by bsg / i1_bsg / t19l64 / t19r if not in use by an Orion test);
   otherwise the first of 02-16 / 02-11 when V7b / V7d finish. One GPU job per node; GPU-free check as fixed at 10:29.
3. Keep V7a2 / V7c2 (damp 2) and V7b / V7d (damp 1) running. Results to caltech_results/V207/V7a and V7c (overwrite the stopped ones).
V8 stays on hold.

-- orion (for the USER)
