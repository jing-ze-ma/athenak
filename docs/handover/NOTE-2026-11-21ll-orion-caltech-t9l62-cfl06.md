# NOTE 2026-11-21ll orion -> Caltech: USER — t9l62 on head: stop now, relaunch from t9l62.00045 at cfl 0.6

USER (11-21 ~23:20 PST), after your 23:09-23:11 CW notes (head 462210d51 at cfl 0.9: stage fallbacks 1-3 at c33376-33391 and a
backward-Euler dt_retry with "last change inf", 4 so far, ~40 cycles after the restart):

1. Stop the head run now (graceful STOP file).
2. Relaunch the SAME head binary (md5 2144006b) from t9l62.00045 with the same command-line keys as NOTE jk EXCEPT
   `time/cfl_number=0.6` (keep clean_conv_max=1e-5, clean_stag_flim=false, clean_stag_fwb=true, clean_stag_fcone=true,
   clean_stag_fwb_stage=true, output5/dt=21600, VGD_XFE=0).
3. Post the first ~200 cycles as before (CW / stage fallbacks / dt_retry, dt, s/cycle, jh noise).
4. If it still falls back in c33370-33460 at cfl 0.6 (or any "last change inf"), stop and revert to B_j3 (athena_j3_a5257b7ee)
   from t9l62.00045 with the rst's own keys — B_j3 got through this stretch (CW 5 / dt_retry 3 / FATAL 0).

Also: your smoke gate on the head gave pass 1 = 4.934e-06 "NOT converged after 50" (old calibration 4.020e-06). Orion is
re-running the same gate on CPU + A100 with the head and the old binary to confirm the new value is the head's, not the build's;
we post the verdict as soon as it is in.

-- orion (coordinator)
