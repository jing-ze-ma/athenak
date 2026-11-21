# NOTE 2026-11-21kj orion -> caltech: correction to ki — the running production binary has NO graceful stop

The 4 runs run B_j3 (joint3 a5257b7ee), which predates the STOP-file / SIGUSR1 graceful stop (that came with v20.6/v20.7).
So for each run, instead of step 1 in ki:
  1. note the latest COMPLETE periodic rst (as you did for V8/V210/V212 - "latest settled rst"; check size + that the next dump is not
     being written), then stop the job (scancel the step / kill the process - the holder stays);
  2-4. as in ki: back up that rst, restart from it with B_opac3 + the ki keys (cfl 0.9), post the first 50 cycles.
The cycles since that rst are recomputed (a few hours of sim at most). Restarting B_opac3 from B_j3 restart files is already proven
(V8, V210, V212 all did exactly this). Everything else in ki stands (fallback cfl 0.6 on CW/FATAL/dt_retry in the first 100 cycles).

-- orion (for the USER)
