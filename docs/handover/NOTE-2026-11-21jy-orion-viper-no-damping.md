# NOTE 2026-11-21jy orion -> viper: USER decision — no pass damping from now on

USER (2026-11-21 ~11:40 PST): do not use rad_m1/clean_pass_damp from now on, in any run (production, successors, tests).
- Never name clean_pass_damp (the code default 0 = off). Remove it from any deck / key file / test spec you hold.
- Evidence: at clean_conv_max 1e-4, T9 / m35 / T19 / HG run 0 CW without it and with fewer passes (V8 damp-0 twins: 9-22 % faster than
  production); damp 1 costs +55-100 % solves there. The RSG stall it used to cure is being fixed at the root
  (inexact inner solve + xdonor donor/cap flips: candidate keys clean_conv_ew_max 0.01 + clean_stag_xdonor_freeze 8, under test on viper).
- Pending damp tests are cancelled (the V7a / V7c T19 damp-0 reruns are no longer needed either).


-- orion (for the USER)
