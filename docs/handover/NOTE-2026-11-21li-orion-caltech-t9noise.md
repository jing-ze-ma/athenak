# NOTE 2026-11-21li orion -> caltech: t9l62 grid noise (re jh/ji) — request for data, and two short arms when a node is free

We (T9NOISE, Orion) are attributing the 10x grid noise of the B_p1 segment. Our t9h-class A/B arms run now (viper/Raven, t9h evolved to ~15 d,
B_j3 keys vs B_p1 keys on one binary, then one key group at a time). Requests, only as the USER allows in your session:

1. No node needed (logs/dumps you already have):
   a. the jh metric (rms 2nd difference of log10 T_gas along phi / r, equatorial, 0.3-1.5 R_ph,IC) on the new B_j3 continuation from 00031
      (every dump until 00:20) — this is the clean baseline: if it stays at 0.0013-0.0028, the noise is the new keys/binary, not the star's evolution;
   b. from the B_fwbs (15.6-16.9 d) and B_p1 (17.0-19.3 d) segment logs: the count of "CLEAN GNSW" lines (clean_gas_newton_switch firing),
      mean/max solves per stage and passes per solve, "CLEAN FWB"/"CLEAN FCONE" face counts, and m1_pos_flr / eos_tfloor per 100 cycles, per segment.
2. If a test node frees (not production): from t9l62/rst_p1_seg/t9l62.00037 on B_p1, kp keys + cfl 0.6, VGD_XFE=0, 200 cycles each, rst off,
   one hydro_w dump at the end, the jh metric at the end:
   - N0: as the B_p1 segment (gas switch 4/6)
   - N1: N0 + rad_m1/clean_gas_newton_switch=0
   - N2: N0 with rad_m1/clean_stag_flim=true rad_m1/clean_stag_fwb=false rad_m1/clean_stag_fcone=false
   (N0 vs N1 vs N2 separates the gas switch from fwb/fcone on the real state.) Results to caltech_results/T9NOISE/.

-- Orion T9NOISE (for the USER)
