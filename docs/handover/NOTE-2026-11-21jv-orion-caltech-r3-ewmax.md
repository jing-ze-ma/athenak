# NOTE 2026-11-21jv orion -> caltech: replace R3 (wb_local alone stalls on viper) by E1 = tighter linear-solve forcing

viper result: clean_wb_local + damp 0 still stalls on rsg.01510 (CW at the first window). So R3 as specified is not needed.
1. Branch caltech-test-pos5-wblocal has moved to c2eb74fac (= 6653fdef4 + new key rad_m1/clean_conv_ew_max, default 0.1 = bitwise).
   Rebuild from it (CUDA 12.9 + g++13, tiny BSG smoke gate).
2. Run E1 instead of R3: same RSG arm as R2 (your rsg restart, 150 cycles, rst off, damp 0) + rad_m1/clean_conv_ew_max = 0.01.
   (Keep R2 running: it is the control on your restart.)
3. Priority unchanged otherwise (ju): R2 > V8a0/V8b0 > E1 > V8c0/V8d0 > V7a/V7c reruns.

-- orion (for the USER)
