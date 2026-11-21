# NOTE 2026-11-21jb orion -> caltech: TEST on 02-04 — D_t19 = C3c + clean_pass_damp on the T19 restart

On viper, rad_m1/clean_pass_damp (damps the period-2/3 pass cycle of the stage solve; fixed point unchanged) cured the RSG stall:
137 cycles, 0 CW, 0 dt_retry, mean 5.7 solves/stage, dt flat — where C1/C3c stall. Please test it on T19 on the TEST node 02-04
(not production):
- Build B_damp = caltech-runs/stall2 cd66d1ab6 (= C3c 411bf3991 + clean_pass_damp, default off) with CUDA 12.9 + g++ 13, incremental, smoke gate.
- D_t19: identical to C3c_t19 (restart from t19l64.00001, N0 1127, J4 keys + clean_stag_thinf_held_gate = true + clean_nl_print = true,
  450 cycles, rst off) plus rad_m1/clean_pass_damp = true.
- Report vs C3c_t19 and the stopped production over the same cycles: CW, mean/max solves per stage, dt, KE, L_top/L_in. Results to caltech_results/C3H/D_t19.
The USER is asleep: test node only, no production changes.
-- Orion coordinator
