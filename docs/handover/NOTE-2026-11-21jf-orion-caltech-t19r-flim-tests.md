# NOTE 2026-11-21jf orion -> caltech: T19 stall tests — replace the two damp tests with three short ones (TEST nodes only)

The t19r damp tests answered their question (v1 7 CW, v3 10 CW). Please stop t19r_damp_test (02-16) and t19r_damp3_test (02-04)
(test runs; keep their dirs), and run three short tests, each from the same t19r rst (43.2 ks), P_t19r keys + thinf_held_gate = true,
~400 cycles (past the first CW at c~1102), rst off:
- 02-04, T1_diag: B_damp3 (2a016d0da) + rad_m1/clean_pass_damp = true + rad_m1/clean_nl_print = true. Push the last 30 CLEAN pass /
  XDONOR / XT / FLIM lines of the FIRST failing stage, plus all CW lines.
- 02-16, T2_v206: build B_v206d = caltech-runs/v206damp 6241cce78 (= v20.6: refined clean_stag_flim, C3c, colflux, stop; + clean_pass_damp
  integer key), CUDA 12.9 + g++ 13, smoke gate; keys + rad_m1/clean_pass_damp = 1 (the deck keeps clean_stag_flim = true -> refined flim).
- 02-11 (idle since i1_bsg ended), T3_v206: B_v206d with rad_m1/clean_pass_damp = 0 (refined flim alone, control).
Report per test: CW count, mean/max solves per stage, dt, KE, L_top/L_in, FATAL. Results to caltech_results/T19S/. Test only; the USER
decides production.
-- Orion coordinator
