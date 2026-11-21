# NOTE 2026-11-21je orion -> caltech: TEST on 02-04 — same t19r damp test with damping v3

Damping v1 (B_damp cd66d1ab6) has a weakness (an E<0 iterate can be held negative), fixed in v3. Please run, on the idle TEST node 02-04,
the same test as t19r_damp_test but on v3, in parallel (keep t19r_damp_test on 02-16 running as the v1 comparison):
- build B_damp3 = caltech-runs/stall2v3 2a016d0da (= C3c + clean_pass_damp v3, key clean_pass_damp = true), CUDA 12.9 + g++ 13, smoke gate;
- restart from the same t19r rst (43.2 ks), P_t19r keys + rad_m1/clean_stag_thinf_held_gate = true + rad_m1/clean_pass_damp = true,
  to t = 1.0e5 s or the wall END, rst every 1e4 s (KEEPRST 3), dir prod_fresh_1121b/t19r_damp3_test.
Report alongside t19r_damp_test in the PROD table (CW, solves/stage, dt, KE, L_top/L_in, FATAL). Test only; USER decides production.
-- Orion coordinator
