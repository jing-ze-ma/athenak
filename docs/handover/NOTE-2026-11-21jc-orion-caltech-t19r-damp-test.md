# NOTE 2026-11-21jc orion -> caltech: TEST on the idle node 02-16 — t19r restart with C3c + clean_pass_damp, past the FATAL point

t19r died (FATAL at c4429, t 7.95e4 s). As a TEST (not production; the USER decides production when awake), please run on 02-16:
- binary B_damp (caltech-runs/stall2 cd66d1ab6, already built and gated on 02-04; copy it),
- restart from t19r's newest rst (t19r.00001 at t 43.2 ks, or whatever is newest before the crash), with the P_t19r production keys +
  rad_m1/clean_stag_thinf_held_gate = true + rad_m1/clean_pass_damp = true,
- run until t = 1.0e5 s (past the FATAL at 7.95e4) or the wall END, rst every 1e4 s (KEEPRST 3), into prod_fresh_1121b/t19r_damp_test.
Report every ~20 min in the PROD table (CW, solves/stage, dt, KE, L_top/L_in, any FATAL) next to t19r's history.
Node 02-04: still free for a test — keep idle unless we send one.
-- Orion coordinator
