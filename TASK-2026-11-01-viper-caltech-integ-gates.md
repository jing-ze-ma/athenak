# TASK 2026-11-01 viper -> Caltech (USER): gate the consolidated defaults merge, fork/integ-test-1101 @ e32a104ea

Follow-up to the heads-up NOTE f74644f05. Please build **e32a104eaa17eb552090ec0484f6e0a7f7909af4** (H200, fmad off) and run (a) and (b) in parallel. Your PASS on both, together with our CPU gates, triggers the push of THIS commit to fin51-r-fix (fast-forward). That push is your r1m launch build.

## Contents (on fin51-r-fix 8cf6db628)
- the tracers (8ccae663) + `ambipolar = eos` default;
- cs_seam_fc_nolimit default;
- the production defaults.

**USER change: use_rkg_sts is NOT defaulted.** rkg_lean = true and rkg_dt_safety = 0.5 become defaults only when STS is on (your decks set use_rkg_sts = true, so r1m gets both).

All defaults apply on fresh starts only (recipe::NewDefaults). Restarts without the key keep the old behaviour, and explicit keys win.

## (a) Key-free bitwise vs 8cf6db628
- A developed r1m restart (e.g. keep/off_0p50.rst) with your current r1m keys, ~150 cycles. Run your 8cf6db628 binary vs e32a104e.
- PASS:
  - rst bytes after `<par_end>` identical;
  - dhj.mhd.hst identical;
  - dhj.user.hst floor columns within 1e-12 (GPU atomics).

## (b) Fresh r1m start with no special keys, ~0.1 rot
- Keys: planet, grid, EOS/ck/IC, outputs, ohmic_resistivity eos, max_eta 1e15, valfven_max_cgs 2e6, use_rkg_sts true, problem/roche true.
- Remove from the deck: xe_tracer, floor_tracer, valfven_noinfall_fac, ambipolar, rkg_lean, rkg_dt_safety, top_bc, roche_wb_tan, cs_seam_fc_nolimit, cs_seam_flux.
- PASS:
  - clean run (no FATAL or nan; event counters like p_f4);
  - the start-up defaults block lists as set: xe_tracer, floor_tracer, valfven_noinfall_fac 4, ambipolar eos, cs_seam_fc_nolimit, rkg_lean true, rkg_dt_safety 0.5, top_bc outflow1, roche_wb_tan true (plus the existing am_transport / dt_retry).
- Please paste the block into your NOTE.
- Optional: the same deck with all of those keys explicit, ~100 cycles, should be bitwise to the no-key run (rst data).
