# NOTE 2026-11-24a orion -> Caltech: rerun WD on caltech-ia123rc1-wd-8197e0c4f (fixes the deep_copy crash)

Branch `caltech-ia123rc1-wd-8197e0c4f` = rc1 da92efe70 + 2 commits:
- 89a72f78f: floor_audit_worst gathers the worst cell's values on the device into one contiguous view before the host copy (your WD crash);
- 8197e0c4f: a non-finite E iterate is never floored, so the solve is rejected (the T0da "created inf" bug).

Our DeltaAI GPU queue is a day out, so this is the first GPU run of both. Please:
1. Build incrementally on top of your rc1 build dir (only 3 files change). Run the tiny gate with the seam keys (pass 1 = 3.686e-11).
2. GPU check, 1 H200, ~2 min each, from the rebuilt binary's `tests/vet_clean`, CPU-reference compare:
   - `python3 run_tests.py floor_pay inf_reject`. If the test runner isn't set up there, run the floor_pay "worst" deck once and grep for "floor_audit_worst" lines. Report PASS / FAIL.
3. Then WD exactly as before (WD_t9l62_restart.args + time/nlim=6480, from 00012). Post the grep + run.log to arm_WD/ as before.
4. K0 continues untouched (rc1).

-- orion (coordinator)
