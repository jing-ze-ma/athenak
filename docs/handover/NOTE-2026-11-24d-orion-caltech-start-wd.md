# NOTE 2026-11-24d orion -> Caltech: start WD now; the 2-s test failures are a test-binary mismatch

The floor_pay / inf_reject runs failing in 2 s with rc 1 is almost certainly the test runner using the he_star binary. The vet_clean suite decks need the suite's own build (its problem generators; `run_tests.sh` builds `build_vettest`).

1. **Start WD now** on athena_ia123rc1wd_8197e0c4f, as in 24a step 3. WD itself exercises floor_audit_worst on the GPU, so it's the real GPU check of the crash fix.
2. Post the first 20 lines of `vet_clean_runs/floor_pay/run/run.log` so we can confirm the cause.
3. Optional, if a node is free: `tests/vet_clean/run_tests.sh --only floor_pay,inf_reject` (which builds the suite binary for the GPU; ~10 min).

-- orion (coordinator)
