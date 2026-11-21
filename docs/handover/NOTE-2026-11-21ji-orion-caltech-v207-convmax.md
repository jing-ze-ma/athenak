# NOTE 2026-11-21ji orion -> caltech: v20.7 pushed (USER); clean_conv_max = 1e-4 for production (USER); T19 tests on the idle nodes

1. USER: origin/clean-stagger2-latest = af9e8d420 (code 5ba73ecab) = v20.7 (v20.5 + graceful stop + refined clean_stag_flim + C3c +
   colflux + xdonor_freeze + L_z-ledger/floorcap rollback fix on BE stage fallbacks + clean_pass_damp integer key; all new keys default off).
   Please build B_v207 = caltech-runs/v207 5ba73ecab (CUDA 12.9 + g++ 13, smoke gate) for future (re)starts.
2. USER: production decks use rad_m1/clean_conv_max = 1.0e-4 from their next (re)start (1e-5 sits below the ~2e-5 pass noise floor and
   caused the T19 stage stalls; T4_conv with 1e-4 had 0 CW). Running production is not restarted for this; apply it at the successor
   restarts after the 17:25 END (successors still need the USER's confirmation in your session).
3. TESTS on the idle nodes (not production), B_v207, ~600 cycles each, rst off, all with clean_conv_max = 1e-4:
   - 02-04 V7a: t19r rst (43.2 ks), P_t19r keys + thinf_held_gate = true, clean_pass_damp = 0
   - 02-16 V7b: same + clean_pass_damp = 1
   - 01-11 V7c: t19l64.00001, J4 keys + thinf_held_gate = true, clean_pass_damp = 0
   - 02-11 V7d: V7c + clean_pass_damp = 1
   Report: CW, mean/max solves per stage, dt, KE, L_top/L_in. Results to caltech_results/V207/. These decide whether damping is needed
   once the tolerance is 1e-4.
-- Orion coordinator (for the USER)
