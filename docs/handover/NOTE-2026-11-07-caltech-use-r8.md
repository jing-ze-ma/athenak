# NOTE 2026-11-07 Caltech -> viper (USER: "tell viper to use r8 settings")

Use the **r8** key set for the M1 runs on Raven: HG, RSG, R136a1 B and the BSG tests. The USER wants similar settings in all runs.

r8 = seamE_1106 arm r8 = r7 + 3 opacity-lag passes (full command line: `$J/seamE_1106/runs/r8/arm.info`). Add these on the command line or in the deck, on top of each run's current production keys:
```
rad_m1/cs_cross_mono=true
rad_m1/implicit_opac_newton=false
mesh/cs_seam_rho_guard=4
rad_m1/implicit_predictor_tcap=1.1
rad_m1/time2_k1_redist=0.1
rad_m1/implicit_fb_local=true
rad_m1/implicit_fb_local_rounds=4
rad_m1/implicit_opac_lag_passes=3
```
- Name every key explicitly at each restart; do not rely on defaults.
- R136a1 keeps `rad_m1/time2_adm_poslim=true`.
- Full tensor stays as already agreed: `vet_gd_tangential` on sp, cs-tan on cs, if its STEP 0 check passes.

**HG result (from HG rst 139, 6 H200)**
- r7 vs the old keys: BE redos 1 vs 153, NOT POSITIVE 1 vs 184, wall −32 %, budget closed.
- r8 = r7 with wall +3 %, 2 redos.
- The opacity-lag error falls ~9x (1.7e-4 -> 1.9e-5 at dt 105). This lag was the cause of the layout excess.

**Known gaps; check them in the first segments**
1. **fb_local:** no restart or layout gate yet. Run one short restart pair and one layout pair with r8 before the long segment. Use bin_double and GPU; seamE REPORT lists the recipe. The rst_state_ext fix in prod-1107 makes hesdirk2 restarts bitwise now.
2. **RSG:** with fb_local, the floor count rose in rsgdiag_1107 (unexplained). Watch NOT POSITIVE, NOT ADMISSIBLE and the floor counts against the RSG r7 test: per cycle, NOT POSITIVE 0.21 and NOT ADMISSIBLE 0.07. If floors grow, report back before continuing.
3. `time2_k1_redist` moves ~15 % L gross through redistribution. It is part of r8 per the USER; keep printing the k1_redist report periodically.
4. r8 was run on HG only, about 1 d sim. The RSG, R136a1 and BSG decks have no r8 run yet, so treat the first segment of each as the test and post the counters.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
