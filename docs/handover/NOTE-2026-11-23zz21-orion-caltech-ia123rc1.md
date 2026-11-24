# NOTE 2026-11-23zz21 orion -> Caltech: IA12.3-rc1 = caltech-ia123rc1-da92efe70. ONE binary for all new tests

## Binary
- Branch `caltech-ia123rc1-da92efe70`: the cleaned release trunk (m1 -> vet rename) + IA12.3, incl. T9FLOOR's `hydro/floor_audit_worst`.
- Build with CUDA 12.9 + g++ 13, `PROBLEM=he_star_vet` (`he_star_m1` is an alias). Tiny BSG gate first: pass 1 = 3.686e-11.
- Deck block is `<rad_vet>`; old `<rad_m1>` decks and rsts run through the aliases.
- Decks: `docs/handover/ia123rc1_decks/{tile,notile}/IA123_*.athinput`, `IA123_t9l62_restart.args`, `WD_t9l62_restart.args`.

## Arms now (in parallel where nodes allow)
1. **WD** (t9l62 diagnostic, decides the outburst fix), from t9l62.00012.rst:
   - `WD_t9l62_restart.args` (every key named + rad_vet/clean_wb_mcap=true, rad_vet/clean_stag_fcone_emax=true, hydro/floor_audit_worst=true);
   - `time/nlim=6480 job/final_restart=false`, rst/bin dt 1e30.
   - Check the log echo for floor_from_ekin = true, floor_audit_worst = true, clean_gas_secant = max_stage2.
   - Post `grep -E "floor_audit|DT_RETRY|CLEAN WARNING|FALLBACK"` of log.0 + run.log to caltech_results/IA121_T9L62/arm_WD/.
2. **K0** (rc1 key-off check + the H200 speed of the new PCR kernel):
   - rc1 binary with the IA12.2 deck `ia122_decks/notile/IA122_m35mms.athinput` + `rad_m1/vet_gd_theta_seam_sym=true rad_m1/vet_gd_wall_sym=true`, from the latest m35mms production rst, 200 cycles, outputs off;
   - compare against the IA12.2 binary on the same deck + keys + rst: hst bitwise expected, and post s/cyc for both.
3. **No K1 yet.** Two speed keys (clos1_budget, bud_maxnorm) fail a free-streaming accuracy test and are being fixed. K1 follows when they're clean.
4. Fresh bsg / t9l62 production stays HELD (USER), until the θ-seam fix lands.
5. m35mms production on IA12.2 B: unchanged.

-- orion (coordinator)
