# NOTE 2026-11-24e orion -> Caltech: IA12.3-rc2 = caltech-ia123rc2-5b493acfa: GPU key tests + fresh bsg / t9l62 (USER-approved production start if clean)

Decks in this repo: `docs/handover/ia123rc2_decks/` (replace merge123/caltech_ia123rc2/ in the text below with that path).

**Also:** your K0 already measured the new PCR kernel on H200: rc1 3.600 vs IA12.2 3.778 s/cyc (-4.7 %), bitwise. Thanks.

**Production start (b):** the USER OK'd bsg / t9l62 production on rc2 if (a) passes and the 200-cycle arms are clean (CW 0, no FATAL). For t9l62 the outburst fix (T9FLOOR's FOFC arms WD0 / WD2, next NOTE) comes later, via a restart long before t ~4.8e5.

## Code
- Branch ia123-rc2, tip 5b493acfa. Pushed as caltech-ia123rc2-5b493acfa.
- It is rc1 da92efe70 plus:
  - VCOV: the cs ceiling FOFC phi_hat fix, the cs floorcap metric norm, and the ceiling_metric test;
  - WALLSYM: CleanSeamImage, which brings the theta-seam residual to round-off, and the z-turned sets;
  - LAYC2: a2b44d59e, the sdirk2 u0 ghost fix for dfloor_flux_limit;
  - T9FLOOR: 89a72f78f, the device gather for floor_audit_worst (WD's H200 crash), and 6caed05b8, a non-finite E is rejected and never floored;
  - HALF12: the FwbProject counters fix and the clean_stag_mg_levels key;
  - visitor skips for the cache / one-shot members cln_pcrl, cln_sg*, cln_dbginf.
- Build with CUDA 12.9 + g++13, twice from the SAME commit:
  - PROBLEM=he_star_vet, for the stars;
  - PROBLEM=built_in_pgens, for the test suite of arm (a).
- Run the tiny gate first. Deck: merge123/caltech_ia123rc2/IA123_gate_bsg_smoke.athinput, plus `mesh/nx2=32 mesh/nx3=32 meshblock/nx2=32 meshblock/nx3=32 time/nlim=2`. It must give pass 1 3.686e-11.

## Arms
- **(a) GPU key tests on 1 H200** (the built_in_pgens binary, from the branch's tests/vet_clean):
  - command: `python3 tests/vet_clean/run_tests.py --exe <bin> --only floor_pay,inf_reject,wedge_sym,layout_tile,ceiling_metric,conservation_lz,restart_sp,dt_retry --jobs 4 --threads 1`;
  - post its PASS / FAIL lines;
  - MI300A reference (viper apudev, the same tests on rc2 5b493acfa built_in_pgens HIP): 9/9 PASS. floor_pay worst_lines 199 (bitwise), inf_reject rejected 1, wedge_sym fth 5.5e-15, layout_tile dfl 0.
- **(b) fresh bsg and t9l62**: decks merge123/caltech_ia123rc2/{tile,notile}/IA123_{bsg,t9l62}.athinput (tile is the production form), 200 cycles each, 4 H200.
  - The decks name, in release names:
    - vet_gd_wall_sym + vet_gd_theta_seam_sym;
    - floor_from_ekin;
    - diff + max_stage2 (already there);
    - clos1_budget / bud_maxnorm / vel_settled named OFF (HALF12 / RWFIX).
  - Post: s/cyc; CLEAN WARNING / FATAL / DT_RETRY counts; L_top / L_in at cycle 200; the hst md5.
  - If both are clean, the USER's OK makes these the production start: same decks, tlim as before.
- **(c) R136: NOT possible on rc2 yet.** R136's set (the default geodesic set, z-rotating) is not mirror-symmetric about the equator. vet_gd_theta_seam_sym then FATALs: "the set is not mirror-symmetric about the wedge mid-plane" (Orion prod row, R2S). Without the key, the theta-image FATAL fires. R136 needs a GL / equator-symmetric set first; that is a USER decision.

-- orion (coordinator)
