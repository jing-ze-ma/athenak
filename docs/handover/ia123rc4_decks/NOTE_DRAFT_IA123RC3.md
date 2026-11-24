# NOTE draft (MERGE123 -> COORD to post): IA12.3-rc4 on Caltech H200, one source, GPU key check + production start

## Code
- Branch ia123-rc3, tip af836f5a3. COORD pushes it as caltech-ia123rc3-af836f5a3.
- It is rc2 5b493acfa + RWFIX (the budget tolerance floor is clean_tol_nl when a budget key is named). rc2 = rc1 da92efe70 plus:
  - VCOV: the cs ceiling FOFC phi_hat fix, the cs floorcap metric norm, and the ceiling_metric test;
  - WALLSYM: CleanSeamImage, which brings the theta-seam residual to round-off, and the z-turned sets;
  - LAYC2: a2b44d59e, the sdirk2 u0 ghost fix for dfloor_flux_limit;
  - T9FLOOR: 89a72f78f, the device gather for floor_audit_worst (WD's H200 crash), and 6caed05b8, a non-finite E is rejected and never floored;
  - HALF12: the FwbProject counters fix and the clean_stag_mg_levels key;
  - visitor skips for the cache / one-shot members cln_pcrl, cln_sg*, cln_dbginf.
- Build with CUDA 12.9 + g++13, twice from the SAME commit:
  - PROBLEM=he_star_vet, for the stars;
  - PROBLEM=built_in_pgens, for the test suite of arm (a).
- Run the tiny gate first. Deck: merge123/caltech_ia123rc3/IA123_gate_bsg_smoke.athinput, plus `mesh/nx2=32 mesh/nx3=32 meshblock/nx2=32 meshblock/nx3=32 time/nlim=2`. It must give pass 1 3.686e-11.

## Arms
- **(a) GPU key tests on 1 H200** (the built_in_pgens binary, from the branch's tests/vet_clean):
  - command: `python3 tests/vet_clean/run_tests.py --exe <bin> --only floor_pay,inf_reject,wedge_sym,layout_tile,ceiling_metric,conservation_lz,restart_sp,dt_retry --jobs 4 --threads 1`;
  - post its PASS / FAIL lines;
  - MI300A reference (viper apudev, the same tests + radwave on rc3 (8ae1cce74 = af836f5a3 + inert d45a71cd2) built_in_pgens HIP): 10/10 PASS. floor_pay worst_lines 199 (bitwise), inf_reject rejected 1, wedge_sym fth 5.5e-15, layout_tile dfl 0.
- **(b) fresh bsg and t9l62**: decks merge123/caltech_ia123rc3/{tile,notile}/IA123_{bsg,t9l62}.athinput (tile is the production form), 200 cycles each, 4 H200.
  - The decks name, in release names:
    - vet_gd_wall_sym + vet_gd_theta_seam_sym;
    - floor_from_ekin;
    - diff + max_stage2 (already there);
    - clean_clos1_budget = true + clean_bud_maxnorm = true (HALF12; radwave fixed by RWFIX); clean_vel_settled named false.
  - Post: s/cyc; CLEAN WARNING / FATAL / DT_RETRY counts; L_top / L_in at cycle 200; the hst md5.
  - If both are clean, the USER's OK makes these the production start: same decks, tlim as before.
- **(c) R136: NOT possible on rc2 yet.** R136's set (the default geodesic set, z-rotating) is not mirror-symmetric about the equator. vet_gd_theta_seam_sym then FATALs: "the set is not mirror-symmetric about the wedge mid-plane" (Orion prod row, R2S). Without the key, the theta-image FATAL fires. R136 needs a GL / equator-symmetric set first; that is a USER decision.

-- orion (coordinator)
