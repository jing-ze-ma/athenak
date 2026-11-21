# NOTE 2026-11-21in orion -> caltech: B_v205 gate is a PASS; V5 order; I1_t19 (rebuilt T19) on the next test node

1. B_v205 gate: pass 1 = 2.520e-06 and NOT a trap (trap = 3.5e-08 + converged in 4). B_v205 turns clean_stag_flim on by default, which changes
   the smoke numbers. The 1-TU rebuild is correct: a5257b7ee..b3fd494e8 differs only in src/rad_m1/rad_m1.cpp. Please use B_v205.
2. Order on the 2 test nodes: C2X_hg_c1 / C2X_hg_v7 (NOTE il) first; then I1_t19 (below); then V5a/V5b (NOTE ik) and V5c after all
   (no other GPU can start it soon).
3. I1_t19: rebuilt T19L6.4 IC (viper ic_star; Jiang init_80 with noise rows removed + mass fix), files docs/handover/caltech_runs_1121/icrel/t19l64.athinput +
   ic_t19l64_thin4m_dfl19.txt.gz (gunzip; set problem/he_ic_file). Deck has he_ic_rad_atm_relax = 0.25 (USER-approved; the atmosphere march
   needs it), clean_init_dtfac 1e6, flim, thinf_held, floorcap. On B_j3 (or B_v205), production layout, one free TEST node, tlim 2.0e4 s;
   progress lines as for I1_t9.

-- Orion coordinator
