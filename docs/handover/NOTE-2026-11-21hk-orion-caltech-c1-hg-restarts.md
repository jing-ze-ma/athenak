# NOTE 2026-11-21hk orion -> caltech: evolved HG/RSG restart test of the new beam fix C1 (one build + 6 short arms)

1. Build B_br = caltech-runs/beamreal e78275364 (= v20.4 207fab564 + clean_init_dtfac default 1 + clean_stag_flim + clean_stag_thinf_held
   "C1"; all new keys default off) with CUDA 12.9 + g++ 13 (incremental from the B_v204n tree, own ccache_cu129), smoke-gate (tiny BSG pass 1 = 4.020e-06).
2. On B_br, 100 cycles each, same restarts and decks as CM4_hg / V203L_rsg (hegiant.00139, rsg.01295), report passes/cycle, s/cycle, CW, min dt:
   - C1R_hg_ref : v20.3 keys (THINF off, maxit 200, conv_max 1e-4)
   - C1R_hg_v7  : as CM4_hg (THINF on, maxit 800, conv_max 1e-4)
   - C1R_hg_c1  : as CM4_hg + rad_m1/clean_stag_thinf_held = true
   - C1R_rsg_ref / C1R_rsg_v7 / C1R_rsg_c1 : the same three on the RSG restart.
   Keys in deck copies (in-place replacement). Test nodes only. Results to caltech_results/C1R/.
Purpose: v7 broke restarted HG (CM4_hg: 139 passes/cycle, 8 CW, dt to 6.8); C1 is designed to remove that.

-- Orion coordinator
