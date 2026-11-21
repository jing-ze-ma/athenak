# NOTE 2026-11-21ho orion -> caltech: full-resolution star test of the joint fix (after/alongside the C1R arms)

Build B_j3 = caltech-runs/joint3 a5257b7ee (= beamreal + flim v3 + FloorGasCap in the clean path; new keys default off) with CUDA 12.9 +
g++ 13, incremental from the B_br tree, own ccache_cu129. Gate: tiny BSG smoke, healthy = pass 1 ~3.5e-06 and NOT converged (trap = 3.5e-08).
Then 4 arms, one node each (use the 2 nodes free now first, the rest as C1R arms finish), production layout, the WN-GL smoke decks
(*_wngl_smoke) with these keys added IN the deck copies (in-place replacement):
  rad_m1/clean_init_dtfac = 1.0   rad_m1/clean_stag_flim = true   rad_m1/clean_stag_thinf_held = true
  floorcap_rho = 1e-14   floorcap_v = 3e7   floorcap_report = 25   (in the block the binary reads them from — check the never-read list)
  J3_t19 (t19l64_wngl_smoke, tlim 2.0e4 s), J3_t9 (t9l62_wngl_smoke, tlim 1.0e4 s), J3_m35 (m35mms_wngl_smoke, nlim 500), J3_bsg (bsg_wngl_smoke, tlim 3e5 s)
Report every ~20 min: dt, KE_int, L_bot/L_in, L_top/L_in, CW, the CLEAN FLIM / floorcap lines, the never-read list. Results to caltech_results/J3/.
Purpose: confirm at full resolution that the floor-gas runaway (dt collapse at ~1.4 ks in T19/T9, m35 start FATAL) is gone. No production.

-- Orion coordinator
