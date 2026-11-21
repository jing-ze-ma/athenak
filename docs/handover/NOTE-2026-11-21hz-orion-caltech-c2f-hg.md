# NOTE 2026-11-21hz orion -> caltech: 2 HG arms with the full joint fix (after J3D/C2B free a node)

On B_j3 (joint3 a5257b7ee), hegiant.00139, 100 cycles, keys in <rad_m1> of deck copies:
  C2F_hg_c1 : C1R_hg_c1 keys + rad_m1/clean_stag_flim = true + rad_m1/floorcap_rho = 1e-14 + rad_m1/floorcap_v = 3e7 + rad_m1/floorcap_report = 25
  C2F_hg_v7 : C1R_hg_v7 keys + the same flim/floorcap keys (control)
Report passes/cycle, s/cycle, CW, min dt, never-read list, FLIM/floorcap lines. Results to caltech_results/C1R/.
Reason: the HG stall cells are thin floor cells filled from E 1e-4 to ~55 in a step; the C1R arms ran without flim/floorcap.

-- Orion coordinator
