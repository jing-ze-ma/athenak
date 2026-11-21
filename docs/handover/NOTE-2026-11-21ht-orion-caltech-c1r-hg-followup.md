# NOTE 2026-11-21ht orion -> caltech: C1R HG follow-up (3 arms on B_br, parallel, alongside J3)

C1R_hg_c1 collapsed (dt 106 -> 3.40 at cycle 625004, 13 CW, all t2s 2 stage solves stalling near 1e-5 under strict retry).
1. Please commit the FULL a_run.log of C1R_hg_c1 for cycles 624983-624985 and 625000-625006 (or the whole log gz) to caltech_results/C1R/.
2. On B_br, hegiant.00139, keys as C1R_hg_c1 plus (in <rad_m1> of the deck copy):
   - C2D_hg_c1diag : + rad_m1/clean_nl_print = true + rad_m1/clean_diag_cells = true, 45 cycles (per-pass residuals + worst cells)
   - C2B_hg_c1     : + rad_m1/clean_retry_budget = true, 100 cycles
   - C2B_hg_v7     : C1R_hg_v7 keys + rad_m1/clean_retry_budget = true, 100 cycles (control)
   Report passes/cycle, s/cycle, CW, min dt, and the never-read list (to confirm the keys were read). Free nodes only; J3 has priority.
Results to caltech_results/C1R/.

-- Orion coordinator
