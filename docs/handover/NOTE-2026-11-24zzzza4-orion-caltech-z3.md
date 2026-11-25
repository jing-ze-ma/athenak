# NOTE 2026-11-24zzzza4 orion -> Caltech: M35ZAMS deck corrected to the Z3 layout

Thanks: the 16-block / tile-32 layout in candidate_1124 was orion's error (an old IA12.2 deck). **The USER's layout is Z3: 4 blocks of 384x64x64, layout_tile 64.** The deck in candidate_1124 is corrected (MD5SUMS now has its line); wall_sym stays on as before (it does not need 16 blocks). Preflight M35ZAMS again (100 cycles, same checks); bsg and m35mms preflights stand as posted. S3G passes (floor_pay 0.959 x, L_top 1.9e-4, zero CW): noted. S5G next, then report.

-- orion (coordinator)
