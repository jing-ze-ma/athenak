# NOTE 2026-11-24zzzze9 orion -> Caltech: bsg's L_z: ledger against residual (read-only, CPU)

Lz_tot of the hst = gas + radiation + clean_lz_ledger. The residual (dLz_tot + int Lz_hflx) is round-off before c6050 and jumps after it. From the c1 bsg log (and T2's, c6061 on): the lines `CLEAN LZB cycle N t2s 3 ledger X` (and `floorcap` lines with the cap's dL_z, floor_pay / floor_audit / FOFC / vceil / fallback / BE counters):
1. Per cycle in c6050-7000 and c8500-9000: d(ledger), the residual of that cycle, nfb change if printed, floorcap dL_z, dfloor_flux_limit and floor counters. Correlate: is the residual (a) equal and opposite to d(ledger) (a booking bug), (b) uncorrelated with it (source outside the clean solve: merge, c2p floors, FOFC, vceil), (c) at cycles with stage fallback / BE redo / DT_RETRY only? Give the correlation coefficients and the 10 largest residual cycles with every column.
2. State what prints are missing. Table in caltech_results/BSG_LZ/lzb.txt. Low priority, one core.

-- orion (coordinator), asked by KEYTRIM
