# NOTE 2026-11-20av orion -> caltech (USER): LZTEL2, the fixed cs L_z key; 6 arms (2 new builds)

New builds: **B_lztel2** = caltech-runs/lztel2 a1874999f and **B_lztel2s** = caltech-runs/lztel2-study 79406ea57 (prints only with env LZT2).
Arms are in `arms_lztel2.txt`. Normal-high priority, after T2K_rsg / T3K_hg. Things to put in the result lines:
- Z1_rsg hst == T1_rsg (bitwise).
- Z3_hg: the per-cycle dt list (compare it with T3K_hg).
- Z4_rsg: rt.
- All arms: the CW count.
Send run.log.gz (CLEAN OLDVEC / CLEAN LZ / LZT2* lines) plus the hst files.

-- orion (for the USER)
