# NOTE 2026-11-20ae orion -> caltech (USER): WITHDRAW the LZTEL arms; run the LZCS tel gate instead (new build lzcs-tel)

1. **Cancel all LZTEL arms** (NOTE 2026-11-20ac, `arms_lztel.txt`, build caltech-runs/lztel 5f18abdc7) that have not started.
   It was a duplicate. If one is already running, let it finish and send it.
2. New build **B_lzcstel** = caltech-runs/lzcs-tel, commit 61a578136 (v20j plus clean_cs_lz_tel, default off, plus a ledger
   booking). Arms: `docs/handover/caltech_runs_1121/arms_lzcs_tel.txt`. Normal priority, after the v20j gate and the SPG A/B.
   - T1_rsg: key off, 6 cycles. Its hst must equal V20JR_rsg (v20j, 4 cycles) over the common cycles; please add that xcheck.
   - T2_rsg: key on, 300 cycles, hst every cycle.
   - T3_hg: key on, 50 cycles.
   - T4_rsg: key on, rt 3 + 3, bitwise.
   Send back: run.log.gz (CLEAN LZ lines, CW) plus the hst files. If the key FATALs as "not found", send the line and skip.

-- orion (for the USER)
