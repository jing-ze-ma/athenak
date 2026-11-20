# NOTE 2026-11-20do orion -> caltech (USER): v20.3 GPU gate (5 short arms)

New build **B_v203** = caltech-runs/v203, commit e688cb7bc65c6650d8221454f062f1369093bdf7 (v20.2 + the EMA-wb key + an FTIE fix + the LZTEL3 L_z fix; the USER adopted clean_cs_lz_tel=true for production). Build it when a slot frees.
The arms are in `arms_v203.txt` (4 cycles each, 1 node, the V202 recipe and restarts), **after the current long arms**:
1. V203_rsg and V203_hg (rt, floor_pay on): hst should equal V202F_*/a bitwise.
2. V203L_rsg and V203L_hg (rt, + clean_cs_lz_tel=true): pass means rt 5/5 and 0 CW.
3. V203LP_rsg (poison 2): hst should equal V203L_rsg.
Do **not** run the EW (ema_wb) pairs on this binary; they stay on B_emawb.

-- orion (for the USER)
