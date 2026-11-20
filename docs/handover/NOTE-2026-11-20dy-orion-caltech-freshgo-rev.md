# NOTE 2026-11-20dy orion -> caltech (USER): FRESHGO revision (f1b6421eb): guarded armg.sh, 4 THINF smoke arms, 2 timing arms

This supersedes the arm list in NOTE dx; use the files as now on bsg-files-1009 in `docs/handover/fresh_prod_1120/freshgo/`.
1. **armg.sh md5 fcc1860b.** It now refuses an empty TAG / BIN / DECK / NC / OUT, a '/' in TAG, and a shallow OUT before the rm -rf. Use this version.
2. **Add 4 THINF smoke arms** on the **same B_v203**, with the keys inside the deck. Add **no extra command-line keys**; an unknown key FATALs.
   - T_t9l62 (decks/t9l62_smoke_thinf.athinput), T_bsg, T_t19l64, T_m35mms: 120 cycles each (a = 120 straight; b = 60 + rst + 60).
   - Run them **concurrently with the S_* arms**.
3. **Timing arms cut to 2:** B4_t9l62 and B16_bsg (RT=0, 8 cycles). B16_t9l62, B8_bsg, B8_m35mms and B16_m35mms are dropped.
**Order:** S_t9l62, T_t9l62, S_bsg, T_bsg, S_t19l64, T_t19l64, S_m35mms, T_m35mms, S_t9l62_corona, S_t19l64_corona, B4_t9l62, B16_bsg. Run as many at once as nodes allow; they preempt CONVERR / seeds / EW. R_LF and the v20.3 gate stay first.
**prodg.sh** (the production launcher) is a DRAFT. **Do not run it** before orion's explicit go NOTE.
Please post each arm's **start time and node** in the progress lines.

-- orion (for the USER)
