# NOTE 2026-11-21er orion -> caltech (USER): 4 stretched-grid arms right AFTER the F_* / MB_bsg62m smokes (this releases the hold in NOTE ep / eq)

Production does not wait for these, so the F_* / MB smokes and the go-NOTE launch come first. The decks are in 6bcac7065 under docs/handover/fresh_prod_1120/freshgo/decks/ (copy them into fresh_prod_1120). Run on B_icbal, 1 node each, NR=4, with your armg.sh copy, RT=1, 120 cycles, no extra keys:
1. `RT=1 armg.sh ST_t9l62 $B_icbal decks/t9l62_st_smoke.athinput 120` (WB, no balance)
2. `RT=1 armg.sh ST_t19l64 $B_icbal decks/t19l64_st_smoke.athinput 120`
3. `RT=1 armg.sh STD_t9l62 $B_icbal decks/t9l62_stbd_smoke.athinput 120` (WB + discrete-force balance)
4. `RT=1 armg.sh STD_t19l64 $B_icbal decks/t19l64_stbd_smoke.athinput 120`
Please keep the **full a/run.log** for these 4 and commit it gzipped (caltech_results/<TAG>/a_run.log.gz) with the res / hst as usual. We need the dt-limit lines if the code prints them; if the full log is too big, the tail is fine.

-- orion (for the USER)
