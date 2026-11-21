# NOTE 2026-11-21el orion -> caltech (USER): WE_t9l62/WE_t19l64 at the head right after WN_*; STOP R_t9l62 and R_t19l64

1. **Stop R_t9l62 (02-09) and R_t19l64 (02-11) now**, if they are still running; these are our own test arms. They are built on the *_wb base, which blew up: all four W_t9l62 / W_t19l64 / W_*_corona arms reached KE ~1e43-2e44 by cycle 120, with L_bot/L_in at -116 / +36, so they cannot judge rotation. We will rebuild them on the base that wins. **R_m35mms and R_bsg62 continue.**
2. Use the freed nodes for these two arms, right after WN_t9l62 / WN_t19l64 (NOTE ek). The decks are in d2736fa5b under docs/handover/fresh_prod_1120/freshgo/decks/ (copy them into fresh_prod_1120). Run on B_v203, 1 node each, NR=4, with your armg.sh copy, RT=1, 120 cycles, no extra keys:
   - `RT=1 armg.sh WE_t9l62 $B_v203 decks/t9l62_wbeq_smoke.athinput 120`
   - `RT=1 armg.sh WE_t19l64 $B_v203 decks/t19l64_wbeq_smoke.athinput 120`
   (WB + balance + he_ic_equil with the HG production keys; the re-march uses the clean steady-start force.)
Please commit res / hst / log tails per arm, and for stopped arms just a one-line note of when you stopped them.

The *_wb decks for T9, T19 and both coronas are on HOLD. Production still waits for the explicit go NOTE.

-- orion (for the USER)
