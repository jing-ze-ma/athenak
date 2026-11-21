# NOTE 2026-11-21ek orion -> caltech (USER, URGENT): 2 arms at the VERY HEAD (ahead of all R_* rotation arms); T9/T19 production on HOLD

W_t9l62 and W_t19l64 (WB + he_ic_balance) blew up: KE reached 2e44 / 8e42 by cycle 120, and L_bot/L_in went to -112 / +36. These two arms test WB **without** the balance. The decks are in f974e9c18 under docs/handover/fresh_prod_1120/freshgo/decks/ (copy them into fresh_prod_1120). Run on B_v203, 1 node each, NR=4, with your armg.sh copy, RT=1, 120 cycles, no extra keys. Start both on the next two free nodes before any R_* arm that has not started; R_* arms that are already running may continue.
1. `RT=1 armg.sh WN_t9l62 $B_v203 decks/t9l62_wbnb_smoke.athinput 120`
2. `RT=1 armg.sh WN_t19l64 $B_v203 decks/t19l64_wbnb_smoke.athinput 120`
Please commit res / hst / log tails per arm as soon as each finishes.

**HOLD:** do NOT launch T9 or T19 production from the *_wb_fresh decks. Production as a whole still waits for the explicit go NOTE, which will name the T9/T19 decks.

(Letter ek is used because your CONVERR answer already carries "ej".)

-- orion (for the USER)
