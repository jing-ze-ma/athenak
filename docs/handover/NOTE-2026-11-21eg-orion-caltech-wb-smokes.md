# NOTE 2026-11-21eg orion -> caltech (USER FINAL): WB on for all fresh stars -- 6 WB smoke arms at the HEAD

The USER's final decision is WB on for every fresh star. This **supersedes NOTE ef**:
- Drop TC_t9l62_corona and TC_t19l64_corona if they have not started; the WB decks replace them.
- M150_bsg / M150T_bsg may finish if they are already running.

The decks are in eaebb9634 under docs/handover/fresh_prod_1120/freshgo/decks/ (copy them into fresh_prod_1120). Run all six **concurrently, one node each**, on B_v203, NR=4, with your armg.sh copy, RT=1, 120 cycles, no extra keys. Keep R_LF on 02-16 running.
1. `RT=1 armg.sh W_t9l62 $B_v203 decks/t9l62_wb_smoke.athinput 120`
2. `RT=1 armg.sh W_bsg $B_v203 decks/bsg_wb_smoke.athinput 120`
3. `RT=1 armg.sh W_t19l64 $B_v203 decks/t19l64_wb_smoke.athinput 120`
4. `RT=1 armg.sh W_m35mms $B_v203 decks/m35mms_wb_smoke.athinput 120`
5. `RT=1 armg.sh W_t9l62_corona $B_v203 decks/t9l62_corona_wb_smoke.athinput 120`
6. `RT=1 armg.sh W_t19l64_corona $B_v203 decks/t19l64_corona_wb_smoke.athinput 120`
Please commit res / hst / log tails for each arm as soon as it finishes; one commit per arm is fine.

The production decks are now the `*_wb_fresh` decks, not *_fresh_thinf as NOTE ef said. Production still waits for the explicit go NOTE; do not run prodg.sh before it.

-- orion (for the USER)
