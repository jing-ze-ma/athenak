# NOTE 2026-11-21ex orion -> caltech (USER, URGENT): production recipe is now WN-GL -- 5 WG_* smokes at the HEAD

The USER's production recipe is now **WN-GL**: WB on, **no** rebalance, GL 7x18 fixed, everything else as before. The decks are in 8f922dfce under docs/handover/fresh_prod_1120/freshgo/decks/ (copy them into fresh_prod_1120; inputs are files you already have). Run on **B_v203**, 1 node each, NR=4, with your armg.sh copy, RT=1, 100 cycles, no extra keys, concurrently at the HEAD:
1. `RT=1 armg.sh WG_t9l62 $B_v203 decks/t9l62_wngl_smoke.athinput 100`
2. `RT=1 armg.sh WG_t19l64 $B_v203 decks/t19l64_wngl_smoke.athinput 100`
3. `RT=1 armg.sh WG_m35mms $B_v203 decks/m35mms_wngl_smoke.athinput 100`
4. `RT=1 armg.sh WG_t9l62_corona $B_v203 decks/t9l62_corona_wngl_smoke.athinput 100`
5. `RT=1 armg.sh WG_t19l64_corona $B_v203 decks/t19l64_corona_wngl_smoke.athinput 100`

GF arms (NOTE ew):
- **Drop** GF_bsg, GF_t9l62_corona and GF_t19l64_corona if they have not started. If they have started, stop them.
- **Keep** GF_bsgnb. It is exactly the BSG WN-GL deck (bsg_wngl, no balance), and B_icbal with balance off is bitwise to B_v203, so it serves as the BSG WN-GL smoke.
- **Let GF_t9l62, GF_t19l64 and GF_m35mms run to the end.** They are the GL + rebalance comparison.

Nodes (R_LF stays on 02-16): the 5 WG_* arms take precedence over everything except the three running GF arms and GF_bsgnb. If there are not enough free nodes, start the WG_* arms as nodes free (order as listed).

Commit res / hst / log tails per arm as each lands.

Production = the `*_wngl_fresh` decks on **B_v203**, only after the explicit go NOTE.

-- orion (for the USER)
