# NOTE 2026-11-21ef orion -> caltech (USER): 2 corona THINF smoke arms at the head, right after M150_bsg / M150T_bsg

These are the final production setup for the coronas (THINF + clean_maxit 800). The S_*_corona arms logged 2 CW at cycle 0. The decks are in 47a66fb68 under docs/handover/fresh_prod_1120/freshgo/decks/ (copy them into fresh_prod_1120). Run on B_v203, 1 node, NR=4, with your armg.sh copy, RT=1, 120 cycles, no extra keys:
1. `RT=1 armg.sh TC_t9l62_corona $B_v203 decks/t9l62_corona_smoke_thinf.athinput 120`
2. `RT=1 armg.sh TC_t19l64_corona $B_v203 decks/t19l64_corona_smoke_thinf.athinput 120`
Please commit res / hst / log tails as before.

Also:
- F0_bsg: no pushable subset is needed. Keep the npz and bins on /resnick.
- FYI: the USER reversed the WB choice. The production decks are the non-WB `*_fresh_thinf` decks (BSG: bsg_mu616_fresh_thinf). Production still waits for the explicit go NOTE; do not run prodg.sh before it.

-- orion (for the USER)
