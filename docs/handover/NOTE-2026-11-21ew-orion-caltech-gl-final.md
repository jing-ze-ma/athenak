# NOTE 2026-11-21ew orion -> caltech (USER, URGENT): GL 7x18 fixed for all fresh stars -- 7 GF_* smokes at the HEAD, concurrently

The USER has decided on **GL 7x18 fixed rays for every fresh star**. The decks are in f54faabcb under docs/handover/fresh_prod_1120/freshgo/decks/ (copy them into fresh_prod_1120; all inputs are files you already have). Run on **B_icbal**, 1 node each, NR=4, with your armg.sh copy, RT=1, **100 cycles**, no extra keys, all **concurrently at the HEAD**.
- These preempt the v20.4 / ST / STD arms; stop those if they hold nodes.
- Running F_* arms that will finish within a few minutes may complete; otherwise these take priority.
- Keep R_LF on 02-16.
1. `RT=1 armg.sh GF_t9l62 $B_icbal decks/t9l62_fing_smoke.athinput 100`
2. `RT=1 armg.sh GF_t19l64 $B_icbal decks/t19l64_fing_smoke.athinput 100`
3. `RT=1 armg.sh GF_m35mms $B_icbal decks/m35mms_fing_smoke.athinput 100`
4. `RT=1 armg.sh GF_t9l62_corona $B_icbal decks/t9l62_corona_fing_smoke.athinput 100`
5. `RT=1 armg.sh GF_t19l64_corona $B_icbal decks/t19l64_corona_fing_smoke.athinput 100`
6. `RT=1 armg.sh GF_bsg $B_icbal decks/bsg_fing_smoke.athinput 100` (BSG, balance anchored at 2.0e12)
7. `RT=1 armg.sh GF_bsgnb $B_icbal decks/bsgnb_fing_smoke.athinput 100` (BSG, no balance)
- If any arm stops at setup with a `vet_gd_gl_nphi` (or other GL) FATAL, **commit that line immediately**; do not wait for the others.
- Otherwise commit res / hst / log tails per arm as each finishes.

Production will launch on the `*_fing_fresh` decks on B_icbal, but only after the explicit go NOTE.

-- orion (for the USER)
