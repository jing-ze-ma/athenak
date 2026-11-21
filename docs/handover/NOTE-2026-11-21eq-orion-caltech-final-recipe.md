# NOTE 2026-11-21eq orion -> caltech (USER): FINAL recipe smokes, 5 arms at the HEAD together with MB_bsg62m

The USER has decided the FINAL recipe for every fresh star: WB + discrete-force he_ic_balance on **B_icbal** (caltech-runs/icbal 4a7f72334, NOTE ep). The decks are in 12eb791a4 under docs/handover/fresh_prod_1120/freshgo/decks/ (copy them into fresh_prod_1120; all inputs are files you already have).
Run on B_icbal, 1 node each, NR=4, with your armg.sh copy, RT=1, 120 cycles, no extra keys, **concurrently** at the HEAD together with MB_bsg62m. MB_bsg62m is the BSG final deck (bsg_fin_smoke is the same deck plus a header line), so it counts as F_bsg; do not run bsg_fin_smoke separately.
1. `RT=1 armg.sh F_t9l62 $B_icbal decks/t9l62_fin_smoke.athinput 120`
2. `RT=1 armg.sh F_t19l64 $B_icbal decks/t19l64_fin_smoke.athinput 120`
3. `RT=1 armg.sh F_m35mms $B_icbal decks/m35mms_fin_smoke.athinput 120`
4. `RT=1 armg.sh F_t9l62_corona $B_icbal decks/t9l62_corona_fin_smoke.athinput 120`
5. `RT=1 armg.sh F_t19l64_corona $B_icbal decks/t19l64_corona_fin_smoke.athinput 120`
Please commit res / hst / log tails per arm the moment each finishes.

The 4 stretched-grid arms stay HELD even if a node is free.

Production will fire on the `*_fin_fresh` decks on B_icbal, but **only** after the explicit go NOTE.

-- orion (for the USER)
