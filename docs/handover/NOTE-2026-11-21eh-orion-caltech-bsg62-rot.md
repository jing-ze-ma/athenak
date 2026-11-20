# NOTE 2026-11-21eh orion -> caltech (USER): W62_bsg (BSG production candidate) at the head + 4 rotation smokes after it

The decks are in 85ced0780 under docs/handover/fresh_prod_1120/freshgo/decks/ (copy them into fresh_prod_1120). Run on B_v203, 1 node, NR=4, with your armg.sh copy, RT=1, 120 cycles, no extra keys.

**(a) HEAD:** run this on the first node that frees up; it goes before anything else queued. It is the BSG production candidate (mu 0.62); W_bsg (mu 0.616) is now only a side arm.
- `RT=1 armg.sh W62_bsg $B_v203 decks/bsg62_wb_smoke.athinput 120`

**(b) After W62_bsg and the running W_* arms, on free nodes, in parallel.** These are the rotation smokes: HG/RSG production rotating rays + EMA 128 / step / ema_wb. The USER wants them in parallel; they must not delay the fixed-ray production launch.
- `RT=1 armg.sh R_t9l62 $B_v203 decks/t9l62_wbrot_smoke.athinput 120`
- `RT=1 armg.sh R_bsg62 $B_v203 decks/bsg62_wbrot_smoke.athinput 120`
- `RT=1 armg.sh R_t19l64 $B_v203 decks/t19l64_wbrot_smoke.athinput 120`
- `RT=1 armg.sh R_m35mms $B_v203 decks/m35mms_wbrot_smoke.athinput 120`

Please commit res / hst / log tails per arm as before. Production still waits for the explicit go NOTE. The BSG production deck will be bsg62_wb_fresh, or a balance build if the USER agrees; the go NOTE will say which.

Thanks for the profc.sh relative-KEYSFILE catch (absolute path is right).

-- orion (for the USER)
