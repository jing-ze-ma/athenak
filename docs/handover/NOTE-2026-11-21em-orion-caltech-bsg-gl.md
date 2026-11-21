# NOTE 2026-11-21em orion -> caltech (USER): 2 BSG GL-ray arms, right after the WN_* / WE_* arms

The USER asked for these. They must not delay production. The decks are in bede0abcb under docs/handover/fresh_prod_1120/freshgo/decks/ (copy them into fresh_prod_1120). Run on B_v203, 1 node each, NR=4, with your armg.sh copy, RT=1, 120 cycles, no extra keys:
1. `RT=1 armg.sh G_bsg62 $B_v203 decks/bsg62_wbgl_smoke.athinput 120`: GL 7x18 product set, fixed.
2. `RT=1 armg.sh GR_bsg62 $B_v203 decks/bsg62_wbglrot_smoke.athinput 120`: GL 7x18 rotating about z (ring + EMA 128 + ema_wb).
Please commit res / hst / log tails per arm as before.

Thanks for 3ea85be5d: we now read the R_t9l62 / R_t19l64 rc 137 as your 16:15:53 stop, not an OOM.

-- orion (for the USER)
