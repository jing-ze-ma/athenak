# NOTE 2026-11-21ep orion -> caltech (USER): build B_icbal; MB_bsg62m (the USER's BSG candidate) at the HEAD

1. **Build B_icbal** from caltech-runs/icbal = 4a7f723342c9a24f18e0de73e857a7756c4e0389. It is v20.3 + he_ic_balance for the Linhao profile with the discrete -dPr/dr force, plus he_ic_balance_aref=file|dpr. Its default behaviour is bitwise to v20.3 (CPU-checked). Please check the .prov / source digest as usual.
2. **HEAD** (first free node, after B_icbal builds): run on B_icbal, 1 node, NR=4, with your armg.sh copy, RT=1, 120 cycles, no extra keys. The deck is in 8ddd358a9 (copy it into fresh_prod_1120; its inputs are the BSG tables you synced plus bsg/ic_ma2026_f0.txt and linhao_profile_full.txt, which you already have):
   - `RT=1 armg.sh MB_bsg62m $B_icbal decks/bsg62mb_wb_smoke.athinput 120`
   This is BSG with the MESA table, mu 0.62, WB, he_ic_balance anchored at 2.0e12, THINF, maxit 800, conv_max 1e-5, tfloor 1000 K.
Please commit res / hst / log tails as soon as it finishes.
3. FYI: the stretched-grid decks (t9l62_st / t19l64_st / t9l62_stb / t19l64_stb, *_smoke) are on the branch too. **Do not run them yet**; their arms come in a later NOTE after a CPU setup check.

Production still waits for the explicit go NOTE.

-- orion (for the USER)
