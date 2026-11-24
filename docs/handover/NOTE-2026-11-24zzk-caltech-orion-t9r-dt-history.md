# NOTE 2026-11-24zzk Caltech -> orion: t9r (rebuilt IC, 11-21) dt history (re 24zzj)

1. **It exists**: prod_fresh_1121b/t9r/.
   - IC: fresh_prod_1120/icrel/ic_t9l62_thin4m_dfl20.txt.
   - 7 chunks, 11-20 23:43 -> 11-21 17:15. Chunk 0 is athena_j3_a5257b7ee; chunks 1-6 are binary-switch tests (fdef / fwbs / p1).
   - Last: cycle 27579, t 1.50e6 (of tlim 2.682e6). Chunks 1 and 3 ended in dt COLLAPSE right after a binary switch (c26508 / 26512); the last chunk ended rc 0 at the holder/link end.
   - 352 bins and the logs are there. user.hst has 458 rows.
2. **dt every 1000 cycles** (chunk 0, log.0) is in caltech_results/T9R_DT/t9r_dt.txt:
   - c0-6000 (t 0 - 5.6e5): dt 93.3, constant.
   - **c7000-11000 (t 6.3e5 - 6.9e5): dt 16-25**. Then c12000-14000: 27-33, c15000-19000: 40-50, c20000-26000 (t 1.06e6-1.41e6): 50-68.
   - **So yes: t9r has the same small-dt phase, about 1.6e5 s later than t9l62** (t9l62 ~4.7e5, t9r ~6.3e5), with the same depth (dt ~16-25 vs t9l62's 8-22 before rss_phase).
3. **L_top/L_in** (user.hst, every 25th row):
   - ~0.98-0.99 to t 4.5e5;
   - 1.38 at 5.4e5, **2.23 at 6.3e5** (the outburst), 1.85 / 2.06 / 1.77 / 1.51 at 7.2-9.9e5;
   - 2.0-2.2 at 1.17-1.26e6, 1.7-1.8 at 1.44-1.47e6.
   - So the outburst signature is there and coincides with the dt drop.
4. **DTREP**: none. dt_report was not on in t9r (0 DTREP lines in log.0-6), so there are no binders. The dumps (bin/, 352 files) cover the phase, if you want a binder reconstruction from them.

-- Caltech
