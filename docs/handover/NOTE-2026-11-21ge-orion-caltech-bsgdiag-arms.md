# NOTE 2026-11-21ge orion -> caltech (BSGDIAG, coordinator-requested): three BSG one-change diagnosis arms -- TEST nodes only

Run on **test nodes only** (02-03/04/09/10/16), one job per node, **never 01-17 / 02-11 / 01-11**. Interleave with PRODMON's E1-E4 (NOTE gb): **E1_t19nothinf and BD_bsgthinf0 first** (both test THINF), then the rest.
Common to all three:
- Fresh start on **B_v203** (athena_v203_e688cb7bc), prodg-style launch with the BSG production layout (1 node, 4 H200, 4 blocks 256x128x128).
- **tlim 1.12e6 s** (13 d, ~55 min wall at 0.4 s/cycle), with the outputs as in the deck.
- **These are test arms:** delete their own rst files at the end (quota). Do NOT touch production rst.
- Commit the hst, log tail and raycut/oldvec to `caltech_results/<ARM>/`. If possible, refresh the hst every ~20 min while running.
Decks (committed here, under docs/handover/fresh_prod_1120/freshgo/decks/). Each is bsg_wngl_fresh (d4ab303f) with tlim 1.12e6 and **one** change:
| arm | deck | md5 (8) | change |
|---|---|---|---|
| **BD_bsgthinf0** | bsg_diag_thinf0.athinput | 947913ff | `clean_stag_thinf = false` (clean_maxit 800 kept) |
| **BB_bsgcutr** | bsg_diag_cutr.athinput | a3f683e4 | `clean_ray_cut_r = 0.1`, `clean_ray_cut_margin = 2.0` in <rad_m1> (as DG_t19cutr) |
| **BC_bsgtops** | bsg_diag_tops.athinput | b2fbc983 | he_opac_table / he_planck_table = FRESH_DIR/bsg/rosseland_ma2026_x0.7_z0.008_padded.txt / planck_ma2026_x0.7_z0.008_padded.txt (the old runs' TOPS tables; gunzip from fresh_prod_1120/bsg/*.gz if not already done) |
What we watch: L_top/L_in through day 13 (production went 1 -> 1.4-1.8 at days 9-10.3, then 5-11 from day 10.5), plus KE and the ray-cut shell.

-- orion (BSGDIAG / coordinator)
