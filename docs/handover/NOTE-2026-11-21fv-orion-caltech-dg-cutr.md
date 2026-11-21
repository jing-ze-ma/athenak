# NOTE 2026-11-21fv orion -> caltech (PRODMON): DG_t19cutr on the next free TEST node, ahead of DD / DE

DF_t19nocut (ray cut off) blew up at the wall within ~10 cycles. DG tests whether keeping the sc ray cut near the photosphere stops the 18-30 Rsun collapse.
- Deck: **docs/handover/fresh_prod_1120/freshgo/decks/t19l64_wngl_cutr_diag.athinput** (md5 starts dad46e45). It is t19l64_wngl_fresh plus `clean_ray_cut_r = 0.1` and `clean_ray_cut_margin = 2.0` written into <rad_m1>. Copy it into your fresh_prod_1120 decks dir; the keys are in the deck, not on the CLI.
- Run: `RT=0 KEEPBIN=1 armg.sh DG_t19cutr $B_v203 decks/t19l64_wngl_cutr_diag.athinput 520` with only the common output overrides (`time/nlim=520 output1/dt=600 output2/dt=14400 output3..7/dt=1e30`).
- Test nodes only (02-03/04/09/10/16), one job per node, **ahead of DD_t19conv / DE_t19noseed**.
- Commit the hst, raycut.txt, eint<0 summary, profiles and log tail, with the deck md5 in info.txt.

-- orion (PRODMON)
