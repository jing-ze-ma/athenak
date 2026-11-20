# NOTE 2026-11-20h orion -> caltech (USER): CONV batch, 9 arms of 20 cycles (pass budget / convergence stop)

(Date fix: our NOTEs 21/21b/21c/21d/21e were really written on 2026-11-20.)

New build **B_conv** = caltech-runs/conv, commit 9e8234b9b04eaa0b936681fd909f8128d0b68bb0. It is v20f plus clean_xt_local plus a
study key clean_conv_df, which defaults to 0 and is then bitwise. Same cmake line as before.
Arms: `docs/handover/caltech_runs_1121/arms_conv.txt` (RSG rsg.01295, HG hegiant.00139, 20 cycles, FINALRST=0, 4 ranks).

Priority, after the v20g gate and batch-1 P1-P3: CV_ref_rsg, CV_tight_rsg, CV_estdf_rsg, CV_est_rsg first, then CV_df2_rsg and
CV_k01_rsg, then the 3 HG arms. Start the RSG arms of this group together on comparable nodes, because they are compared with
each other. This group can go ahead of batch 2/3.

Send back: full run.log gzipped (the 'CLEAN pass' / 'CLEAN DF' lines), `*.user.hst` + `*.hydro.hst`, res.txt.

-- orion (for the USER)
