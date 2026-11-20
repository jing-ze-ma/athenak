# NOTE 2026-11-20j orion -> caltech (USER): CONV batch amendment, +1 arm and a new order

One arm is added to `arms_conv.txt` (same B_conv 9e8234b9b):
`FINALRST=0 CV_k01_hg B_conv HG 20 hg_v20.keys time/dt_retry_host=m1 rad_m1/clean_nl_print=true rad_m1/clean_xt_local=true rad_m1/clean_conv_df=1 rad_m1/clean_conv_kappa=0.1`

New CONV order:
1. CV_ref_rsg, CV_tight_rsg, CV_estdf_rsg, CV_est_rsg, CV_k01_rsg
2. **HG trio CV_ref_hg / CV_k01_hg / CV_tight_hg** (together, on comparable nodes). This trio answers a question the USER gets
   in the morning, so please run it before the rest.
3. CV_estdf_hg, CV_df2_rsg at the end of the CONV group.

-- orion (for the USER)
