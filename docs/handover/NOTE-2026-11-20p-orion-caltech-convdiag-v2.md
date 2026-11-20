# NOTE 2026-11-20p orion -> caltech (USER): CONV diag arms changed (kappa 0.01, shorter), +CEX_hg; ahead of all other CONV work

`arms_convdiag.txt` is replaced. All of its arms now add `rad_m1/clean_conv_kappa=0.01`, to reproduce CV_tight_hg's stage-solve
failure at cycle 624962:
- CD_hg runs 2 cycles.
- CE_hg and CXS_hg run 6 cycles.
- New arm CEX_hg (6 cycles) sets ew_off and xt_sticky together.
If the old CD_hg (kappa 1, 3 cycles) already started, let it finish and run CD2_hg too. Otherwise skip CD2_hg.
Run these right after the v20g3 arms, ahead of the CONSERV ablation remainder and batch 2/3.

-- orion (for the USER)
