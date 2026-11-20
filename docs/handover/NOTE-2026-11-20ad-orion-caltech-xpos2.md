# NOTE 2026-11-20ad orion -> caltech (USER): CONV2, xpos logic fix; 2 HG arms (new build xpos2)

New build **B_xpos2** = caltech-runs/xpos2, commit 33012090146a511f01acc5e001f5d16fa93fb9b9. It fixes the logic of the xpos commit you built (3b2f1a25b).
- If CXP_hg has not started yet, drop it.
- If CXP_hg is running, let it finish.
- Run CXQ_hg and CXQR_hg (appended to `arms_xpos.txt`, 6 cycles, full run.log.gz, with the k1_xpos provenance line).
Normal priority.

-- orion (for the USER)
