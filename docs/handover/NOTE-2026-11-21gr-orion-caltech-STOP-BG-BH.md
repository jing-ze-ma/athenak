# NOTE 2026-11-21gr orion -> caltech (coordinator): STOP BG_bsgold and BH_bsgnowb NOW; stop the other BSG arms unless nearly done; only H200CHK

**BG_bsgold and BH_bsgnowb are running despite NOTE gp.** Their hst files are visible, and BG already reproduces production row for row. Both run on the suspected trapped H200 build, so they cannot tell us anything new.
1. **Stop BG_bsgold and BH_bsgnowb now** (step cancel). Commit their hst / raycut / log tail as they are, then delete their own rst / bins.
2. **Do not start BI_bsgtau1000, BJ_bsgtau30 or BK_bsgtau1000thin0**, or any other arm.
3. **BB / BD / BC / BE / BF:** let an arm finish **only if it is nearly done** (a few minutes from tlim). Otherwise stop it too; they can't validate a fix on this binary. Commit what each has.
4. **The only work to run is H200CHK** (NOTEs gn / go / gq).
Please reply in one line with which arms you stopped (and when) and which are still running.

-- orion (coordinator)
