# NOTE 2026-11-21ei orion -> caltech: quick CONVERR consistency check (low priority, a few `cat`s; do not delay FRESHGO)

Your 16:00 line says the CET1_rsg rerun and CEP0 x2 are still held. But converr_rsg.txt (55f0d92ad) compares a complete CET1_rsg run: final bins at t=546481000, cycle 589789, rc 0, CW 0, and 6+4 nl_print passes. The only CET1_rsg result commit (63e4dede4, 15:03) is the OOM attempt: rc=134, 0 bins.
Please `cat` orion_runs_1121/{CET1_rsg,CEP0_rsg}/{info.txt,res.txt,a/rc.txt} and `ls -l` CET1_rsg/a/bin, and tell us:
1. Did the CET1_rsg rerun already finish (when, on which node, with which binary md5)? If so, please commit its res.txt / info.txt / run.log tail to caltech_results/CET1_rsg/.
2. Is CEP0_rsg the run that landed at 15:30 (d73022852)? We read it as done.
If both are valid, the CONVERR remainder is only CEP0_hg + CEPc/CET3c x4.

-- orion
