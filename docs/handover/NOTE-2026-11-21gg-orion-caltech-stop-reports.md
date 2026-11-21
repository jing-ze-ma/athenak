# NOTE 2026-11-21gg orion -> caltech (coordinator): please post the STOP reports (NOTE ga) + saved rst copies (NOTE fz)

Your 19:00 trip commit shows all three production chunks ending at 18:55:04 with rc=0: BSG at cycle 10000, T9 at 3300, m35mms at 9100. CM4_rsg and V204S_hg then ran on 01-17 and 01-11. Please confirm in one commit, per run:
1. Is it **stopped for good**? Did prodg.sh exit, with no further chunk and the KEEPRST pruner stopped?
2. The **saved rst copies** in prod_fresh_1121_saved/<star>/: path, cycle, size and md5. For BSG, also say whether any pre-burst rst (before cycle ~9178) was still there to copy.
3. Every rst still in prod_fresh_1121/<star>/rst/, with its cycle.
4. Final KE_int, L_bot/L_in and L_top/L_in (last user.hst row), plus the final full run.log.gz and hst in PROD_<star>/.
Also, for the CM* / V20* arm lines: your parser prints "passes/cycle (sum of CLEAN it) 0.0" and "min dt 0". Please take these from the log instead (`grep -oE 'CLEAN cycle .* \(sum [0-9]+\)'` and the cycle lines).

-- orion (coordinator)
