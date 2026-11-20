# NOTE 2026-11-20di orion -> caltech (USER): RSGDIE replay arms (R_LF ahead of the seeds)

There are no restart dumps near the deaths, so a deterministic replay with per-pass prints is the way to get the failing step. Two arms in `arms_rsgdie.txt`, on **B_v201e** (the same binary as LF_rsg):
1. **R_LF**: exactly LF_rsg (rsg.01295, NR=4, floor_pay=rad_column) plus rad_m1/clean_nl_print=true, run until it dies (~cycle 592698, ~2.3-2.6 h).
   - Please **write restart dumps at cycles 590230, 592500 and 592690 and keep them** on Caltech. Use an extra rst output with a cycle trigger, or split the run at those cycles, whichever is easiest.
   - First check: hst bitwise to LF_rsg/a over the common rows (nl_print should be diagnostic only). If it is not bitwise, say so and **keep running anyway**.
2. **R_off**: L1_rsg_off (key off) plus clean_nl_print, cycles 589788 -> 590195 (407 cycles; covers off's 6-retry cascade at 590188). FINALRST=1.
**Priority:** R_LF goes ahead of the LFa/L0a/LFb/L0b seeds (NOTE dh) if nodes are short; the seeds are still wanted. R_off is short and can fit anywhere.
Please return the full run.log.gz for both (it will be large, so split it if needed) and the **dump paths**.

-- orion (for the USER)
