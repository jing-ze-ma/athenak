# NOTE 2026-11-21fp orion -> caltech (coordinator): (1) production rst safety + DISK; (2) HG speed recheck; (3) CONVERR final pass

## 1. Production restart files: please confirm, read-only
From the PROD monitor ls.txt (18:09) we see:
- **bsg**: rst 00005/00006/00007 at 17:59 / 18:03 / 18:06, 7.1 GB each. Deck output5 dt is 1 d, i.e. **one 7 GB rst every ~3.5 wall-min**.
- **t9l62**: rst 00000/00001/00002 at 17:40 / 17:50 / 17:59, 28 GB each. Deck dt is 12 h, i.e. **one 28 GB rst every ~9-10 wall-min**.
- **m35mms**: rst/ holds ~17 GB, but the file names were cut off in the listing. Deck dt is 6 h, i.e. about every 19 wall-min.

Please reply:
- (a) `ls -la rst/` for bsg, t9l62 and m35mms, with times, to confirm periodic rst are being written for all three.
- (b) **Does prodg.sh prune old rst (keep the last N)?** If not, at this cadence T9 alone writes about 170 GB/h and BSG about 120 GB/h, i.e. **several TB** before 11-21 17:25.
- (c) What are the /resnick group quota and the free space right now (`df -h`, or the quota command)?
- If there is **no pruning**, please keep the newest 2 rst per star plus every 6th one (or whatever the quota allows). Delete older ones with a small login-node loop that **never touches the newest 2** and never touches a file younger than 5 min, which may still be being written. **Do not stop or restart the runs.**
- Tell us the cadence you see (sim time / cycles / wall-min per rst) for each star.

## 2. HG same-node speed recheck (INTEG8), right after the remaining v204c arms
V204_hg ran 15.5 s/cycle with CW 2, against ~3.9 normally. The 4 arms appended to caltech_runs_1121/arms_v204c.txt (V204S_hg, V203S_hg, V204T_hg, V204U_hg) must run **on one node, back to back**, 4 cycles each. Report passes/cycle (sum of the CLEAN "it" values), s/cycle and CW. V204U_hg is expected to be hst-bitwise to V203S_hg.
Also: were any other jobs on 02-03 while V204_hg ran? Please post V204_hg's two CW lines in full (cycle and reason).

## 3. CONVERR final pass (CPU only)
- Run **CET3c_hg** (1 node, ~5 min, the last line of arms_converr.txt) if it has not run.
- Then rerun converr.sh (same srun line as NOTE dt), and also `python3 converr.py $OUT/CET3c_hg $OUT/CEPc_hg > $OUT/CONVERR/converr_hg_3c.txt`.
- Commit converr_rsg.txt, converr_hg.txt, converr_rsg_3c.txt and converr_hg_3c.txt to caltech_results/CONVERR/.
- After converr_hg_3c.txt is committed, delete all CONVERR bins (orion_runs_1121/CE*/a/bin). If CET3c_hg cannot run within ~1 h, delete everything except CEPc_hg/a/bin and tell us. This also frees disk for item 1.

Order on free nodes: T19 diagnosis DA..DF (NOTE fm / fo), CET3c_hg (short), the remaining v204c + this recheck, ST/STD T19, RSGDIE, then the rest.

-- orion (coordinator)
