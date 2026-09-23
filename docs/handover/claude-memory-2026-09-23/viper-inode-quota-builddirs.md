---
name: viper-inode-quota-builddirs
description: viper GPFS INODE quota hit 09-17 (worktree creation failed while 200 MB writes worked); fixed by deleting 76 build dirs in frozen worktrees; test headroom with a 3000-empty-file loop
metadata:
  type: project
---

09-17 ~10:00: "Disk quota exceeded" on `git worktree add` while `dd` of 200 MB succeeded =
the INODE quota (mmlsquota is not installed; `quota` prints nothing). Diagnose with
`mkdir .itest; for i in $(seq 1 3000); do : > .itest/f$i || break; done` (failed at 708).
Fix: rm build*/ in the frozen worktrees (wt_m3acc had 26, wt_rgbox 27; ~734 files each) ->
20 000 files of headroom. Kept builds only in wt_merge, wt_he4, wt_he4_adi.
Bytes were ALSO full the same morning: prod_w7 rst 00044 was truncated mid-write; user approved
deleting prod_w7 rst 00000-40 (keep newest 3), prod_w5/prod_w6 rst dirs, he4_rg_gate (~50 GB).

**How to apply:** before any new worktree/build or big test, run the inode loop; each build
dir costs ~800 inodes, each handover memory copy ~350. See [[rt-integration-branch]] (GPFS
inode quota note), [[save-tokens-everywhere]].
- ptmp2 RULE (09-23): /viper/ptmp2/jinma is scratch (not backed up, may be purged): only regenerable build/run trees there; patches, READMEs, results on u2 / in the repo.
