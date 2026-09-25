---
name: merge-finished-branches
description: user 09-23 - merge each agent branch into rt-integration as soon as its gates pass (default-off / bitwise), keep rt-integration the latest version
metadata:
  type: feedback
---

User (09-23 evening): "remember to merge updates to rt-integration once they finish."
Every agent works on its own branch (worktrees on /viper/ptmp2/jinma/wt_*). When an agent reports
its gates passed (new behaviour behind a switch, default off or bitwise; or a gated bug fix):
merge the branch into rt-integration right away (git merge --no-edit, message listing the switch),
check `git diff <branch> HEAD --stat` is only docs, and remove untracked README copies that collide.
Pushing to the fork (git push fork rt-integration) was approved 09-23; push after merges.

**Why:** the user wants one latest version; stale side branches were confusing ("clean-rt-in?").
**How to apply:** never merge a branch whose gate failed or is still PENDING (e.g. GPU bitwise not yet run);
defaults changes get their own gated branch. Productions keep their own binaries.
Related: [[use-rt-integration-branch]].
