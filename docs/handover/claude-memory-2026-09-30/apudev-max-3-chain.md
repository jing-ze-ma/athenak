---
name: apudev-max-3-chain
description: At most 3 chained jobs on apudev (user rule, restated 09-29); never submit longer apudev chains
metadata:
  type: feedback
---
User 09-29 ("don't chain them. remember the rule is at most 3 chains on apudev"): an apudev afterany chain may hold at most
3 links. An agent had submitted a 14-link chain for the He presn r1 continuation; cancelled.
**Why:** apudev is the shared short-job partition (2 nodes); long chains monopolise it.
**How to apply:** put "at most 3 chained apudev jobs" in every agent brief that may use apudev; longer runs go to apu
(tight --time) or are resubmitted in batches of <= 3 after checking. Related: [[apudev-for-short-jobs]].
