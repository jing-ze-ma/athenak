---
name: agent-collision-lesson
description: An agent that hands back with "background work still running" can RESUME itself later; before starting a replacement agent on the same branch, stop the old one (TaskStop) or confirm it is not listed as running
metadata:
  type: feedback
---

**What happened, 2026-09-26.** The first ck-next agent handed back an interim report while its jobs were still running. I started a fresh agent on the same branch and worktree. The old agent then resumed when its jobs finished: it committed cn14 and started diagnostic branches, all in the same scratch directory. The fresh agent detected this and stopped itself.

**Why it matters.** Two agents on one worktree, one scratch directory and one GPU budget overwrite each other's work.

**How to apply**
- Before spawning a replacement agent for the same branch, either:
  - TaskStop the old agent (it may appear as completed "with background work still running", yet still resume); or
  - resume the old agent itself with SendMessage.
- Give each agent its own scratch subdirectory when several could overlap.
- Also stop finished agents that stay listed as running: the Kokkos 5 agent lingered for hours.
