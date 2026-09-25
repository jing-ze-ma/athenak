---
name: housekeeping-minimal-tokens
description: "user 09-23 - spend minimal tokens on housekeeping (copies, moves, disk cleanup); no progress checks, act only on completion notices"
metadata:
  node_type: memory
  type: feedback
  originSessionId: bf0d9db5-2135-4897-9191-80542c9c35dd
  modified: 2026-09-23T06:26:17.991Z
---

User (09-23): "you shouldn't spend too much token on copying" -- and "the agent also shouldn't" (added to .claude/agents/worker.md). Housekeeping (u2 -> ptmp2 copies,
worktree cleanup, quota) gets a single launch + a single action on completion; no polling, no
throughput diagnosis unless asked.

**Why:** tokens go to the science/code work; copies finish on their own.
**How to apply:** launch background copy with a completion notification, then do the swap
(delete original, symlink) only when notified. User pre-approved the swap of verified dirs
(u2_bench_archive, VERIFIED.txt). See [[save-tokens-everywhere]], [[no-heavy-monitoring]].
