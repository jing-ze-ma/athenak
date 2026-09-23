---
name: no-heavy-monitoring
description: User (09-15): do not monitor runs that much - no rolling watches on production chains; one startup check (PCR line, ranks, first cycles, FATAL) is enough, then check only when asked or at the next planned decision point
metadata:
  type: feedback
---

After the B-star and He chains were verified at start (solver line, ranks, first turnover), I kept arming 30-min watches and per-event notifications; the user said "I don't think you need to monitor that much".
**Why:** each watch costs tokens and attention and the chains are restart-aware with a STOP-on-FATAL; the user checks in themselves.
**How to apply:** verify a launch once (stored script BIN, startup lines, first cycles), record it in memory, and stop. Re-check only when the user asks or when a decision depends on it (link break, gate). Same for test arms: one bounded wait, no rolling heartbeat. See [[save-tokens-everywhere]], [[sbatch-snapshots-script]].
