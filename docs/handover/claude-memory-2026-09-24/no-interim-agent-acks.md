---
name: no-interim-agent-acks
description: User 09-23 "you don't need a monitor" - do not reply to interim agent notifications (agent stopped/waiting, no report yet); speak only when a result arrives or the user is needed
metadata:
  type: feedback
---

User 2026-09-23, after a string of one-line "X still hasn't reported" replies to interim task
notifications: "you don't need a monitor".
**Why:** those replies add noise and tokens without information.
**How to apply:** when a task notification says the agent has not reported yet / is waiting on its own
background work, produce no user-facing update (at most a minimal internal no-op). Report only
(a) a delivered result, (b) something that needs the user's decision, (c) a failure. See
[[no-heavy-monitoring]], [[save-tokens-everywhere]].
