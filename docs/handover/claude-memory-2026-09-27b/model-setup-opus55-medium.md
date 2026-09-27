---
name: model-setup-opus55-medium
description: From 2026-09-22 the session and all agents run Opus 5.5 at medium effort (trial); supersedes the Opus 5 / Sonnet delegation rules
metadata:
  type: feedback
---

User 2026-09-22: switched the session model to Opus 5.5 (claude-opus-5-5), medium effort, saved as default.
Agents: ALWAYS launch with subagent_type "worker" (.claude/agents/worker.md: model claude-opus-5-5, effort medium; user order 09-22). If that type is missing, it needs a session restart.
(do not pass sonnet/haiku/fable; "opus" may resolve to an older Opus).

**Why:** trial of Opus 5.5 as the default work set, replacing Fable 5.1; 2.5x cheaper per token than Fable.
**How to apply:** supersedes [[delegate-simple-tasks-to-opus]], [[delegate-heavy-work-to-opus]], [[delegate-to-sonnet-when-adequate]]
for model choice (still delegate, still verify every diff and number). Fable 5.1 only if the user switches back
with /model. Watch for retractions and wasted runs; that is the trial's measure.
