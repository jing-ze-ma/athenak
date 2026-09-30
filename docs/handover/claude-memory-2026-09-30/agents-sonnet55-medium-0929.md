---
name: agents-sonnet55-medium-0929
description: ALL agents = Opus 5.5 medium effort via worker (user 09-29 ~06:00); supersedes the Sonnet 5.5 worker setting of earlier 09-29
metadata:
  type: feedback
---
User 09-29 ~06:00: "for all opus 5.5 medium effort". .claude/agents/worker.md now has model: claude-opus-5-5, effort: medium.
Every delegated task uses subagent_type "worker" (pass model "opus" as well if the session predates the edit).
**Why:** the Sonnet worker mislabelled inputs and left causes open on the He presn IC gate; the user replaced it.
**How to apply:** no Sonnet split for "simple" tasks any more; supersedes [[simple-agents-sonnet-0929]] and the Sonnet lines in [[delegate-to-sonnet-when-adequate]].
