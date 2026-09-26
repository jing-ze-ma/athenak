---
name: save-tokens
description: User rule 2026-09-11 - save tokens, Fable above all but also every subagent; never spend tokens on work that is not clearly needed
metadata:
  type: feedback
---

Save tokens. Fable (the main session) is the most expensive; Opus subagents are not cheap either
(one 7-call "read the handover" agent cost ~80k tokens). Do not spend tokens when the work is
known to be unnecessary.

**Why:** the user raised it on 2026-09-11 after a broad reconnaissance agent was launched for a
next step the user had not chosen, and had to be killed.

**How to apply:**
- Do not launch agents for work the user has not chosen; ask/report first when direction is open.
- Narrow prompts: name files and lines, cap the answer length, tell agents to grep/head instead of
  reading whole files, forbid sub-sub-agents.
- Model by task: haiku/sonnet for lookups; opus for build/test/job cycles; Fable only thinks.
- A single small file read is cheaper than any agent (fixed ~15-20k startup per agent).
- One agent per question; never "be exhaustive" unless the user asked for it.
- Related: [[delegate-simple-tasks-to-opus]] (delegate, but only what is needed).
