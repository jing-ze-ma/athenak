---
name: delegate-to-agents
description: Delegate work to Opus 5.5 (medium effort) agents; AVOID simple models (haiku/sonnet) -- user 2026-09-25; never raise effort to high
metadata:
  node_type: memory
  type: feedback
  originSessionId: 1f4629f7-1667-4eeb-be90-3eb8315950c0
  modified: 2026-09-26T05:20:40.157Z
---

Distribute tasks to subagents (Agent tool, default model = Opus 5.5 at medium effort). **Avoid simple models (haiku, sonnet).** Do not switch on high effort (no high-effort guidance docs in agent prompts).

**History:**
- 2026-09-25: the user first asked to save tokens by sending easy tasks (reading, writing, submitting jobs) to cheaper models.
- A Haiku inventory of speed-up ideas across the memory notes came back shallow and garbled. The user then said "you probably need an opus 5.5 agent to read them not haiku", and then "try to avoid simple models".

**How to apply:**
- Use the default (Opus 5.5) for agents.
- Save tokens by delegating instead: narrow one-deliverable prompts, reports of results only, no file dumps. Do trivial single commands myself rather than spawning any agent.
- Never add effort-raising instructions to agent prompts.
