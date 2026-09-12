---
name: save-tokens-everywhere
description: "STANDING RULE 2026-09-11: save tokens, ESPECIALLY on the main (Fable) model but also in agents; never spend tokens on anything known to be unnecessary"
metadata:
  type: feedback
---
**What the user said (09-11 06:10):** "remember we need to save tokens especially for fable and
also in general. so we should not waste tokens when we know it's not necessary."

**Why:** the main model is the expensive one; agents are cheaper but one over-broad agent still
burned 270k tokens (see [[delegate-heavy-work-to-opus]]).

**How to apply:**
- Main model: no orientation greps beyond the one needed to write a brief; no re-reading files an
  agent will read; replies short; no memory rewrites for trivia; do not answer stale
  agent-timer notifications with anything but a one-line acknowledgement (or nothing).
- Agents: one deliverable each, name exactly what to read, no exploratory surveys, no
  bisection-by-rebuild unless asked, "STOP after the report, leave no timers armed",
  report length capped. Use sonnet/haiku for pure polling or file copying.
- Before any tool call ask: does the outcome change what happens next? If not, skip it.
