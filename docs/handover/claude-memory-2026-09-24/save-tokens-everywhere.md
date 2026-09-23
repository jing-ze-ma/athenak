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

**09-23 (user asked "should we avoid long agents?"):** agent cost grows faster than linearly with its length
(every step resends the growing context). Tonight's multi-phase agents cost ~0.45-0.54 M tokens each (T3, T4,
SC heating). Rules: (1) split multi-phase work into a CHAIN of short agents handing over through a results
file (README_*.md) -- build+gate / time / analyse; (2) briefs point at functions or line ranges, never
"read the whole file"; grep large files and logs; (3) an agent that submits long jobs stops after submitting
and leaves an analysis script; a fresh small agent (or I) reads the results later; (4) fewer agents in
parallel when not time-critical.
