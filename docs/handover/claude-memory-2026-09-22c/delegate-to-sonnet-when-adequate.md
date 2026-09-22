---
name: delegate-to-sonnet-when-adequate
description: User 2026-09-22 - to save tokens, give tasks to Sonnet agents when I judge Sonnet can handle them; Opus for judgement-heavy work
metadata:
  type: feedback
---

User, 2026-09-22 03:15: "to save tokens, maybe from now on you can also distribute tasks to sonnet if
you decide they can handle it."

**Why:** Opus agents cost several times more, and the API was overloaded (two 529 failures in a row);
many tasks are mechanical with a binary gate that catches mistakes.

**How to apply:** default to `model: sonnet` for: mechanical edits guarded by bitwise gates, read-only
inventories/greps, plotting and table extraction from hst/dumps, input/script staging, job submission
prep, cleanup runs. Keep `opus` for: new numerical schemes, diagnoses with several hypotheses, anything
where a wrong deletion would not show in a gate (restart readers, shared helpers), and briefs that need
the agent to push back. State the model choice in one line when launching. Verification of every diff
and number by me stays unchanged ([[delegate-heavy-work-to-opus]] still applies for the heavy classes).
