---
name: no-blocking-waits
description: user rejects tool calls that block for minutes waiting on jobs (until-loops, long sleeps); check once quickly and report
metadata:
  type: feedback
---

Do not run blocking waits on Slurm jobs (until ! squeue ...; sleep loops, `sleep 240`): the user interrupted two such calls on 2026-09-19 and asked "how are runs" instead.
**Why:** it freezes the conversation; the user wants to steer while jobs run. Same spirit as [[no-heavy-monitoring]].
**How to apply:** submit, verify startup with ONE quick non-blocking look (or a run_in_background waiter), then return control; read results when the user asks.
