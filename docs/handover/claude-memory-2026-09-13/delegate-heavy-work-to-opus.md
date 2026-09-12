---
name: delegate-heavy-work-to-opus
description: "STANDING RULE, re-affirmed 2026-09-10 (apply from the FIRST tool call of EVERY session, no reminder needed): hand ALL mechanical work -- reading big code regions, editing, building, launching runs, MONITORING/polling jobs, bulk analysis -- to Opus 5 subagents (Agent tool, model opus). The main model keeps only diagnosis, the precise brief, and VERIFYING every diff and number that comes back."
metadata:
  type: feedback
---

**What the user said:** "Can you save some token by getting ideas and inferring possible
bugs but you let the real hard work be done by opus 5?" and "like easy tasks such as reading
and editing, but you should make sure the editing is correct", then "remember to hand over
heavy tasks to opus 5".

**2026-09-10 01:35, mid-session:** "remember you can hand over simple tasks to opus 5" -> "like
monitoring" -> "you should remember this everytime we come back. it should be in your memory".
I had armed Monitor tasks myself; the user wants even job polling done by an Opus agent that
then runs the first analysis when the event fires.  Do not use Monitor/Bash polling here.

**Why:** token economy on the main model; Opus agents are cheaper for mechanical work and
run in parallel.  The user still wants correctness owned here.

**How to apply:**
- Delegate: multi-file reads, code edits, builds, sbatch, POLLING/WATCHING jobs (the agent
  sleeps 60 s between checks and does the first analysis on the event), dump analysis, plots.
- At session start: read the state note, then spawn the agents in ONE message; never poll
  or grep the logs myself beyond the first orientation look.
- Keep: root-cause reasoning, the spec (exact anchors, success criteria, numbers to
  report), reviewing the returned diff (`git diff` of the touched region) and sanity
  checking the reported numbers before acting on them.
- Brief format that worked: measured context first, numbered tasks, explicit verification
  with pass/fail thresholds, "report under N lines with tables", the traps (module load in
  the same command as make; never run athena on the login node; command-line override of
  a missing parameter is fatal; pin the binary per run; 90-column style).
- Continue an existing agent with SendMessage when the follow-up is in the same code --
  it keeps the context and saves the re-brief.
- Coordinate concurrent agents on the same file/build tree explicitly (targeted edits only,
  wait for another `make`), see [[two-stream-rt-three-kernels-trap]] for what one missed
  kernel cost.

**2026-09-11 06:00 addendum (user: "even opus agent is using a lot of tokens. why?")**: a
single agent given audit + 30-file inventory + 4 builds + bisection + note burned 270k tokens
/ 200 tool calls and kept re-firing stale wait timers after finishing. Rules: ONE deliverable
per agent; name the exact hunks/files to read ("do not survey the tree"); no bisection by
rebuild unless asked; tell every agent "when the report is delivered STOP, leave no timers or
background jobs armed"; consider sonnet for pure polling.
