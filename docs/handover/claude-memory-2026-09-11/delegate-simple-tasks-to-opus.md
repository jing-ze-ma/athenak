---
name: delegate-simple-tasks-to-opus
description: User wants Fable 5.1 to reserve its own tokens for thinking, bug-finding and design; simple tasks (reading, editing, running jobs, monitoring) go to Opus 5 subagents
metadata:
  type: feedback
---

Fable 5.1 should do the thinking (ideas, finding bugs, inventing solutions) itself, and hand
simple mechanical work to subagents running on Opus 5: reading files and summarising, applying
edits that are already specified, building, launching and chaining Slurm jobs, and monitoring
runs or logs.

**Why:** Token economy. The user said on 2026-09-09: "let's save tokens by letting fable 5.1
think about ideas, find bugs, invent solutions but distribute simple tasks like reading, editing,
running jobs, monitoring to opus 5".

**How to apply:** START EVERY SESSION THIS WAY, from the first tool call, not after the user
reminds you: the user repeated this on 2026-09-09 ("remember you can distribute simple tasks to
opus 5. I told you that last time ... I would prefer that to save token") after Fable had done
builds, merges, smoke tests and log reads itself. Reading a handover note, git merge/build/test
cycles, creating run directories from a spec, and waiting on jobs are ALL subagent work. Use `Agent` with `model: "opus"` (general-purpose or Explore) for file
reads across many files, spec'd edits, job submission and polling. Keep the analysis and the
decision in the main session; give the subagent a precise, self-contained brief and ask it to
return only the conclusion, not file dumps. Do not fork (forks inherit the Fable model). Single
known-file lookups are still fine to do directly. See [[validate-the-instrument]] for what a
subagent's "clean" report must be checked against.
