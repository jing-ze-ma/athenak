---
name: index-standing-rules
description: Full original text of the standing user rules and working-practice directives kept short in MEMORY.md.
metadata:
  type: reference
---

## Standing rules — full lines

- **[VERIFY SPECTRAL VERDICTS (user 09-15): look at the maps; argmax of a broad spectrum is not 'box-bounded'; pulsation maps are not convection scales](verify-spectral-verdicts.md)**
- **[NO BLOCKING WAITS (user 09-19): no until/sleep loops on jobs in the foreground; submit, one quick look, return control](no-blocking-waits.md)**
- **[NO HEAVY MONITORING (user 09-15): verify a launch once, then stop; no rolling watches](no-heavy-monitoring.md)**
- **[SBATCH SNAPSHOTS THE SCRIPT: cancel + resubmit after editing submit.sh; verify with scontrol write batch_script (09-15 stale-binary incident)](sbatch-snapshots-script.md)**
- **[SAVE TOKENS](save-tokens.md) — no unneeded work; narrow prompts, cheap models for lookups**
- **[DELEGATE simple tasks to OPUS 5 subagents](delegate-simple-tasks-to-opus.md) — from the first tool call; Fable decides, Opus executes**
- **[STANDING RULE: delegate to Opus 5 agents](delegate-heavy-work-to-opus.md) — EVERY session from the first tool call: Opus agents do reading/editing/building/launching/MONITORING/analysis; I do diagnosis, briefs, VERIFY every diff and number; ONE deliverable per agent, name the hunks to read, "STOP, no timers" (token blow-up 09-11)**
- **[apudev FOR SUB-15-MIN JOBS (user rule 09-15)](apudev-for-short-jobs.md)**
- **[SAVE TOKENS (standing rule 09-11): especially on Fable; no unnecessary tool calls, no surveys, agents one deliverable each, no timers left armed](save-tokens-everywhere.md)**
- **[COMPARE WITH v_MLT AND F_conv/F, NOT WITH box_w4 (user 09-17): w4 was 19x v_MLT; He w9 already 0.6 v_MLT at 3 turnovers](compare-with-mlt-not-old-runs.md)**
- **[viper INODE QUOTA (09-17): worktree add fails while big writes work; delete build dirs in frozen worktrees; test with a 3000-file loop](viper-inode-quota-builddirs.md)**
- **[USE rt-integration FOR EVERYTHING incl. He4 (user 09-21): he4-presn-global merged as e688efd4 and frozen; main checkout = rt-integration](use-rt-integration-branch.md)**
