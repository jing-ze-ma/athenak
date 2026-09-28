---
name: login-node-ram-cap-pool-size
description: viper login node user cgroup is 108 GB; Pool(24) on 6 GB dhj snapshots OOM-killed the whole session 09-25 00:44 (2nd time)
metadata:
  type: feedback
---
The viper login node caps the user slice at 108 GB (memory.max 108017184768). If exceeded, EVERY process is killed,
including Claude itself (core dumps land in the cwd, 5-8 GB each).

- 09-21: session crash (see [[rad-m1-design]]).
- 09-25 00:44: the prod4-vs-hyd4 field-relaxation agent ran /viper/ptmp2/jinma/mhdrelax_0924/extract.py with
  Pool(24), each worker ~6-8.6 GB -> ~160 GB -> everything died; 50 python cores (163 GB) + 4 claude cores (22 GB).

**Why:** agents pick worker counts from core count, not RAM.
**How to apply:** every analysis brief states the budget: total RSS <= 60 GB on the login node
(workers x per-snapshot size), or run it as a CPU batch job on a compute node. Check per-worker size before scaling up.

**Guard installed 09-25 (user approved):** ~/.claude/hooks/memcap.sh (PreToolUse Bash hook in ~/.claude/settings.json,
viper only) runs every Bash command in systemd user slice claude-work.slice (MemoryMax=70G, unit file
~/.config/systemd/user/claude-work.slice), shared by all sessions and agents. A runaway command is killed alone; Claude
survives (tested with an 80 GB allocation). Side effect: `cd` does not persist between Bash calls. Also: after a crash,
check for orphaned multiprocessing parents (Pool respawns workers) and old sessions' watcher shells.
