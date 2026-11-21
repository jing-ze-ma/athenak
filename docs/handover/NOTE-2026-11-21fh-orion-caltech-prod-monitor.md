# NOTE 2026-11-21fh orion -> caltech (coordinator / PRODMON): production monitoring files, refreshed every ~30 min. READ-ONLY on the runs.

Please do NOT touch the production runs; this is read-only. For each of the 5 production runs (bsg, t19l64, t9l62, m35mms, t19l64_corona in /resnick/groups/carnegie_poc/jingze/prod_fresh_1121/<star>), commit into **docs/handover/caltech_results/PROD_<star>/**:
1. the full `*.hydro.hst` and `*.user.hst`;
2. `run_tail300.txt`: the last 300 lines of the run log;
3. `grep.txt`: `grep -nE 'DT_RETRY|retry|CLEAN WARN|floor|FATAL|NaN' <run.log> | tail -100`;
4. `ls.txt`: `ls -la` of the run dir plus the sizes of its bin/ (and rst/) entries, so we can see the latest hydro_w bin and any radial-profile outputs.
**Refresh every ~30 min** at the same paths (overwrite), one commit per refresh titled "caltech -> orion: PROD monitor <time>", until we say stop.
If a run trips (FATAL, a CW burst, a dt collapse or NaN), commit at once instead of waiting for the next refresh.

-- orion (coordinator)
