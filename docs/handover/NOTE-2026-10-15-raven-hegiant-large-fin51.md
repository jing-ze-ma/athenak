# NOTE Raven -> Caltech (10-15): He giant LARGE -> redo the gate on FIN5.1 (USER); drop the FIN4 C2/C3 BGs

USER decisions:
1. The large He giant runs on **FIN5.1** (fin51-1014 @ 21c4839d; build per TASK-2026-10-15-caltech-fin51-build.md, nofma HOPPER90),
   not FIN4. Reason: the FIN4 fresh-start radiation IC hits the extent_int overflow on the 1697-cell block (Viper FIN4 smoke 12206133
   printed `ic_rad: no fused twin intensities`; twin array 19.8 GB/rank), so gate C on d505d21f started from a wrong IC and the
   C2/C3 separating BGs on FIN4 would be confounded too. **Cancel/skip C2 and C3** (if already run, report them anyway, labelled FIN4).
2. **dt_retry**: production decks set `<time>/dt_retry = true` explicitly (FIN5.1 default on).
Gate on FIN5.1, same deck hegiant_large.athinput (09a3b446 after @FILES@; add `time/dt_retry = true` explicitly; record new md5):
- A smoke + B restart BITWISE (as before); confirm run.log has NO "no fused twin" line and peak GB (FIN5.1 should be well below 77.8 GB).
- C background (seed 0, ~2 h) with the same band analysis vs E0 -> PASS/FAIL as before.
- If PASS: start production (hgLw links, 7 d, tlim 1e9, rmax guard RLIM 525, UCX_TLS=^xpmem). If FAIL: report, hold.
Also: grep the 200-Rsun FIN4 He giant fresh start (4489484 run.log) for the same line (asked in a1f3a175) and report.
