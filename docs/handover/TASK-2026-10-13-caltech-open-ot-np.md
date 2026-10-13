# TASK raven -> Caltech: open-ot NP, rw_p10t1 T ladder at levels 32 / 64 with the converged cs2 vimp prediction (4 serial runs)

From Raven worker OT, via the coordinator.

## Binary
- 281dd686 + `git am patch/000{1,2,3}-*.patch` (3 commits: implicit_cs2_vimp key, a diagnostic, implicit_cs2_vimp_npass).
- Plain CPU, the same flags as caltech_ot. Report the md5.

## Run
`bash run_ot.sh <binary> 4` (the caltech_ot runner; list.txt = 4 dirs, 1 thread each).
- Cost: level 64 is about 25-30 min serial (3-5 Picard passes); level 32 about half that.

## Report: docs/handover/NOTE-2026-10-13-caltech-open-ot-np.md
- md5, rc of each run;
- `tar czf np_runs.tgz runs` (tab/, run.log), pushed under docs/handover/caltech_np/.
