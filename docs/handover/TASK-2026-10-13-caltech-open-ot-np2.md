# TASK raven -> Caltech: open-ot NP2, rlatt_om5_T levels 8 / 16 with the npass / rebuild / prediction-off arms (8 runs)

From Raven worker OT. Use the same binary as NP (281dd686 + patch 0001-0003, md5 0d486e6a).
`bash run_ot.sh <NP binary> 8` (list.txt: 4 runs at 4 threads (level 8), 4 at 8 threads (level 16)).
- Level 16 took about 12 min at 4 threads on Raven.
- Arms:
  - N3 / N5: implicit_cs2_vimp_npass 3 / 5;
  - VR: time2_vet_col rebuild;
  - CS: implicit_cs2_dbg 2.

Report in docs/handover/NOTE-2026-10-13-caltech-open-ot-np2.md: rc, then `tar czf np2_runs.tgz runs` (tab/ bin/ run.log), pushed under
docs/handover/caltech_np2/.
