# NOTE 2026-11-12 Orion -> viper: rw6-paper (IDORT paper linear-wave ladders, amp 1e-6, FINAL vs PROD)

Builds: PROD 7d233ac30 (md5 423209affc92f7e157c24c3194a97357), FINAL 84cd20888 (8ec46383d524802271739678056f13c4); Release, MPI, -ffp-contract=off.
Lists orion_paper_{FINAL,PROD}.txt as posted, run in reverse order (slow cases first). One rank per case, 16 per arm, p.shared orion003.
Wall: PROD 27939 s, FINAL 28148 s; the slowest were xy N128 tau 0.01, ~7.7 h each.

**Result: FINAL 54/54 ok, PROD 54/54 ok, 0 FAILED, 0 FATAL.**

results.json of every case: docs/handover/scripts/rw_amp1e6/results_paper/{FINAL,PROD}/<case dir>/results.json (108 files).
Logs: results_paper/runs_paper_FINAL.out, results_paper/runs_paper_PROD.out.

(A backup job re-running the 4 slow xy cases was cancelled when the main job finished in time. A 16-rank variant of those
cases is not possible: vet_x1_periodic needs a single MeshBlock.)
