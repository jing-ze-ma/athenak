# NOTE 2026-10-20 caltech -> raven: RSG H200 restart DIFF arms (TASK e3b5c424) — results

Job 5059481 (1 node x 4 H200, a2713be0 e31ea39e, N=48 mb 24x24, nlim 200; log + scripts + cmdlines in docs/handover/csedge-1019/caltech/rdiff/).
- R2 vs R (both restarted from C rst.00001, cycle 100): hydro.hst + user.hst 100/100 rows BITWISE, all bins (hydro_w, m1, m1_vet) data SAME.
  R2 vs C: DIFF from cycle 141 (user.hst) / 145 (hydro.hst), exactly as R vs C before.
  => the restarted run is deterministic on H200; R == R2 != C => STATE MISSING FROM THE RST (or rebuilt differently on restart), not a race.
- CV (continuous, rad_m1/vet_rotate_avg=1) vs C: differs from cycle 0/1 (expected: different key).
- RV vs CV: INVALID pair as specified — RV restarted from C's rst (written WITHOUT vet_rotate_avg=1), so it starts from C's state, not CV's
  (DIFF from cycle 101). Follow-up submitted: RV2 = restart from CV's own rst.00001 (cycle 100) with vet_rotate_avg=1, compared to CV
  (job 5117623, rv2.sh). Result follows in this NOTE.
Note the onset: cycle 141 = restart + 41; in det1 R-C was also cycle 141/145. With vet_gd_rotate_every 10 x vet_rotate_avg 4 = 40 the rotate-avg
buffer remains the prime suspect; RV2 == CV would confirm it.
Command lines of C and R (det1, job 5008037) are in NOTE-2026-10-19-caltech-csedge-rdiff.md.
