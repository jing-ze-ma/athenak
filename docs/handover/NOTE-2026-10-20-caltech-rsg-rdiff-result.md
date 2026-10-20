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

## RV2 result (job 5117623, 10-20 10:00 PDT): restart from CV's OWN rst.00001 (cycle 100), both with rad_m1/vet_rotate_avg=1
- RV2 vs CV: user.hst DIFF from the FIRST row after the restart (t 3.3186e4 = cycle 101), hydro.hst from t 3.45e4 (~cycle 105);
  all bins (hydro_w, m1, m1_vet) DIFF.
=> The rotate-avg buffer is NOT the cause: without rotate averaging the restart diverges even sooner (cycle 101 instead of 141).
   With vet_rotate_avg=4 the averaging apparently masks/delays it. Combined with R == R2 (deterministic restart), some state
   used from the first post-restart step is missing from (or rebuilt differently from) the rst on this a2713be0 build — likely
   VET-related (m1_vet bins differ; user.hst moves before hydro.hst). Log: docs/handover/csedge-1019/caltech/rdiff/rv2.5117623.out.
Not re-run on newer binaries; tell us if you want the same R/R2/RV2 triplet on 563ab474.

## 563ab474 rerun (reply to c7f32c84): job 5140304 submitted 10-20 10:15 PDT (1 node x 4 H200, script rsg_E_1020/rdiff563/rdiff563.sh).
A: deck E N=32 (548x16x16 x 24 blocks), C / R / R2 from C's cycle-100 rst. B: old rsg_det N=48 + vet_rotate_avg=1, CB / RB / RB2.
Results follow. RSG E 5126654 not held (est. start 10-21 01:29 PDT).
