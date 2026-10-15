# NOTE Caltech -> Raven (10-15): dt_retry fix re-gate (fin51-r-fix b16b66db) -- TEST ONLY

Answers TASK-2026-10-15-caltech-dtretry-regate.md (664339dd). No production job touched.

## Build
b16b66db nofma HOPPER90 (fin51 recipe, builds_fin51, incremental, 197 objects recompiled), BUILD_OK.
`/resnick/home/jingze/ATHENAK/builds/bin/athena_gpu_he_star_m1_b16b66db_nofma`, md5 **2df1a7744c95857e42e4904ccf32bc18**.

## Gate
`rol2u_1015/caltech_rol2u/gate_cr2_fix.sh` = copy of gate_cr2.sh with BQ = b16b66db for B, C, L1/L2, F (P = production binary as
before: hg d505d21f, E1 5b7e1bc5), output in `<case>/gate_fix/j<job>`, extra E1 arm F_all = F + `time/dt_retry_views=all`,
first differing hst row printed. Same seeds (hegiant.00060, agcar3d.00052), debug QoS, both on hpc-sm-01-10.
Logs: `rol2u_1015/logs/gatefix_hg.4586986.out`, `gatefix_e1.4586987.out`.

## He giant (4586986, 4 H200, 30 cycles) -- GATE PASS
| gate | rst | hydro.hst | user.hst | bins |
|---|---|---|---|---|
| GATE-BASE P=B | BITWISE | BITWISE | BITWISE | BITWISE |
| GATE-RESTART C=L2 | BITWISE | BITWISE | BITWISE | BITWISE |
| GATE-RETRY C=F (1 DT_RETRY) | BITWISE | BITWISE | BITWISE | BITWISE |

## E1 (4586987, 2 H200, 20 cycles)
| gate | rst | hydro.hst | user.hst | bins |
|---|---|---|---|---|
| GATE-BASE P=B | DIFF | DIFF | DIFF | 9 differ |
| GATE-RESTART C=L2 | BITWISE | BITWISE | BITWISE | BITWISE |
| GATE-RETRY C=F (1 DT_RETRY) | **BITWISE** | BITWISE | BITWISE | BITWISE |
| extra C=F_all (views all, 1 DT_RETRY) | BITWISE | BITWISE | BITWISE | BITWISE |

- **The 4552183 C != F failure is gone with b16b66db**: C = F and C = F_all bitwise. (Old 1bd0c2ce C vs F: first diff hydro row 9,
  user row 6.) C itself is unchanged between 1bd0c2ce and b16b66db (hydro.hst identical): straight runs bitwise, as Raven found.
- **GATE-BASE P != B is the expected consequence of the fix, not a new problem.** P (5b7e1bc5) and B (b16b66db) both take two
  NATURAL retries on this seed with production keys: cycle 38828 (t 2.60017e7, dt 340.426 -> 172.264, resid 0.0211) and cycle 38841
  (NON-CONVERGED, dt 331.063 -> 167.863). Rows 1-5 are bitwise; the first differing row is exactly the retried step:
  - hydro.hst row 6/20, t = 2.6001874391181756e+07 (dt 172.264, the retried cycle 38828): col 6 (3-mom) 8.7875502689635318e+33 vs 8.7875502689641221e+33 (only column).
  - user.hst row 6/20, same t: cols 4, 5, 20 (e.g. col 5 1.2632912460411525e+38 vs 1.2632912460411670e+38).
  So P (old restore, the bug) and B (fixed restore) part at the first restore, at round-off level, and the rst differs at the end.
  P is reproducible (identical hst to the 4552183 P arm). B also differs from the 1bd0c2ce B for the same reason.
  The retry reason label also changed: P prints `reason=DIVERGED`, B `reason=resid>implicit_resid_fatal` (same resid).
- Since F_all = C and F = C, the minimal View set is now sufficient here; no non-View state left on this test.
- Memory: F_all peak 49.9 GiB/GPU vs F 34.7 GiB.

## Verdict
He giant: full pass. E1: GATE-RESTART and GATE-RETRY pass (the earlier mid-block C != F is cured by b16b66db); GATE-BASE cannot
be bitwise against a pre-fix binary on a seed that retries naturally, since the fix changes exactly the post-restore state
(first diff = the retried step, round-off size). A P=B check on E1 would need a seed window without natural retries.
Switching agcE1r / hgR2w to b16b66db at a link boundary is the coordinator's proposal / user decision; nothing was changed here.
