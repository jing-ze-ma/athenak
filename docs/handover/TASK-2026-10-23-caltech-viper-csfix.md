# TASK 2026-10-23 caltech -> viper (USER): test + merge csfix-1023 into fin51-r-fix, then report back

Branch `csfix-1023` = ONE commit 958cac78 on fin51-r-fix c1681dbe (pushed). USER: "merge that into production and make that default",
"let viper do the test", "let viper do the merge and get back".

## What it fixes (RSG E dt collapse at 230.1 d, 340 s -> 2 s; diagnosis in Caltech rsg_E_1020/collapse_test_1023/RESULTS.txt)
1. REAL BUG (rad_m1_vetgd.cpp, rad_m1_vetgd_cs.cpp): a dt_retry failed attempt on the LAST step of a vet_rotate_avg block left the
   vet_gd direction tables on the NEXT angle set, while the checkpoint restored vgd_alpha to the current one; the retried step skipped
   VetGdTables (al == vgd_alpha) -> crash "compact lists: N remote band reads not in the read mask", or (compact band off) a silently
   wrong closure. Also: a cubed-sphere compact-list entry rebuilt by the failed attempt kept a stale key. Fix: tables of the restored
   angle remade (same numbers), rebuilt entries dropped (generation counters outside the checkpoint).
2. DEFAULT CHANGE: `<time>/dt_retry_t2_be` 2 -> 0 (only acts with time2_reject = retry AND dt_retry on): a non-admissible hesdirk2
   stage is redone with backward Euler at the cycle dt instead of cutting dt (2 cuts trapped dt at the ~2.3 s stage-2 limit).

## Caltech evidence (on ff1f52bc, binary ..._ff1f52bc_csfix2_nofma md5 d9d1192a, + t2_be=0 explicit)
- arm A1 (old binary + t2_be=0): FATAL at cycle 32893 (above message). A4 (cache fix only): same fatal bitwise. A6 (both fixes): passes,
  dt 1000 s, 4.4 d / 27-min segment, ~19 % BE cycles, 2 retries, no NC. A5 (t2_be=1): no fatal but dt 300-1000, 116 retries.
- RSG E production on csfix2 + t2_be=0 since 230.6 d: now 263 d, dt 1000 s steady, no fatal.

## Please
1. Gates as usual (build + bitwise default for a deck WITHOUT time2_reject = retry / dt_retry, restart, cpplint). A deck that uses
   time2_reject = retry + dt_retry and does not name dt_retry_t2_be changes by design.
2. If you can: an RSG-like retry test that hits a rotation-block-end retry (message "<rad_m1> dt_retry (csfix-1023): vet_gd tables back").
   Caltech can supply A6/rsg.00093.rst (230.6 d, 24 MB x 548x64x64) if wanted.
3. Merge into fin51-r-fix, push, and post a NOTE viper->caltech with the merge hash + gate results. Caltech will rebuild the RSG E/HG
   production binaries from it (RSG keeps the csfix2 binary until then; Caltech did NOT run its own gate, user handed it to viper).
