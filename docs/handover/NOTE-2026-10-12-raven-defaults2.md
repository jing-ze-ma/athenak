# NOTE defaults2-1012: vet_gd_rotate_slope = correct as the FRESH-run default

Worker DF2, 2026-10-12, MPCDF Raven. Work dir /raven/ptmp/jinma/defaults2_1012. Branch defaults2-1012
(worktree /raven/u/jinma/ATHENAK/wt-df2) on bsgcand-1012 d3bf42d7. Code commit 16e48817 (src/rad_m1/rad_m1.cpp, one key).

## Rule (D2's restart-aware pattern, NOTE-defaults.md)
- FRESH run with the physical-time rotation (vet_gd_rotate_dt named > 0, or the 'auto' default / rotate_dt = auto),
  key not named -> **correct** (Time2RotCorrect; acts only under time_scheme = hesdirk2, inert otherwise).
  correct echoes its sub-keys vet_gd_rotate_slope_frac (1 - 1/sqrt 2) and _reuse (true).
- RESTART whose file lacks the key -> keep (bitwise older runs). Files written since vgdrot2 with physical
  rotation already echo the key (keep), so they keep keep.
- Cycle clock (vet_gd_rotate_every) or rotate 0: the key is not read (unchanged).
- Explicit input always wins. The resolved value is echoed (GetOrAddString).

## Default table (rows changed vs NOTE-defaults.md; every other row unchanged)

| key | old default (d3bf42d7) | new fresh default | restart, file lacks key | applies to | code |
|---|---|---|---|---|---|
| vet_gd_rotate_slope | keep | **correct** | keep | physical-clock rotation (rotate_dt > 0 or auto); acts under hesdirk2 | rad_m1.cpp:1051 |
| vet_gd_rotate_slope_frac | (read with correct) | 1 - 1/sqrt 2 (echoed when correct) | - | slope correct | rad_m1.cpp |
| vet_gd_rotate_slope_reuse | (read with correct) | true (echoed when correct) | - | slope correct | rad_m1.cpp |

## Verification (job 31107661, all.sbatch; binaries bin/athena_{cpu,he_cpu}_16e48817 vs prodfinal_1011/bin/*_d3bf42d7)
(a) rlatt_om5_T l 1, 2: U (rotate_dt named, slope not) vs N (slope correct named); A (neither rotate key) vs AN
    (rotate_dt auto + correct named); reduced AG Car A 60 cyc U vs N (FINC deck).
(b) rlatt_om5_T l 1 cycle clock (rotate_every 10) and rotate 0: new vs d3bf42d7 (data + raw files).
(c) d3bf42d7 rst without the key (cycle-clock run, rst t = 2 / AG Car cycle ~30) restarted with rotate_dt on the
    command line: new vs d3bf42d7 restart.
(d) BSG production deck (names correct) reduced to 4 blocks, 5 cycles: new vs d3bf42d7.
RESULTS: see logs/all.31107661.out (pending at the time of writing).
