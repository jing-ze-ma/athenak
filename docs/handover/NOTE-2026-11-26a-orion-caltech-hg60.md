# NOTE 2026-11-26a orion -> Caltech: HG old / new over a longer window, in the same environment

Thank you: e8d3e64d9 is gated on all three platforms (suite, restart chain, HG). Two last arms until 08:00, on 02-03 + 02-09 (T6 environment, hegiant.00205), priority over everything else you have:
1. Old HG binary 1d1972805 (its own keys, HGO) 60 cycles with T6, and e8d3e64d9 (HGN's command) 60 cycles with T6: s/cyc mean and median, dt sequence, DT_RETRY / fallback / BE, wall of the main loop, wall per 1000 s of model time, L_top at the matched end time. The old collapses to dt 22.5 after cycle ~40: this is the number for the goal.
2. If nodes remain: e8d3e64d9 on ONE node (4 H200), 60 cycles: same numbers, peak memory.
No more tests after that. At 08:00 PST stop everything of ours (by step id) and post one line. Production untouched.

-- orion (coordinator)
