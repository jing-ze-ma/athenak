# TASK raven -> DeltaAI: swap the running BSG chain to the restart-fixed OLD-method binary (d3rst-1013 @ 9ceb9cdc)

USER decision: DeltaAI's BSG run CONTINUES on the old method (d3bf42d7) as the comparison run; it is NOT switched to FIN3.
It only gets the restart fix, so that its links restart bitwise from now on.

## Why
- In d3bf42d7, implicit_gas_src = conservative2 under hesdirk2 does not store its cell-F slope (csk1, csk1_ok, gcs_first) in the rst.
- So every link restart cold-starts that slope: about 1e-8 relative in L_top from the first step after a restart, on every
  GPU/CPU tested. Continuous runs are unaffected.
- Raven (worker RS) fixed it without any method change.

## Branch: `d3rst-1013` @ `9ceb9cdc` = d3bf42d7 + three changes
- the conservative2 slope in the rst;
- the he_star_m1 Picard window in the rst;
- link-end outputs keep the output schedule. A link that ends on nlim or wall time no longer appends an extra hst row or dump,
  so chained runs reproduce the continuous run's outputs row for row. If any of your scripts read "the last hst row of a link",
  use the scheduled rows instead.

There is NO method change. Continuous runs are bitwise to d3bf42d7.

## Raven verification (reduced BSG, old deck)
| check | result |
|---|---|
| continuous vs d3bf42d7 | BITWISE |
| restart mid-block / at block boundary | BITWISE |
| an OLD d3bf42d7 rst read by the new binary vs d3bf42d7 itself | BITWISE (cold start, exactly as before) |
| chained runs (AG Car A, 3 chain types) | match the continuous run row for row |

One-way only: the new binary reads old rst files, but d3bf42d7 can NOT read rst files written by the new binary.
Once swapped, do not go back to d3bf42d7.

## Steps
1. Build d3rst-1013 @ 9ceb9cdc with exactly the flags of your d3bf42d7 build (athena_hes_gpu_d3bf42d7...). Report md5.
2. Quick check before the swap (one short job, <= 1 h, 1 node, in a SCRATCH dir, NOT your run dir, no run lock):
   - BSG deck (your current deck, unchanged), fresh t = 0, 25 cycles continuous with an rst at cycle 12 (`output5/dt=1050`);
   - restart in a new srun from that rst to cycle 25;
   - compare the cycle-25 dumps and the user.hst rows 13-25: must be BITWISE.
   - Also confirm the continuous 25 cycles match your earlier d3bf42d7 smoke (same deck), i.e. no method change.
3. Swap: between two links, point your link script / config at the new binary and its md5. Deck, run dir, lock, 12 h batch job
   and everything else stay as they are. The next link restarts from your newest rst. That one restart is still a cold start of
   the slope (the old rst has no slope data); every later restart is bitwise.
4. If step 2 is not bitwise, do NOT swap. Keep running d3bf42d7 and report.

## Report: NOTE-2026-10-13-deltaai-bsg-d3rst.md (push to bsg-files-1009)
- build md5;
- the step-2 verdict;
- the link at which the swap happened;
- then per link as usual.
