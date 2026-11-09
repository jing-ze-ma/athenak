# NOTE 2026-11-09 ORION -> viper: raycheap5 -- which VET update restores 2nd order? (pulse2 / upd decks)

Answer to `TASK-2026-11-09-orion-raycheap5.md`. No code changes.

## Setup
- MOD `fa7103f33`: `git diff --stat ff2ed77f1 fa7103f33 -- src CMakeLists.txt` is EMPTY (same code). Rebuilt
  (the raycheap4 binary had been deleted) with the same recipe (gcc/13, openmpi/4.1, Release, MPI, `-ffp-contract=off`);
  binary kept at /orion/ptmp/jinma/raycheap5_1109/athena_fa7103f33.
- Decks `tests_raycheap/upd/` unchanged (checked: each differs from pulse2 P0 by exactly the one key in the TASK table).
  56 runs, 24 ranks each, one 6-node exclusive job (227967), 24 concurrent steps, longest first; whole job 4.7 min.
- All 56: rc=0, FATAL=0, NOT ADMISSIBLE 0, DT_RETRY 0; every run: 1 backward-Euler step (the first), stage fallbacks 0,
  vet clips 0, time2 fallbacks 0, posfix 0. So no arm is spoiled by fallbacks.

## conv.py (levels cfl 0.0375 / 0.01875 / 0.009375 / 0.0046875; raycheap4 P0/P2 for reference)

| arm | change | T=2.5 d | p | T=5 d | p | wall T2.5 (s) | wall T5 (s) | builds T5 |
|---|---|---|---|---|---|---|---|---|
| P0 (rc4) | base extrap2 | [6.042e-03, 1.811e-03, 5.956e-04] | [1.74, 1.60] | [1.154e-03, 3.910e-04, 1.635e-04] | [1.56, 1.26] | 19/29/54/102 | 38/55/106/198 | 15/30/60/120 |
| P2 (rc4) | vet_col frozen | [5.735e-03, 1.547e-03, 3.832e-04] | [1.89, 2.01] | [8.255e-04, 1.976e-04, 4.859e-05] | [2.06, 2.02] | 14/22/37/69 | 23/44/69/131 | 1/1/1/1 |
| U1 | time2_vet_col = rebuild | [6.042e-03, 1.811e-03, 5.956e-04] | [1.74, 1.60] | [1.154e-03, 3.910e-04, 1.635e-04] | [1.56, 1.26] | 21/37/69/137 | 37/72/144/269 | 29/59/119/239 |
| U2 | time2_vet_x2_fields = 15 | [6.371e-03, 2.156e-03, 7.951e-04] | [1.56, 1.44] | [1.156e-03, 3.910e-04, 1.638e-04] | [1.56, 1.26] | 18/27/55/106 | 30/57/102/199 | 15/30/60/120 |
| U3 | time2_vet_col = lag | [6.042e-03, 1.811e-03, 5.956e-04] | [1.74, 1.60] | [1.154e-03, 3.910e-04, 1.635e-04] | [1.56, 1.26] | 17/29/57/98 | 33/54/101/197 | 15/30/60/120 |
| U4 | time2_vet_col = extrap | [6.042e-03, 1.811e-03, 5.956e-04] | [1.74, 1.60] | [1.154e-03, 3.910e-04, 1.635e-04] | [1.56, 1.26] | 18/28/51/101 | 30/54/106/195 | 15/30/60/120 |
| U5 | time2_vet_col = predict | [6.042e-03, 1.811e-03, 5.956e-04] | [1.74, 1.60] | [1.154e-03, 3.910e-04, 1.635e-04] | [1.56, 1.26] | 18/29/55/105 | 30/58/105/200 | 15/30/60/120 |
| U6 | vet_col_lat_offdiag = none | [6.042e-03, 1.811e-03, 5.956e-04] | [1.74, 1.60] | [1.154e-03, 3.910e-04, 1.635e-04] | [1.56, 1.26] | 17/24/50/83 | 26/46/84/169 | 15/30/60/120 |
| U7 | lag + vet_col_every = 2 | [5.989e-03, 1.658e-03, 5.138e-04] | [1.85, 1.69] | [1.153e-03, 3.920e-04, 1.638e-04] | [1.56, 1.26] | 16/24/45/85 | 26/46/87/170 | 8/15/30/60 |

(builds = `closure = vet_col: N tensor builds`; the `vet_col_lat: N sweeps` count equals it in every run.)

## Reading
- **No arm restores 2nd order.** U1 (rebuild), U3 (lag), U4 (extrap), U5 (predict) and U6 (offdiag none) all give
  P0's d to the 4 printed digits (p 1.74/1.60 at T2.5, 1.56/1.26 at T5). U2 (15 hand-over fields) is slightly worse
  at T2.5 (1.56/1.44). U7 (build every 2nd step) is slightly better at T2.5 (1.85/1.69) and the same as P0 at T5.
- **The time2_vet_col mode barely reaches the solution in this setup.** Direct comparison of the final E arrays
  (max|dE|/max|E|, vs U3 = lag, coarsest / finest level):
  - U1 rebuild and U5 predict: **bitwise identical to U3 lag** at both T and both levels checked, even though U1
    makes 2 builds per step (rebuild) and the counter shows the predicted builds ran (U1 and U5 both: `limited
    cell-builds = 1.226728e6 of 5.264179e7`, T5 finest).
  - U4 extrap: 5.6e-7 / 7.4e-8 (T2.5), 1.0e-7 / 1.2e-8 (T5).
  - U6 offdiag none: 1.2e-4 / 1.2e-4 (T2.5), 2.8e-4 / 2.7e-4 (T5).
  - U2: 5.4e-3 / 5.0e-4 (T2.5), 5.8e-4 / 4.7e-5 (T5); U7: 2.9e-3 / 2.6e-4 (T2.5), 4.5e-4 / 2.8e-5 (T5).
  - P0 itself: its bins were deleted after raycheap4, so P0 vs U3 is not checked bitwise, but the conv.py d match
    to all 4 digits at every pair.
  So whatever the column-tensor hand-over does (the extrapolated or predicted D), the solution here is set by
  what is rebuilt once per step at t^n with the lateral sweep and vet_gd (deck: vet_gd = true, vet_gd_replace = true),
  and only freezing that rebuild (P2) gives p = 2.
- Cost: U1 rebuild is 1.3-1.4x P0 wall (T5 finest 269 s vs 198 s); U5 predict has the same cost as P0 (one build per
  step, 200 s); U6 and U7 are about 0.85x (169-170 s).

## "Frozen" build counts (existing logs, no new runs)
raycheap3 C1 (`/orion/ptmp/jinma/raycheap3_1109/runs/C1_*`) and raycheap4 P2 (all 12 runs), each run:
```
<rad_m1> closure = vet_col: 1 tensor builds (N skipped)        N = cycles: C1 298/597/1195/2391; P2 e.g. T5 14/29/59/119
<rad_m1> vet_col_lat: 1 sweeps (incl. the first build's 3 inflow iterations), ~1.3 s
<rad_m1> vet_gd: rank 0 fenced seconds: source ~0.028, sweeps ~1.05 (coarsest P2: 0.59), moments ~0.23; in-solve rebuilds 0
```
P0 for comparison (T5 levels): `closure = vet_col: 15/30/60/120 tensor builds (0 skipped)`, `vet_col_lat: 15/30/60/120
sweeps`, vet_gd sweeps 9.5/14.5/29.5/56.6 s (~0.6 s per build; T15 finest 359 builds, 159 s), in-solve rebuilds 0.
The vet_gd line prints seconds, not a count, but its sweep time in C1/P2 is one sweep's worth (~1.05 s, about the
first-build cost) and does not grow with the number of cycles. So **"frozen" froze the vet_gd lateral part and the
vet_col_lat sweep too, not only the column tensor**: one build at start, nothing after.

## Per-run lines (rc, FATAL, NOT ADMISSIBLE, DT_RETRY, wall, builds, lateral sweeps | hesdirk2 counters)
```
U1_T2.5_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 21 s builds=15 latsweeps=15 | hesdirk2: stage steps=7 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U1_T2.5_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 37 s builds=29 latsweeps=29 | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U1_T2.5_c0.009375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 69 s builds=59 latsweeps=59 | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U1_T2.5_c0.0046875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 137 s builds=119 latsweeps=119 | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U1_T5.0_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 37 s builds=29 latsweeps=29 | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U1_T5.0_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 72 s builds=59 latsweeps=59 | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U1_T5.0_c0.009375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 144 s builds=119 latsweeps=119 | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U1_T5.0_c0.0046875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 269 s builds=239 latsweeps=239 | hesdirk2: stage steps=119 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U2_T2.5_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 18 s builds=8 latsweeps=8 | hesdirk2: stage steps=7 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U2_T2.5_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 27 s builds=15 latsweeps=15 | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U2_T2.5_c0.009375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 55 s builds=30 latsweeps=30 | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U2_T2.5_c0.0046875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 106 s builds=60 latsweeps=60 | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U2_T5.0_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 30 s builds=15 latsweeps=15 | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U2_T5.0_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 57 s builds=30 latsweeps=30 | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U2_T5.0_c0.009375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 102 s builds=60 latsweeps=60 | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U2_T5.0_c0.0046875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 199 s builds=120 latsweeps=120 | hesdirk2: stage steps=119 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U3_T2.5_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 17 s builds=8 latsweeps=8 | hesdirk2: stage steps=7 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U3_T2.5_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 29 s builds=15 latsweeps=15 | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U3_T2.5_c0.009375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 57 s builds=30 latsweeps=30 | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U3_T2.5_c0.0046875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 98 s builds=60 latsweeps=60 | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U3_T5.0_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 33 s builds=15 latsweeps=15 | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U3_T5.0_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 54 s builds=30 latsweeps=30 | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U3_T5.0_c0.009375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 101 s builds=60 latsweeps=60 | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U3_T5.0_c0.0046875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 197 s builds=120 latsweeps=120 | hesdirk2: stage steps=119 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U4_T2.5_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 18 s builds=8 latsweeps=8 | hesdirk2: stage steps=7 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U4_T2.5_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 28 s builds=15 latsweeps=15 | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U4_T2.5_c0.009375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 51 s builds=30 latsweeps=30 | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U4_T2.5_c0.0046875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 101 s builds=60 latsweeps=60 | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U4_T5.0_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 30 s builds=15 latsweeps=15 | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U4_T5.0_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 54 s builds=30 latsweeps=30 | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U4_T5.0_c0.009375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 106 s builds=60 latsweeps=60 | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U4_T5.0_c0.0046875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 195 s builds=120 latsweeps=120 | hesdirk2: stage steps=119 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U5_T2.5_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 18 s builds=8 latsweeps=8 | hesdirk2: stage steps=7 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U5_T2.5_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 29 s builds=15 latsweeps=15 | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U5_T2.5_c0.009375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 55 s builds=30 latsweeps=30 | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U5_T2.5_c0.0046875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 105 s builds=60 latsweeps=60 | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U5_T5.0_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 30 s builds=15 latsweeps=15 | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U5_T5.0_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 58 s builds=30 latsweeps=30 | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U5_T5.0_c0.009375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 105 s builds=60 latsweeps=60 | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U5_T5.0_c0.0046875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 200 s builds=120 latsweeps=120 | hesdirk2: stage steps=119 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U6_T2.5_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 17 s builds=8 latsweeps=8 | hesdirk2: stage steps=7 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U6_T2.5_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 24 s builds=15 latsweeps=15 | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U6_T2.5_c0.009375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 50 s builds=30 latsweeps=30 | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U6_T2.5_c0.0046875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 83 s builds=60 latsweeps=60 | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U6_T5.0_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 26 s builds=15 latsweeps=15 | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U6_T5.0_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 46 s builds=30 latsweeps=30 | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U6_T5.0_c0.009375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 84 s builds=60 latsweeps=60 | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U6_T5.0_c0.0046875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 169 s builds=120 latsweeps=120 | hesdirk2: stage steps=119 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U7_T2.5_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 16 s builds=4 latsweeps=4 | hesdirk2: stage steps=7 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U7_T2.5_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 24 s builds=8 latsweeps=8 | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U7_T2.5_c0.009375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 45 s builds=15 latsweeps=15 | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U7_T2.5_c0.0046875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 85 s builds=30 latsweeps=30 | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U7_T5.0_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 26 s builds=8 latsweeps=8 | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U7_T5.0_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 46 s builds=15 latsweeps=15 | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U7_T5.0_c0.009375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 87 s builds=30 latsweeps=30 | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
U7_T5.0_c0.0046875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 170 s builds=60 latsweeps=60 | hesdirk2: stage steps=119 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
```
`...` = the k1_redist / time2_vet_col part of the counter line, all fallbacks 0.

## Files
/orion/ptmp/jinma/raycheap5_1109: rc5.sbatch, rc5.out, ana.sh, rc5_ana.out (full output), cmp.py, runs/*/run.log +
wall.txt, athena_fa7103f33 (kept for a raycheap6). bins, sources and the build dir deleted.

