# NOTE 2026-11-09 ORION -> viper: raycheap9 -- the common O(dt) term (Ar4, evolving gas)

Answer to `TASK-2026-11-09-orion-raycheap9.md`. No code changes.

**Summary**
- **Part 0 (HG cs small, cycle 20, predict run):** max|F2|/max|F1| = **3.12e-4**, max|F3|/max|F1| = **2.96e-4**;
  max|F2_PR - F2_RB|/max|F1| = **2.96e-4** (F3: 2.77e-4, F1: 3.43e-4). Confirmed: the lateral flux is ~3e-4 of F1, so
  "0.95 of its own max" is a small absolute change, of the same size as the F1 change.
- **No arm reaches p ~ 2** (t_end 10, finest pair): C4 0.84, C5 0.72, C6 (BE) 1.07, C2 0.50, C3 0.51; C0 0.87.
- **C1 (`dbg_hydro_off = true` on the default hesdirk2 path) is inert**: driver.cpp:739 has
  `nohyd = pm1->dbg_hydro_off && m1mr`, so the key only acts when the multi-rate / coupling_split path is on.
  C1 reproduces C0's d to all 4 digits (8.178e-4, 4.472e-4, 2.443e-4). The hesdirk2 R-alone test cannot be run with this key.
- Added **C1s = coupling_split=strang + dbg_hydro_off=true** (R alone through the multi-rate SDIRK2): p 0.32 / 0.49,
  not 2nd order either. Caveat for C1s, C2, C3: the multi-rate R stage 2 is NOT ADMISSIBLE (min(E,T) < 0) in the FIRST
  3 / 4 / 5 cycles at cfl 0.15 / 0.075 / 0.0375 and is redone with backward Euler; never at 0.01875. That is a BE stretch
  over the same early time window at the three coarser levels, which by itself caps p near 1 and makes the level pairs
  inconsistent, so C1s / C2 / C3 do not decide the question cleanly.
- Reading: the O(dt) term survives without the radiative force on the gas (C4), without the clump (C5), and with a
  different coupling placement (C2 / C3). The clean R-alone check is blocked (C1 inert on hesdirk2; C1s spoiled by the
  early BE fallbacks). A clean next step: C1s with those fallbacks avoided (e.g. a smaller first dt or a later start),
  or a hydro switch on the hesdirk2 path.

## Binary
raycheap2-1109 head is still `4f227b777`; `git diff --stat f872899da 4f227b777 -- src CMakeLists.txt` is EMPTY, so the
raycheap6 binary /orion/ptmp/jinma/raycheap6_1109/athena_f872899da was used. Part 0: the hglagt he_star_m1 binary
/orion/ptmp/jinma/hglagt_1109/wt_new/build/src/athena (7fb7af05a).

## Part 0
The raycheap8 HG end bins had been deleted and the logs do not carry F2/F3 maxima, so the two 20-cycle runs were RERUN
with the raycheap8 recipe (on.athinput copies with time2_vet_col = predict / rebuild, 4 ranks, --mem=48G, nlim 20, m1 bins
every cycle). Both rc=0 (1079 s / 1114 s), t = 2857.14 at cycle 20, 0 NOT ADMISSIBLE / DT_RETRY.
Cycle 20, predict: max E 4.513e9, max F1 4.882e13, max F2 1.525e10, max F3 1.444e10.

| quantity | value |
|---|---|
| max F2 / max F1 | 3.123e-4 |
| max F3 / max F1 | 2.957e-4 |
| max abs(F2_PR - F2_RB) / max F1 | 2.958e-4 (= 0.947 of max F2) |
| max abs(F3_PR - F3_RB) / max F1 | 2.766e-4 (= 0.935 of max F3) |
| max abs(F1_PR - F1_RB) / max F1 | 3.426e-4 |
| max abs(E_PR - E_RB) / max E | 2.335e-9 |

(The raycheap8 Part 1 numbers are reproduced exactly: F1 3.43e-4, F2 0.947, F3 0.935, E 2.3e-9.)

## Part 1 (Ar4 = raycheap8 full_decks evoAr4 with time2_vet_col = predict; t_end 10; cfl 0.15 / 0.075 / 0.0375 / 0.01875)
Keys absent from the deck were written into deck copies (sed) since fresh starts refuse absent command-line keys;
C5 = the same deck without the clump_* lines; C6 = `time_scheme = be`. All runs rc=0, FATAL 0, DT_RETRY 0.

| arm | d | p | wall (s) | NOT ADMISSIBLE | hesdirk2 stage/BE/fallback | builds |
|---|---|---|---|---|---|---|
| C0 (= raycheap8 PR_T10; bins deleted, quoted) | [8.178e-04, 4.472e-04, 2.443e-04] | [0.87, 0.87] | 227/329/484/812 | 0/0/0/0 | 38/1/0 77/1/0 155/1/0 311/1/0 | 39/78/156/312 |
| C1 dbg_hydro_off (INERT, see below) | [8.178e-04, 4.472e-04, 2.443e-04] | [0.87, 0.87] | 230/339/498/804 | 0/0/0/0 | 38/1/0 77/1/0 155/1/0 311/1/0 | 39/78/156/312 |
| C2 coupling_split=strang | [9.120e-04, 5.922e-04, 4.178e-04] | [0.62, 0.50] | 247/365/526/865 | 3/4/5/0 | 0/0/0 0/0/0 0/0/0 0/0/0 | 42/82/160/310 |
| C3 coupling_split=alternate | [8.836e-04, 5.922e-04, 4.164e-04] | [0.58, 0.51] | 255/379/524/904 | 3/4/5/0 | 0/0/0 0/0/0 0/0/0 0/0/0 | 42/82/160/310 |
| C4 dbg_gas_force=false | [8.061e-04, 4.367e-04, 2.447e-04] | [0.88, 0.84] | 174/245/376/615 | 0/0/0/0 | 38/1/0 77/1/0 155/1/0 311/1/0 | 39/78/156/312 |
| C5 no clump | [8.357e-04, 5.060e-04, 3.078e-04] | [0.72, 0.72] | 224/335/490/716 | 0/0/0/0 | 38/1/0 77/1/0 155/1/0 311/1/0 | 39/78/156/312 |
| C6 time_scheme=be | [2.774e-03, 1.355e-03, 6.459e-04] | [1.03, 1.07] | 173/274/419/638 | 0/0/0/0 | 0/39/0 0/77/0 0/155/0 0/309/0 | 39/77/155/309 |
| C1s strang + dbg_hydro_off (R alone; added) | [8.987e-04, 7.199e-04, 5.120e-04] | [0.32, 0.49] | 229/332/470/710 | 3/4/5/0 | 0/0/0 0/0/0 0/0/0 0/0/0 | 37/72/141/271 |

- C2 / C3 / C1s run the multi-rate path, so their hesdirk2 counter line reads 0/0/0; their NOT ADMISSIBLE lines (printed as
  `hesdirk2 stage 2 ... redone with backward Euler`) are in cycles 0-2 (c0.15), 0-3 (c0.075), 0-4 (c0.0375), none at c0.01875.
- C6: `hesdirk2(be_only)` counter line, BE steps 39/77/155/309.
- Finest-level max|dE|/max|E|: C1s vs C2 0.13 (the hydro matters in this test); C4 vs C1 (= C0) 6.3e-3; C5 vs C1 0.31 (different problem).

## Per-run lines (rc, FATAL, NOT ADMISSIBLE, DT_RETRY, wall, builds | hesdirk2 counters)
```
C1_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 230 s builds=39 | hesdirk2: stage steps=38 backward-Euler steps=1 stage fallbacks=0
C1_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 339 s builds=78 | hesdirk2: stage steps=77 backward-Euler steps=1 stage fallbacks=0
C1_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 498 s builds=156 | hesdirk2: stage steps=155 backward-Euler steps=1 stage fallbacks=0
C1_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 804 s builds=312 | hesdirk2: stage steps=311 backward-Euler steps=1 stage fallbacks=0
C2_c0.15 rc=0 FATAL=0 NOTADM=3 DTRETRY=0 wall 247 s builds=42 | hesdirk2: stage steps=0 backward-Euler steps=0 stage fallbacks=0
C2_c0.075 rc=0 FATAL=0 NOTADM=4 DTRETRY=0 wall 365 s builds=82 | hesdirk2: stage steps=0 backward-Euler steps=0 stage fallbacks=0
C2_c0.0375 rc=0 FATAL=0 NOTADM=5 DTRETRY=0 wall 526 s builds=160 | hesdirk2: stage steps=0 backward-Euler steps=0 stage fallbacks=0
C2_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 865 s builds=310 | hesdirk2: stage steps=0 backward-Euler steps=0 stage fallbacks=0
C3_c0.15 rc=0 FATAL=0 NOTADM=3 DTRETRY=0 wall 255 s builds=42 | hesdirk2: stage steps=0 backward-Euler steps=0 stage fallbacks=0
C3_c0.075 rc=0 FATAL=0 NOTADM=4 DTRETRY=0 wall 379 s builds=82 | hesdirk2: stage steps=0 backward-Euler steps=0 stage fallbacks=0
C3_c0.0375 rc=0 FATAL=0 NOTADM=5 DTRETRY=0 wall 524 s builds=160 | hesdirk2: stage steps=0 backward-Euler steps=0 stage fallbacks=0
C3_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 904 s builds=310 | hesdirk2: stage steps=0 backward-Euler steps=0 stage fallbacks=0
C4_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 174 s builds=39 | hesdirk2: stage steps=38 backward-Euler steps=1 stage fallbacks=0
C4_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 245 s builds=78 | hesdirk2: stage steps=77 backward-Euler steps=1 stage fallbacks=0
C4_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 376 s builds=156 | hesdirk2: stage steps=155 backward-Euler steps=1 stage fallbacks=0
C4_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 615 s builds=312 | hesdirk2: stage steps=311 backward-Euler steps=1 stage fallbacks=0
C5_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 224 s builds=39 | hesdirk2: stage steps=38 backward-Euler steps=1 stage fallbacks=0
C5_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 335 s builds=78 | hesdirk2: stage steps=77 backward-Euler steps=1 stage fallbacks=0
C5_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 490 s builds=156 | hesdirk2: stage steps=155 backward-Euler steps=1 stage fallbacks=0
C5_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 716 s builds=312 | hesdirk2: stage steps=311 backward-Euler steps=1 stage fallbacks=0
C6_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 173 s builds=39 | hesdirk2(be_only): stage steps=0 backward-Euler steps=39 stage fallbacks=0
C6_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 274 s builds=77 | hesdirk2(be_only): stage steps=0 backward-Euler steps=77 stage fallbacks=0
C6_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 419 s builds=155 | hesdirk2(be_only): stage steps=0 backward-Euler steps=155 stage fallbacks=0
C6_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 638 s builds=309 | hesdirk2(be_only): stage steps=0 backward-Euler steps=309 stage fallbacks=0
C1s_c0.15 rc=0 FATAL=0 NOTADM=3 DTRETRY=0 wall 229 s builds=37 | hesdirk2: stage steps=0 backward-Euler steps=0 stage fallbacks=0
C1s_c0.075 rc=0 FATAL=0 NOTADM=4 DTRETRY=0 wall 332 s builds=72 | hesdirk2: stage steps=0 backward-Euler steps=0 stage fallbacks=0
C1s_c0.0375 rc=0 FATAL=0 NOTADM=5 DTRETRY=0 wall 470 s builds=141 | hesdirk2: stage steps=0 backward-Euler steps=0 stage fallbacks=0
C1s_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 710 s builds=271 | hesdirk2: stage steps=0 backward-Euler steps=0 stage fallbacks=0
```

## Jobs and files
Job 228050 (5 exclusive nodes; Part 0 + C1-C6 all at once, 19 min) and job 228059 (1 node, C1s, ~12 min).
/orion/ptmp/jinma/raycheap9_1109: decks/ (C1-C6, C1s per level), p0deck/, rc9*.sbatch, ana.sh, rc9_ana.out, runs/*/run.log +
wall.txt, p0/*/log. bins and rst deleted.
