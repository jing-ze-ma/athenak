# NOTE 2026-11-09 ORION -> viper: raycheap6 -- gas-driven closure change, pulse2 re-test, vet_source = gas

Answer to `TASK-2026-11-09-orion-raycheap6.md`. No code changes.

**Summary**
- **PASS: none.** Best finest-pair p in (b) with the tangential tensor is 1.68 (PRt = RBt); X9t is 0.95.
  No arm reaches 1.8.
- In (b) the arms fall into TWO families that agree to ~2e-8 at the finest level:
  {L, X9, sg_gasX9} with p 0.80/0.95 (1st order), and {X15, PR, RB} with d ~10x smaller and p -0.2/1.67.
  **RB is bitwise identical to PR**, and **RBt to PRt**, at both levels checked. The tangential tensor only adds a
  ~6.7e-5 shift (X9t vs X9, PRt vs PR) and leaves p unchanged.
- If a cheapest arm in the better family is wanted: PRt (4243 s finest), the same cost as X9t (4262 s). RBt costs
  1.35x more (5727 s) and gives bitwise the same result.
- (a) pulse2 re-test with the new binary: **U1 (rebuild) and U5 (predict) are still bitwise identical to U3 (lag)**
  at coarsest and finest, both t_end; U2 vs U3 is 5.4e-3 / 5.0e-4 (T2.5) and 5.8e-4 / 4.7e-5 (T5), the same as in raycheap5.
  The test inner flux at t_n + dt in the predicted build therefore does not reach the solution.
- (c) **vet_source = gas changes nothing measurable**. The key is read (no fatal; the `relaxed_esc (default)` log
  line disappears), but sg_gasX9 = gas_X9 bitwise at the finest level (1.6e-8 at the coarsest),
  sg_gasX9t = gas_X9t bitwise (finest), sg_gasPRt = gas_PRt bitwise (both levels), sg_pulseP0 d = P0 d to all
  4 digits (1.0e-7 from U3), and **sg_forcedB0 d = 1.805e-4, 8.932e-5, 4.258e-5 is raycheap3 B0's d digit for digit**
  (p 1.01/1.07). So the source-relaxation suspect is not confirmed, at least not through this key in this setup
  (vet_gd_replace = true; whether the vet_gd path reads the vet_source switch at all is a question for the code).

## Setup
- MOD `f872899da`: `git diff --stat a1481b185 f872899da -- src CMakeLists.txt` is empty. Fresh build (gcc/13, openmpi/4.1,
  Release, MPI, `-ffp-contract=off`), binary kept at /orion/ptmp/jinma/raycheap6_1109/athena_f872899da.
- 88 runs, 24 ranks each, in one 12-node exclusive job (227986), 48 concurrent steps, ordered by tlim/cfl
  (all 48 long runs started at once, the 40 short ones filled in behind). Whole job 1 h 36 min; critical path gas_RBt_c0.01875 (5727 s).
- **All 88: rc=0, FATAL=0, NOT ADMISSIBLE 0, DT_RETRY 0; every run 1 backward-Euler step (the first), stage fallbacks 0,
  time2 fallbacks 0, posfix 0.** No arm is spoiled.

## (b) gasmod (t_end 400, cfl 0.15 / 0.075 / 0.0375 / 0.01875)

| arm | d | p | wall (s) | vet_col tensor builds |
|---|---|---|---|---|
| L lag | [1.045e-05, 6.017e-06, 3.114e-06] | [0.80, 0.95] | 546/968/1683/3429 | 297/594/1187/2374 |
| X9 extrap2 (9) | [1.061e-05, 6.017e-06, 3.114e-06] | [0.82, 0.95] | 553/978/1698/3449 | 297/594/1187/2374 |
| X15 extrap2 (15) | [9.487e-07, 1.166e-06, 3.652e-07] | [-0.30, 1.67] | 546/965/1683/3475 | 297/594/1187/2374 |
| PR predict | [1.021e-06, 1.166e-06, 3.654e-07] | [-0.19, 1.67] | 543/962/1678/3417 | 297/594/1187/2374 |
| RB rebuild | [1.021e-06, 1.166e-06, 3.654e-07] | [-0.19, 1.67] | 723/1312/2363/4774 | 593/1187/2373/4747 |
| **X9t** | [1.043e-05, 6.017e-06, 3.112e-06] | [0.79, 0.95] | 681/1185/2039/4262 | 297/594/1187/2374 |
| **PRt** | [9.374e-07, 1.166e-06, 3.644e-07] | [-0.31, 1.68] | 669/1180/2019/4243 | 297/594/1187/2374 |
| **RBt** | [9.374e-07, 1.166e-06, 3.644e-07] | [-0.31, 1.68] | 864/1542/2747/5727 | 593/1187/2373/4747 |

PASS check (finest-pair p >= 1.8 with tangential): X9t 0.95 FAIL, PRt 1.68 FAIL, RBt 1.68 FAIL.
Note the PR family's first pair is non-monotone (d rises from 1.0e-6 to 1.2e-6), so its 1.67 is not a clean asymptotic value.

Final-E comparisons (max|dE|/max|E|, coarsest c0.15 / finest c0.01875):
- X9 vs L 1.0e-3 / 2.1e-8; X15 vs L 9.9e-4 / 1.7e-4; PR vs L 1.0e-3 / 1.7e-4
- X15 vs PR 2.3e-4 / 1.7e-8; **RB vs PR 0 / 0 (identical)**; **RBt vs PRt 0 / 0**
- X9t vs X9 6.7e-5 / 6.7e-5; PRt vs PR 1.0e-3 / 6.7e-5
- predict and rebuild counters are identical as well (finest: `limited cell-builds = 4.136048e6 of 1.049739e9` in both)

## (c) srcgas (vet_source = gas)

| arm | d | p | wall (s) | builds |
|---|---|---|---|---|
| sg_gasX9 | [1.061e-05, 6.017e-06, 3.114e-06] | [0.82, 0.95] | 546/972/1692/3475 | 297/594/1187/2374 |
| sg_gasX9t | [1.043e-05, 6.017e-06, 3.112e-06] | [0.79, 0.95] | 685/1198/2051/4281 | 297/594/1187/2374 |
| sg_gasPRt | [9.374e-07, 1.166e-06, 3.644e-07] | [-0.31, 1.68] | 682/1193/2043/4266 | 297/594/1187/2374 |
| sg_forcedB0 | [1.805e-04, 8.932e-05, 4.258e-05] | [1.01, 1.07] | 567/1029/1859/3701 | 299/598/1196/2392 |
| sg_pulseP0 T2.5 (cfl 0.0375..0.0046875) | [6.042e-03, 1.811e-03, 5.956e-04] | [1.74, 1.60] | 19/30/58/102 | 8/15/30/60 |
| sg_pulseP0 T5 | [1.154e-03, 3.910e-04, 1.635e-04] | [1.56, 1.26] | 32/58/103/203 | 15/30/60/120 |

For reference: gas twins above; raycheap3 B0 d = 1.805e-4, 8.932e-5, 4.258e-5, p 1.01/1.07; raycheap4 P0 T2.5 p 1.74/1.60, T5 1.56/1.26.

## (a) pulse2 re-test (cfl 0.0375 / 0.01875 / 0.009375 / 0.0046875)

| arm | d | p | wall (s) | builds |
|---|---|---|---|---|
| U1 T2.5 | [6.042e-03, 1.811e-03, 5.956e-04] | [1.74, 1.60] | 23/37/73/139 | 15/29/59/119 |
| U1 T5.0 | [1.154e-03, 3.910e-04, 1.635e-04] | [1.56, 1.26] | 37/73/142/274 | 29/59/119/239 |
| U2 T2.5 | [6.371e-03, 2.156e-03, 7.951e-04] | [1.56, 1.44] | 19/29/53/101 | 8/15/30/60 |
| U2 T5.0 | [1.156e-03, 3.910e-04, 1.638e-04] | [1.56, 1.26] | 31/56/101/206 | 15/30/60/120 |
| U3 T2.5 | [6.042e-03, 1.811e-03, 5.956e-04] | [1.74, 1.60] | 18/29/54/103 | 8/15/30/60 |
| U3 T5.0 | [1.154e-03, 3.910e-04, 1.635e-04] | [1.56, 1.26] | 31/56/107/195 | 15/30/60/120 |
| U5 T2.5 | [6.042e-03, 1.811e-03, 5.956e-04] | [1.74, 1.60] | 19/31/54/101 | 8/15/30/60 |
| U5 T5.0 | [1.154e-03, 3.910e-04, 1.635e-04] | [1.56, 1.26] | 33/59/104/205 | 15/30/60/120 |

max|dE|/max|E| vs U3 (c0.0375 / c0.0046875): U1 0 / 0 and U5 0 / 0 (bitwise) at both T; U2 5.377e-3 / 4.984e-4 (T2.5),
5.790e-4 / 4.652e-5 (T5). The predict counter shows the predicted builds ran (U5 T5 c0.0375: limited cell-builds 4.17e6 of 6.19e6).

## Per-run lines (rc, FATAL, NOT ADMISSIBLE, DT_RETRY, wall, builds | hesdirk2 counters)
```
gas_L_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 546 s builds=297 | hesdirk2: stage steps=296 backward-Euler steps=1 stage fallbacks=0
gas_L_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 968 s builds=594 | hesdirk2: stage steps=593 backward-Euler steps=1 stage fallbacks=0
gas_L_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1683 s builds=1187 | hesdirk2: stage steps=1186 backward-Euler steps=1 stage fallbacks=0
gas_L_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 3429 s builds=2374 | hesdirk2: stage steps=2373 backward-Euler steps=1 stage fallbacks=0
gas_X9_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 553 s builds=297 | hesdirk2: stage steps=296 backward-Euler steps=1 stage fallbacks=0
gas_X9_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 978 s builds=594 | hesdirk2: stage steps=593 backward-Euler steps=1 stage fallbacks=0
gas_X9_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1698 s builds=1187 | hesdirk2: stage steps=1186 backward-Euler steps=1 stage fallbacks=0
gas_X9_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 3449 s builds=2374 | hesdirk2: stage steps=2373 backward-Euler steps=1 stage fallbacks=0
gas_X15_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 546 s builds=297 | hesdirk2: stage steps=296 backward-Euler steps=1 stage fallbacks=0
gas_X15_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 965 s builds=594 | hesdirk2: stage steps=593 backward-Euler steps=1 stage fallbacks=0
gas_X15_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1683 s builds=1187 | hesdirk2: stage steps=1186 backward-Euler steps=1 stage fallbacks=0
gas_X15_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 3475 s builds=2374 | hesdirk2: stage steps=2373 backward-Euler steps=1 stage fallbacks=0
gas_PR_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 543 s builds=297 | hesdirk2: stage steps=296 backward-Euler steps=1 stage fallbacks=0
gas_PR_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 962 s builds=594 | hesdirk2: stage steps=593 backward-Euler steps=1 stage fallbacks=0
gas_PR_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1678 s builds=1187 | hesdirk2: stage steps=1186 backward-Euler steps=1 stage fallbacks=0
gas_PR_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 3417 s builds=2374 | hesdirk2: stage steps=2373 backward-Euler steps=1 stage fallbacks=0
gas_RB_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 723 s builds=593 | hesdirk2: stage steps=296 backward-Euler steps=1 stage fallbacks=0
gas_RB_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1312 s builds=1187 | hesdirk2: stage steps=593 backward-Euler steps=1 stage fallbacks=0
gas_RB_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 2363 s builds=2373 | hesdirk2: stage steps=1186 backward-Euler steps=1 stage fallbacks=0
gas_RB_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 4774 s builds=4747 | hesdirk2: stage steps=2373 backward-Euler steps=1 stage fallbacks=0
gas_X9t_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 681 s builds=297 | hesdirk2: stage steps=296 backward-Euler steps=1 stage fallbacks=0
gas_X9t_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1185 s builds=594 | hesdirk2: stage steps=593 backward-Euler steps=1 stage fallbacks=0
gas_X9t_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 2039 s builds=1187 | hesdirk2: stage steps=1186 backward-Euler steps=1 stage fallbacks=0
gas_X9t_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 4262 s builds=2374 | hesdirk2: stage steps=2373 backward-Euler steps=1 stage fallbacks=0
gas_PRt_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 669 s builds=297 | hesdirk2: stage steps=296 backward-Euler steps=1 stage fallbacks=0
gas_PRt_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1180 s builds=594 | hesdirk2: stage steps=593 backward-Euler steps=1 stage fallbacks=0
gas_PRt_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 2019 s builds=1187 | hesdirk2: stage steps=1186 backward-Euler steps=1 stage fallbacks=0
gas_PRt_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 4243 s builds=2374 | hesdirk2: stage steps=2373 backward-Euler steps=1 stage fallbacks=0
gas_RBt_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 864 s builds=593 | hesdirk2: stage steps=296 backward-Euler steps=1 stage fallbacks=0
gas_RBt_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1542 s builds=1187 | hesdirk2: stage steps=593 backward-Euler steps=1 stage fallbacks=0
gas_RBt_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 2747 s builds=2373 | hesdirk2: stage steps=1186 backward-Euler steps=1 stage fallbacks=0
gas_RBt_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 5727 s builds=4747 | hesdirk2: stage steps=2373 backward-Euler steps=1 stage fallbacks=0
sg_gasX9_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 546 s builds=297 | hesdirk2: stage steps=296 backward-Euler steps=1 stage fallbacks=0
sg_gasX9_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 972 s builds=594 | hesdirk2: stage steps=593 backward-Euler steps=1 stage fallbacks=0
sg_gasX9_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1692 s builds=1187 | hesdirk2: stage steps=1186 backward-Euler steps=1 stage fallbacks=0
sg_gasX9_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 3475 s builds=2374 | hesdirk2: stage steps=2373 backward-Euler steps=1 stage fallbacks=0
sg_gasX9t_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 685 s builds=297 | hesdirk2: stage steps=296 backward-Euler steps=1 stage fallbacks=0
sg_gasX9t_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1198 s builds=594 | hesdirk2: stage steps=593 backward-Euler steps=1 stage fallbacks=0
sg_gasX9t_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 2051 s builds=1187 | hesdirk2: stage steps=1186 backward-Euler steps=1 stage fallbacks=0
sg_gasX9t_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 4281 s builds=2374 | hesdirk2: stage steps=2373 backward-Euler steps=1 stage fallbacks=0
sg_gasPRt_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 682 s builds=297 | hesdirk2: stage steps=296 backward-Euler steps=1 stage fallbacks=0
sg_gasPRt_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1193 s builds=594 | hesdirk2: stage steps=593 backward-Euler steps=1 stage fallbacks=0
sg_gasPRt_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 2043 s builds=1187 | hesdirk2: stage steps=1186 backward-Euler steps=1 stage fallbacks=0
sg_gasPRt_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 4266 s builds=2374 | hesdirk2: stage steps=2373 backward-Euler steps=1 stage fallbacks=0
sg_forcedB0_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 567 s builds=299 | hesdirk2: stage steps=298 backward-Euler steps=1 stage fallbacks=0
sg_forcedB0_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1029 s builds=598 | hesdirk2: stage steps=597 backward-Euler steps=1 stage fallbacks=0
sg_forcedB0_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1859 s builds=1196 | hesdirk2: stage steps=1195 backward-Euler steps=1 stage fallbacks=0
sg_forcedB0_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 3701 s builds=2392 | hesdirk2: stage steps=2391 backward-Euler steps=1 stage fallbacks=0
upd_U1_T2.5_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 23 s builds=15 | hesdirk2: stage steps=7 backward-Euler steps=1 stage fallbacks=0
upd_U1_T2.5_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 37 s builds=29 | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0
upd_U1_T2.5_c0.009375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 73 s builds=59 | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0
upd_U1_T2.5_c0.0046875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 139 s builds=119 | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0
upd_U1_T5.0_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 37 s builds=29 | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0
upd_U1_T5.0_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 73 s builds=59 | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0
upd_U1_T5.0_c0.009375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 142 s builds=119 | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0
upd_U1_T5.0_c0.0046875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 274 s builds=239 | hesdirk2: stage steps=119 backward-Euler steps=1 stage fallbacks=0
upd_U2_T2.5_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 19 s builds=8 | hesdirk2: stage steps=7 backward-Euler steps=1 stage fallbacks=0
upd_U2_T2.5_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 29 s builds=15 | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0
upd_U2_T2.5_c0.009375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 53 s builds=30 | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0
upd_U2_T2.5_c0.0046875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 101 s builds=60 | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0
upd_U2_T5.0_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 31 s builds=15 | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0
upd_U2_T5.0_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 56 s builds=30 | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0
upd_U2_T5.0_c0.009375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 101 s builds=60 | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0
upd_U2_T5.0_c0.0046875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 206 s builds=120 | hesdirk2: stage steps=119 backward-Euler steps=1 stage fallbacks=0
upd_U3_T2.5_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 18 s builds=8 | hesdirk2: stage steps=7 backward-Euler steps=1 stage fallbacks=0
upd_U3_T2.5_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 29 s builds=15 | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0
upd_U3_T2.5_c0.009375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 54 s builds=30 | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0
upd_U3_T2.5_c0.0046875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 103 s builds=60 | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0
upd_U3_T5.0_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 31 s builds=15 | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0
upd_U3_T5.0_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 56 s builds=30 | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0
upd_U3_T5.0_c0.009375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 107 s builds=60 | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0
upd_U3_T5.0_c0.0046875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 195 s builds=120 | hesdirk2: stage steps=119 backward-Euler steps=1 stage fallbacks=0
upd_U5_T2.5_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 19 s builds=8 | hesdirk2: stage steps=7 backward-Euler steps=1 stage fallbacks=0
upd_U5_T2.5_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 31 s builds=15 | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0
upd_U5_T2.5_c0.009375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 54 s builds=30 | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0
upd_U5_T2.5_c0.0046875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 101 s builds=60 | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0
upd_U5_T5.0_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 33 s builds=15 | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0
upd_U5_T5.0_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 59 s builds=30 | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0
upd_U5_T5.0_c0.009375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 104 s builds=60 | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0
upd_U5_T5.0_c0.0046875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 205 s builds=120 | hesdirk2: stage steps=119 backward-Euler steps=1 stage fallbacks=0
sg_pulseP0_T2.5_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 19 s builds=8 | hesdirk2: stage steps=7 backward-Euler steps=1 stage fallbacks=0
sg_pulseP0_T2.5_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 30 s builds=15 | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0
sg_pulseP0_T2.5_c0.009375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 58 s builds=30 | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0
sg_pulseP0_T2.5_c0.0046875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 102 s builds=60 | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0
sg_pulseP0_T5.0_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 32 s builds=15 | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0
sg_pulseP0_T5.0_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 58 s builds=30 | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0
sg_pulseP0_T5.0_c0.009375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 103 s builds=60 | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0
sg_pulseP0_T5.0_c0.0046875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 203 s builds=120 | hesdirk2: stage steps=119 backward-Euler steps=1 stage fallbacks=0
```
(time2 and posfix counters are 0 in every run.)

## Files
/orion/ptmp/jinma/raycheap6_1109: rc6.sbatch, rc6.out, ana.sh, rc6_ana.out, cmp.py, runs/*/run.log + wall.txt,
athena_f872899da (kept). bins, sources and the build dir deleted after the comparisons.

