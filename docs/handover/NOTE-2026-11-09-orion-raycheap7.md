# NOTE 2026-11-09 ORION -> viper: raycheap7 -- gas_feedback = true (gasfb decks)

Answer to `TASK-2026-11-09-orion-raycheap7.md`. No code changes.

**Summary**
- **PASS: none.** PRt finest-pair p = 1.68, X9t = 0.95 (need >= 1.8).
- **gas_feedback = true changes nothing measurable against raycheap6 (gas_feedback = false):** all six arms give
  raycheap6's d and p to all 4 digits (gas_X9 / gas_PR / gas_X9t / gas_PRt). The log confirms the key is on
  (`coupling=true gas_feedback=true`). The raycheap6 bins are deleted, so this is a d comparison, not a bitwise one.
- **RB vs PR: no difference.** m1 (E, F) and m1_lat bitwise identical at c0.01875, c0.0375 and c0.075; at c0.15,
  E and m1_lat are bitwise and F1/F3 differ by 3.9e-9 / 1.9e-9 relative (max abs 3.7e-9, cmpbin "bitwise 2/4").
- **PRg vs PR: no difference**, exactly the same pattern (bitwise at the three finer levels, F 3.9e-9 at c0.15).
  **RB and PRg are bitwise identical to each other at c0.15**, so the 4e-9 is in PR's coarsest run, not in what RB/PRg change.
  So with gas_feedback on, neither the stage-2 rebuild nor vet_source = gas reaches the operator in this test.
- No arm stopped: all 24 rc=0, FATAL 0, NOT ADMISSIBLE 0, DT_RETRY 0; every run 1 backward-Euler step (the first),
  stage fallbacks 0, accepted-with-fallback 0.

## Setup
- MOD `8532af080`: `git diff --stat f872899da 8532af080 -- src CMakeLists.txt` is EMPTY, so the raycheap6 binary
  /orion/ptmp/jinma/raycheap6_1109/athena_f872899da was reused (no rebuild). Decks and cmpbin.py from `git archive 8532af080`
  (checked: gfb_* = gas_* + `gas_feedback = true`; PRg = PR + `vet_source = gas`).
- 24 runs, 24 ranks each, one 6-node exclusive job (228007), all 24 at once, longest first. Whole job 1 h 34 min;
  critical path gfb_RB_c0.01875 (5658 s).

## conv.py (t_end 400; cfl 0.15 / 0.075 / 0.0375 / 0.01875)

| arm | d | p | wall (s) | vet_col tensor builds | raycheap6 twin d / p |
|---|---|---|---|---|---|
| X9 | [1.061e-05, 6.017e-06, 3.114e-06] | [0.82, 0.95] | 665/1169/2233/4297 | 297/594/1187/2374 | 1.061e-05, 6.017e-06, 3.114e-06 / 0.82, 0.95 |
| PR | [1.021e-06, 1.166e-06, 3.654e-07] | [-0.19, 1.67] | 654/1154/2238/4276 | 297/594/1187/2374 | 1.021e-06, 1.166e-06, 3.654e-07 / -0.19, 1.67 |
| RB | [1.021e-06, 1.166e-06, 3.654e-07] | [-0.19, 1.67] | 832/1503/2905/5658 | 593/1187/2373/4747 | = PR |
| PRg | [1.021e-06, 1.166e-06, 3.654e-07] | [-0.19, 1.67] | 655/1157/2213/4274 | 297/594/1187/2374 | no direct twin (nearest: sg_gasPRt, tangential) |
| **X9t** | [1.043e-05, 6.017e-06, 3.112e-06] | [0.79, 0.95] | 788/1376/2651/5143 | 297/594/1187/2374 | 1.043e-05, 6.017e-06, 3.112e-06 / 0.79, 0.95 |
| **PRt** | [9.374e-07, 1.166e-06, 3.644e-07] | [-0.31, 1.68] | 783/1369/2638/5172 | 297/594/1187/2374 | 9.374e-07, 1.166e-06, 3.644e-07 / -0.31, 1.68 |

PASS check (finest-pair p >= 1.8): X9t 0.95 FAIL; PRt 1.68 FAIL (first pair non-monotone, as in raycheap6).

## Differences (cmpbin.py: per-variable bitwise / max abs; plus max|dE|/max|E|)

| pair | c0.15 | c0.01875 |
|---|---|---|
| RB vs PR | m1 bitwise 2/4 (E bitwise; F1 3.9e-9, F3 1.9e-9 rel.), max abs 3.7e-9; m1_lat 12/12; dE = 0 | m1 4/4, m1_lat 12/12 bitwise; dE = 0 |
| PRg vs PR | identical pattern to RB vs PR; dE = 0 | m1 4/4, m1_lat 12/12 bitwise; dE = 0 |
| RB vs PRg | m1 4/4, m1_lat 12/12 bitwise | - |
| X9 vs PR (info) | max|dE|/max|E| 1.0e-3 | 1.7e-4 |
| PRt vs PR (info) | 1.0e-3 | 6.7e-5 |
| X9t vs X9 (info) | 6.7e-5 | 6.7e-5 |

RB vs PR and PRg vs PR are also m1-bitwise at c0.075 and c0.0375.

## Per-run lines (rc, FATAL, NOT ADMISSIBLE, DT_RETRY, wall, builds | hesdirk2 counters)
```
gfb_X9_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 665 s builds=297 | hesdirk2: stage steps=296 backward-Euler steps=1 stage fallbacks=0
gfb_X9_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1169 s builds=594 | hesdirk2: stage steps=593 backward-Euler steps=1 stage fallbacks=0
gfb_X9_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 2233 s builds=1187 | hesdirk2: stage steps=1186 backward-Euler steps=1 stage fallbacks=0
gfb_X9_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 4297 s builds=2374 | hesdirk2: stage steps=2373 backward-Euler steps=1 stage fallbacks=0
gfb_PR_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 654 s builds=297 | hesdirk2: stage steps=296 backward-Euler steps=1 stage fallbacks=0
gfb_PR_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1154 s builds=594 | hesdirk2: stage steps=593 backward-Euler steps=1 stage fallbacks=0
gfb_PR_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 2238 s builds=1187 | hesdirk2: stage steps=1186 backward-Euler steps=1 stage fallbacks=0
gfb_PR_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 4276 s builds=2374 | hesdirk2: stage steps=2373 backward-Euler steps=1 stage fallbacks=0
gfb_RB_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 832 s builds=593 | hesdirk2: stage steps=296 backward-Euler steps=1 stage fallbacks=0
gfb_RB_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1503 s builds=1187 | hesdirk2: stage steps=593 backward-Euler steps=1 stage fallbacks=0
gfb_RB_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 2905 s builds=2373 | hesdirk2: stage steps=1186 backward-Euler steps=1 stage fallbacks=0
gfb_RB_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 5658 s builds=4747 | hesdirk2: stage steps=2373 backward-Euler steps=1 stage fallbacks=0
gfb_PRg_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 655 s builds=297 | hesdirk2: stage steps=296 backward-Euler steps=1 stage fallbacks=0
gfb_PRg_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1157 s builds=594 | hesdirk2: stage steps=593 backward-Euler steps=1 stage fallbacks=0
gfb_PRg_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 2213 s builds=1187 | hesdirk2: stage steps=1186 backward-Euler steps=1 stage fallbacks=0
gfb_PRg_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 4274 s builds=2374 | hesdirk2: stage steps=2373 backward-Euler steps=1 stage fallbacks=0
gfb_X9t_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 788 s builds=297 | hesdirk2: stage steps=296 backward-Euler steps=1 stage fallbacks=0
gfb_X9t_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1376 s builds=594 | hesdirk2: stage steps=593 backward-Euler steps=1 stage fallbacks=0
gfb_X9t_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 2651 s builds=1187 | hesdirk2: stage steps=1186 backward-Euler steps=1 stage fallbacks=0
gfb_X9t_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 5143 s builds=2374 | hesdirk2: stage steps=2373 backward-Euler steps=1 stage fallbacks=0
gfb_PRt_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 783 s builds=297 | hesdirk2: stage steps=296 backward-Euler steps=1 stage fallbacks=0
gfb_PRt_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1369 s builds=594 | hesdirk2: stage steps=593 backward-Euler steps=1 stage fallbacks=0
gfb_PRt_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 2638 s builds=1187 | hesdirk2: stage steps=1186 backward-Euler steps=1 stage fallbacks=0
gfb_PRt_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 5172 s builds=2374 | hesdirk2: stage steps=2373 backward-Euler steps=1 stage fallbacks=0
```

## Files
/orion/ptmp/jinma/raycheap7_1109: rc7.sbatch, rc7.out, ana.sh, rc7_ana.out, cmp.py, runs/*/run.log + wall.txt.
bins deleted after the comparisons; the raycheap6 binary is kept.

