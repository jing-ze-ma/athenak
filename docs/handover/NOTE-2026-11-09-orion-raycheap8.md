# NOTE 2026-11-09 ORION -> viper: raycheap8 -- evolving gas: rebuild vs predict, and the VET update order

Answer to `TASK-2026-11-09-orion-raycheap8.md`. No code changes.

**Summary**
- **Part 1: the stage-2 rebuild DOES reach the operator.** On the HG cs small deck (evolving gas), predict and rebuild are
  bitwise equal after cycle 1 (the backward-Euler first step) and differ from cycle 2 on. After 20 cycles:
  m1_e 2.3e-9, F1 3.4e-4, F2/F3 0.95 (of max|F|), eint 7.6e-8, dens 6.3e-7, m1_vet D 6-8e-2 relative (max|d|/max|v|).
- **Part 2: none of the deck candidates A, B, C runs** (nor the TASK's adjustments e_out 3, gas_rho 1e-4, period 400,
  modulation off): every one dies at cycle 1 with `eint<=0` in ALL 442368 cells (DT_RETRY x8 then FATAL).
  Cause as far as I can tell from the inputs: eint/E is 1.5e-10 (A), 1.5e-7 (B), 8e-10 (C) at T = (E/arad)^(1/4) = 1e-4 (C: 5.6e-5) with
  arad 1e16, so any O(1e-10) work/exchange makes eint negative. A deck with eint ~ E runs cleanly; used **"Ar4" = A with
  gas_rho 1e4, kappa_p/e/f 1e-4, clump_rho 5e5** (ambient chi 1, clump chi 50 as in A; eint/E = 1.5).
- **PASS: none.** With evolving gas every arm is 1st order, including the production X9 and the rebuild RB:
  finest-pair p 0.84-1.00 at t_end 10 and 20. PRt 0.84 (T10) / 0.99 (T20); X9t 0.85 / 0.93.
- **RB vs PR now differ**: 1.0e-4 (T10) and 4.1e-4 (T20) of max E at the finest level; PRg vs PR 8.3e-5 / 2.0e-4.
  These are of the same size as the finest-pair d (2.4e-4 T10, 1.5e-4 T20), but neither arm changes the order.
  So in this regime there is an O(dt) term common to all arms that the VET update choice does not touch.
  (Not identified here; with gas_hold off the radiation-gas-hydro coupling is a candidate.)
- All 48 full-set runs: rc=0, FATAL 0, NOT ADMISSIBLE 0, DT_RETRY 0, stage fallbacks 0; 1 backward-Euler step each.

## Binary
- MOD `4f227b777`: `git diff --stat f872899da 4f227b777 -- src CMakeLists.txt` is EMPTY; reused
  /orion/ptmp/jinma/raycheap6_1109/athena_f872899da (PROBLEM none) for Part 2.
- Part 1 needs PROBLEM=he_star_m1, so it used the he_star_m1 binary of the hglagt gate,
  /orion/ptmp/jinma/hglagt_1109/wt_new/build/src/athena (`7fb7af05a`, hgpos-1109c).

## Part 1 (HG cs small, 4 ranks, 20 cycles)
Deck = /orion/ptmp/jinma/hglagt_1109/deck/on.athinput (hgF4w production recipe; `time2_vet_col = extrap2` in the deck),
copies with that line set to predict / rebuild (sed). Same srun as the hglagt job (`-n 4 --mem=48G`), outputs
`output1/dt=1e-30 output2/dt=1e-30 output3/dt=1e-30 output5/dt=1e30 output6/dt=1e-30 time/nlim=20`.
Both rc=0 (1063 s / 1103 s), t = 2857.143 at cycle 20, no NOT ADMISSIBLE / DT_RETRY; predicted-build counter
`limited cell-builds` 3.407e5 (predict) vs 2.732e5 (rebuild) of 2.055e7.

| cycle | m1_e | m1_f1 | eint |
|---|---|---|---|
| 0, 1 | bitwise | bitwise | bitwise |
| 2 | 6.7e-10 | 4.2e-4 | 2.7e-8 |
| 10 | 1.2e-9 | 4.2e-4 | 4.0e-8 |
| 20 | 2.3e-9 | 3.4e-4 | 7.6e-8 |

(max|predict - rebuild| / max|predict|.) Final cycle 20, all variables, none bitwise: dens 6.3e-7, vel 8e-2 ... 2e-1,
eint 7.6e-8; m1_e 2.3e-9, f1 3.4e-4, f2 0.95, f3 0.94; m1_vet fk/d11 6.3e-2, d22/d33 8.4e-2.

## Part 2 pre-check log (c0.15, PR and RB each; all run concurrently)
| deck | result |
|---|---|
| A, B, C (as given) | rc 1: FATAL `dt_retry: cycle 1 failed 9 times (dt_retry_max 8)`, all retries `reason=eint<=0`; hesdirk2 stage 1 NOT POSITIVE in all 442368 cells -> 4 posfix rounds -> NOT ADMISSIBLE -> BE; A also cuts dt 3x at cycle 0; NC-local max up to 1.9e69 |
| A e_out 3 / gas_rho 1e-4 / mod period 400 / mod amp 0 | same failure (e_out 3: NC-local 1.2e103) |
| A, B, C without the clump; A clump_rho = ambient | same failure (so not the clump) |
| A arad 1 (T = 1 = c: hydro dt drops 100x), with and without clump | same failure |
| Ar2: A + gas_rho 1e2, kappa 1e-2, clump_rho 5e3 (eint/E 1.5e-2), with and without clump | ran cleanly to t = 13 (0 NOT ADMISSIBLE, RB != PR) but dt fell from 0.20 to 0.10; cancelled when Ar4 had passed |
| **Ar4: A + gas_rho 1e4, kappa 1e-4, clump_rho 5e5** (eint/E 1.5), with and without clump | **rc 0 at tlim 20, 79 cycles, 0 NOT ADMISSIBLE / DT_RETRY, dt steady 0.25; RB != PR (1.4e-3 of max E with clump, 1.0e-5 without)** -> used |

Side observation (pgen, not changed): with `gas_teq = true` IEN is set from the ambient `gas_rho` after the clump has set
IDN, so the clump starts at T_eq x gas_rho/clump_rho (rad_m1_tests2.cpp ~1309). Removing the clump did not fix A/B/C.

## Part 2 conv.py (Ar4; cfl 0.15 / 0.075 / 0.0375 / 0.01875)

| arm | T | d | p | wall (s) | tensor builds |
|---|---|---|---|---|---|
| X9 | 10 | [1.155e-03, 7.322e-04, 4.059e-04] | [0.66, 0.85] | 223/330/493/832 | 39/78/156/312 |
| X9 | 20 | [1.038e-03, 5.656e-04, 2.957e-04] | [0.88, 0.94] | 407/602/1030/1781 | 79/158/315/630 |
| PR | 10 | [8.178e-04, 4.472e-04, 2.443e-04] | [0.87, 0.87] | 227/329/484/812 | 39/78/156/312 |
| PR | 20 | [6.175e-04, 3.055e-04, 1.534e-04] | [1.02, 0.99] | 410/590/1011/1753 | 79/158/315/630 |
| RB | 10 | [8.251e-04, 4.482e-04, 2.447e-04] | [0.88, 0.87] | 236/365/587/970 | 77/155/311/623 |
| RB | 20 | [6.296e-04, 3.070e-04, 1.537e-04] | [1.04, 1.00] | 463/661/1184/2049 | 157/315/629/1259 |
| PRg | 10 | [6.837e-04, 4.053e-04, 2.408e-04] | [0.75, 0.75] | 221/326/483/817 | 39/78/156/312 |
| PRg | 20 | [5.956e-04, 3.088e-04, 1.555e-04] | [0.95, 0.99] | 404/592/1019/1766 | 79/158/315/630 |
| X9t | 10 | [1.154e-03, 7.282e-04, 4.033e-04] | [0.66, 0.85] | 289/453/684/1047 | 40/79/157/314 |
| X9t | 20 | [1.035e-03, 5.657e-04, 2.959e-04] | [0.87, 0.93] | 569/793/1354/2241 | 80/159/317/634 |
| PRt | 10 | [8.101e-04, 4.450e-04, 2.478e-04] | [0.86, 0.84] | 292/450/678/1058 | 40/79/157/314 |
| PRt | 20 | [6.178e-04, 3.058e-04, 1.536e-04] | [1.01, 0.99] | 575/783/1355/2252 | 80/159/317/634 |

PASS check (finest pair >= 1.8): PRt 0.84 (T10) / 0.99 (T20) FAIL; X9t 0.85 / 0.93 FAIL.

## Differences (max|dE|/max|E| from cmp.py; cmpbin.py: no variable of m1 bitwise in any pair)

| pair | T10 c0.15 | T10 c0.01875 | T20 c0.15 | T20 c0.01875 |
|---|---|---|---|---|
| RB vs PR | 1.2e-3 | **1.0e-4** | 1.4e-3 | **4.1e-4** |
| PRg vs PR | 2.6e-3 | **8.3e-5** | 8.2e-4 | **2.0e-4** |
| X9 vs PR | 1.8e-3 | 1.2e-3 | 4.3e-3 | 7.4e-4 |
| PRt vs PR | 3.3e-2 | 3.3e-2 | 2.9e-2 | 2.9e-2 |

(The tangential arms take 1-4 more cycles: X9t/PRt 634 vs 630 at T20 finest.)

## Per-run lines (rc, FATAL, NOT ADMISSIBLE, DT_RETRY, wall, builds | hesdirk2)
```
X9_T10_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 223 s builds=39 | hesdirk2: stage steps=38 backward-Euler steps=1 stage fallbacks=0
X9_T10_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 330 s builds=78 | hesdirk2: stage steps=77 backward-Euler steps=1 stage fallbacks=0
X9_T10_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 493 s builds=156 | hesdirk2: stage steps=155 backward-Euler steps=1 stage fallbacks=0
X9_T10_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 832 s builds=312 | hesdirk2: stage steps=311 backward-Euler steps=1 stage fallbacks=0
X9_T20_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 407 s builds=79 | hesdirk2: stage steps=78 backward-Euler steps=1 stage fallbacks=0
X9_T20_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 602 s builds=158 | hesdirk2: stage steps=157 backward-Euler steps=1 stage fallbacks=0
X9_T20_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1030 s builds=315 | hesdirk2: stage steps=314 backward-Euler steps=1 stage fallbacks=0
X9_T20_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1781 s builds=630 | hesdirk2: stage steps=629 backward-Euler steps=1 stage fallbacks=0
PR_T10_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 227 s builds=39 | hesdirk2: stage steps=38 backward-Euler steps=1 stage fallbacks=0
PR_T10_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 329 s builds=78 | hesdirk2: stage steps=77 backward-Euler steps=1 stage fallbacks=0
PR_T10_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 484 s builds=156 | hesdirk2: stage steps=155 backward-Euler steps=1 stage fallbacks=0
PR_T10_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 812 s builds=312 | hesdirk2: stage steps=311 backward-Euler steps=1 stage fallbacks=0
PR_T20_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 410 s builds=79 | hesdirk2: stage steps=78 backward-Euler steps=1 stage fallbacks=0
PR_T20_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 590 s builds=158 | hesdirk2: stage steps=157 backward-Euler steps=1 stage fallbacks=0
PR_T20_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1011 s builds=315 | hesdirk2: stage steps=314 backward-Euler steps=1 stage fallbacks=0
PR_T20_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1753 s builds=630 | hesdirk2: stage steps=629 backward-Euler steps=1 stage fallbacks=0
RB_T10_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 236 s builds=77 | hesdirk2: stage steps=38 backward-Euler steps=1 stage fallbacks=0
RB_T10_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 365 s builds=155 | hesdirk2: stage steps=77 backward-Euler steps=1 stage fallbacks=0
RB_T10_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 587 s builds=311 | hesdirk2: stage steps=155 backward-Euler steps=1 stage fallbacks=0
RB_T10_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 970 s builds=623 | hesdirk2: stage steps=311 backward-Euler steps=1 stage fallbacks=0
RB_T20_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 463 s builds=157 | hesdirk2: stage steps=78 backward-Euler steps=1 stage fallbacks=0
RB_T20_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 661 s builds=315 | hesdirk2: stage steps=157 backward-Euler steps=1 stage fallbacks=0
RB_T20_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1184 s builds=629 | hesdirk2: stage steps=314 backward-Euler steps=1 stage fallbacks=0
RB_T20_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 2049 s builds=1259 | hesdirk2: stage steps=629 backward-Euler steps=1 stage fallbacks=0
PRg_T10_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 221 s builds=39 | hesdirk2: stage steps=38 backward-Euler steps=1 stage fallbacks=0
PRg_T10_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 326 s builds=78 | hesdirk2: stage steps=77 backward-Euler steps=1 stage fallbacks=0
PRg_T10_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 483 s builds=156 | hesdirk2: stage steps=155 backward-Euler steps=1 stage fallbacks=0
PRg_T10_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 817 s builds=312 | hesdirk2: stage steps=311 backward-Euler steps=1 stage fallbacks=0
PRg_T20_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 404 s builds=79 | hesdirk2: stage steps=78 backward-Euler steps=1 stage fallbacks=0
PRg_T20_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 592 s builds=158 | hesdirk2: stage steps=157 backward-Euler steps=1 stage fallbacks=0
PRg_T20_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1019 s builds=315 | hesdirk2: stage steps=314 backward-Euler steps=1 stage fallbacks=0
PRg_T20_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1766 s builds=630 | hesdirk2: stage steps=629 backward-Euler steps=1 stage fallbacks=0
X9t_T10_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 289 s builds=40 | hesdirk2: stage steps=39 backward-Euler steps=1 stage fallbacks=0
X9t_T10_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 453 s builds=79 | hesdirk2: stage steps=78 backward-Euler steps=1 stage fallbacks=0
X9t_T10_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 684 s builds=157 | hesdirk2: stage steps=156 backward-Euler steps=1 stage fallbacks=0
X9t_T10_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1047 s builds=314 | hesdirk2: stage steps=313 backward-Euler steps=1 stage fallbacks=0
X9t_T20_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 569 s builds=80 | hesdirk2: stage steps=79 backward-Euler steps=1 stage fallbacks=0
X9t_T20_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 793 s builds=159 | hesdirk2: stage steps=158 backward-Euler steps=1 stage fallbacks=0
X9t_T20_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1354 s builds=317 | hesdirk2: stage steps=316 backward-Euler steps=1 stage fallbacks=0
X9t_T20_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 2241 s builds=634 | hesdirk2: stage steps=633 backward-Euler steps=1 stage fallbacks=0
PRt_T10_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 292 s builds=40 | hesdirk2: stage steps=39 backward-Euler steps=1 stage fallbacks=0
PRt_T10_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 450 s builds=79 | hesdirk2: stage steps=78 backward-Euler steps=1 stage fallbacks=0
PRt_T10_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 678 s builds=157 | hesdirk2: stage steps=156 backward-Euler steps=1 stage fallbacks=0
PRt_T10_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1058 s builds=314 | hesdirk2: stage steps=313 backward-Euler steps=1 stage fallbacks=0
PRt_T20_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 575 s builds=80 | hesdirk2: stage steps=79 backward-Euler steps=1 stage fallbacks=0
PRt_T20_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 783 s builds=159 | hesdirk2: stage steps=158 backward-Euler steps=1 stage fallbacks=0
PRt_T20_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 1355 s builds=317 | hesdirk2: stage steps=316 backward-Euler steps=1 stage fallbacks=0
PRt_T20_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 2252 s builds=634 | hesdirk2: stage steps=633 backward-Euler steps=1 stage fallbacks=0
```
All stage fallbacks and accepted-with-fallback counts are 0.

## Jobs and files
Pre-check + Part 1: jobs 228038 (2 nodes), 228039 (3 nodes), 228041 (2 nodes) concurrently; full set 228043
(6 exclusive nodes, 24 concurrent 24-rank steps, longest first, 38 min; critical path X9t/PRt_T20_c0.01875 ~2250 s).
/orion/ptmp/jinma/raycheap8_1109: pre/, p1/, runs/ (logs + wall.txt), decks/ and full_decks/ (the adjusted decks),
*.sbatch, ana.sh, rc8_ana.out, cmp.py, p1cmp.py. bins and rst deleted.

