# NOTE 2026-11-09 ORION -> viper: raycheap11 -- stage order of the exchange (F1 E1 N1 T1 R1 R2 R3 on Ar4)

Answer to `TASK-2026-11-09-orion-raycheap11.md`. No code changes.

**Summary (finest-pair p; t_end 10)**
- **T1 (coupling = false, transport alone, predict): p 0.35 / 0.51, NOT 2.** Clean run (no fallbacks). So in this deck
  even the transport with the modulated inner flux and the per-step VET update is not 2nd order without the exchange.
  This breaks the TASK's premise that the exchange is the only remaining term.
- **F1 (frozen closure) p 0.18 / 0.36 and R1 (TR-BDF2 + frozen closure) p -0.08 / -0.13: neither converges.** R1's d is
  flat (3.1e-4, 3.3e-4, 3.6e-4) and F1's nearly so. A frozen closure is built ONCE at the first step, and that
  build differs between levels. This is consistent with viper's raycheap6 point: with `vet_source = relaxed_esc`
  (the default here, log line present) the formal-solution source is relaxed over dt0, so the single build depends
  on the level's dt0 and leaves a level-dependent, non-vanishing difference. Not verified further. So F1 vs R1 cannot
  separate the stage order: R1 vs F1 differ by only 6.5e-3 (c0.15) / 6.0e-5 (c0.01875) of max E.
- **R2 (TR-BDF2, extrap2 update, hydro skipped): p 0.31 / 0.59; R3 (TR-BDF2, full coupled problem): p 0.59 / 0.66.**
  TR-BDF2 does not restore 2nd order (A-X9 with hesdirk2 was 0.47 / 0.89; raycheap8 X9 T10 with hydro 0.66 / 0.85).
  Both have 1-2 stage-2 NOT ADMISSIBLE -> backward-Euler redos at cycles 1 (and 3) on the three coarser levels, none at
  the finest, so the coarse pairs are slightly spoiled (BE steps 3-6 vs 1).
- **N1 (opac_update + gas_newton + opac_newton, predict): identical d and p to raycheap10 A-PR** (9.100e-4, 6.215e-4,
  2.971e-4; p 0.55 / 1.06): the Newton-opacity set runs (B2 done properly, accepted) but changes nothing visible
  with this constant-opacity deck.
- **E1 (eddington + implicit_hr_model): spoiled.** Repeated stage-1 NOT ADMISSIBLE (5-19 per level, at cycles
  5, 8, 13, 22, 39, 56, then every 17) and `time2_auto_be = 16` turn almost every step into backward Euler (BE steps
  61 / 66 / 134 / 269 vs 1-2 hesdirk2 stage steps); c0.15 also has 5 DT_RETRY. p -1.84 / 0.82 is a BE result.
- All 28 runs rc = 0, FATAL 0; no key was refused.

## Binary and decks
- MOD `addb642d0` rebuilt (gcc/13, openmpi/4.1, Release, MPI, `-ffp-contract=off`); kept at
  /orion/ptmp/jinma/raycheap11_1109/athena_addb642d0. The log confirms `time2_tableau = trbdf2: TR-BDF2 implicit stages
  (stage order 2), RK2 hydro with c2 = 2 - sqrt(2)`.
- Base: raycheap8 Ar4 deck (full_decks/evoAr4_c*), t_end via `time/tlim=10`. Absent keys were written into deck copies;
  keys present in the deck (time2_vet_col, closure, coupling) were changed in the copies. All arms carry
  `dbg_hydro_skip_t2 = true` except R3.
- One job 228063 (7 exclusive nodes, all 28 runs at once, 14 min).

## conv.py (cfl 0.15 / 0.075 / 0.0375 / 0.01875)

| arm | d | p | wall (s) | NOT ADMISSIBLE | hesdirk2 stage/BE/fallback | builds |
|---|---|---|---|---|---|---|
| F1 lag + every 1e6 (frozen) | [4.912e-04, 4.335e-04, 3.381e-04] | [0.18, 0.36] | 184/267/365/534 | 0/0/0/0 | 33/1/0 67/1/0 135/1/0 270/1/0 | 1/1/1/1 |
| E1 eddington + hr_model | [8.892e-04, 3.193e-03, 1.808e-03] | [-1.84, 0.82] | 166/61/121/294 | 7/7/11/19 | 1/61/7 2/66/7 2/134/11 2/269/19 | -/-/-/- |
| N1 Newton opacity set + predict | [9.100e-04, 6.215e-04, 2.971e-04] | [0.55, 1.06] | 217/315/454/682 | 0/0/0/0 | 33/1/0 67/1/0 135/1/0 270/1/0 | 34/68/136/271 |
| **T1** coupling=false + predict | [2.148e-03, 1.683e-03, 1.184e-03] | [0.35, 0.51] | 161/237/362/543 | 0/0/0/0 | 33/1/0 67/1/0 135/1/0 270/1/0 | 34/68/136/271 |
| **R1** trbdf2 + frozen | [3.122e-04, 3.292e-04, 3.595e-04] | [-0.08, -0.13] | 183/266/377/537 | 1/1/2/0 | 31/3/1 65/3/1 130/6/2 270/1/0 | 1/1/1/1 |
| R2 trbdf2 (extrap2) | [9.386e-04, 7.577e-04, 5.046e-04] | [0.31, 0.59] | 205/308/454/677 | 1/1/2/0 | 31/3/1 65/3/1 130/6/2 270/1/0 | 35/69/138/271 |
| **R3** trbdf2, hydro ON | [9.743e-04, 6.494e-04, 4.110e-04] | [0.59, 0.66] | 227/342/500/811 | 1/2/2/0 | 36/3/1 72/6/2 150/6/2 311/1/0 | 40/80/158/312 |

References (raycheap10, hesdirk2, hydro skipped): A-PR p 0.55 / 1.06, A-X9 0.47 / 0.89.

Final-E differences (max|dE|/max|E|, c0.15 / c0.01875): R1 vs F1 6.5e-3 / 6.0e-5; R3 vs R2 (hydro on vs skipped) 0.119 / 0.120.

NOT ADMISSIBLE positions (stage@cycle):
- R1, R2: 2@1 (c0.15, c0.075), 2@1, 2@3 (c0.0375), none at c0.01875; R3: 2@1 (c0.15), 2@1, 2@3 (c0.075, c0.0375), none at c0.01875.
- E1: c0.15 2@2, 1@4, 7, 12, 21, 38, 55; c0.01875 2@3, 1@5, 8, 13, 22, 39, 56, 73, ... every 17 to 260.
  E1 finest wall split: clean stage steps 3.7 s, failed stage attempts 71.7 s, their BE redo 16.4 s, auto_be steps 194 s.

Suggested next check (not run): T1 with the closure frozen (coupling false + lag + vet_col_every 1e6) and with
`vet_source = gas`, to see whether the frozen-closure non-convergence is the dt0-dependent source relaxation.

## Per-run lines (rc, FATAL, NOT ADMISSIBLE, DT_RETRY, wall, builds | hesdirk2 counters)
```
F1_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 184 s builds=1 | hesdirk2: stage steps=33 backward-Euler steps=1 stage fallbacks=0
F1_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 267 s builds=1 | hesdirk2: stage steps=67 backward-Euler steps=1 stage fallbacks=0
F1_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 365 s builds=1 | hesdirk2: stage steps=135 backward-Euler steps=1 stage fallbacks=0
F1_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 534 s builds=1 | hesdirk2: stage steps=270 backward-Euler steps=1 stage fallbacks=0
E1_c0.15 rc=0 FATAL=0 NOTADM=7 DTRETRY=5 wall 166 s builds= | hesdirk2: stage steps=1 backward-Euler steps=61 stage fallbacks=7
E1_c0.075 rc=0 FATAL=0 NOTADM=7 DTRETRY=0 wall 61 s builds= | hesdirk2: stage steps=2 backward-Euler steps=66 stage fallbacks=7
E1_c0.0375 rc=0 FATAL=0 NOTADM=11 DTRETRY=0 wall 121 s builds= | hesdirk2: stage steps=2 backward-Euler steps=134 stage fallbacks=11
E1_c0.01875 rc=0 FATAL=0 NOTADM=19 DTRETRY=0 wall 294 s builds= | hesdirk2: stage steps=2 backward-Euler steps=269 stage fallbacks=19
N1_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 217 s builds=34 | hesdirk2: stage steps=33 backward-Euler steps=1 stage fallbacks=0
N1_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 315 s builds=68 | hesdirk2: stage steps=67 backward-Euler steps=1 stage fallbacks=0
N1_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 454 s builds=136 | hesdirk2: stage steps=135 backward-Euler steps=1 stage fallbacks=0
N1_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 682 s builds=271 | hesdirk2: stage steps=270 backward-Euler steps=1 stage fallbacks=0
T1_c0.15 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 161 s builds=34 | hesdirk2: stage steps=33 backward-Euler steps=1 stage fallbacks=0
T1_c0.075 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 237 s builds=68 | hesdirk2: stage steps=67 backward-Euler steps=1 stage fallbacks=0
T1_c0.0375 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 362 s builds=136 | hesdirk2: stage steps=135 backward-Euler steps=1 stage fallbacks=0
T1_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 543 s builds=271 | hesdirk2: stage steps=270 backward-Euler steps=1 stage fallbacks=0
R1_c0.15 rc=0 FATAL=0 NOTADM=1 DTRETRY=0 wall 183 s builds=1 | hesdirk2: stage steps=31 backward-Euler steps=3 stage fallbacks=1
R1_c0.075 rc=0 FATAL=0 NOTADM=1 DTRETRY=0 wall 266 s builds=1 | hesdirk2: stage steps=65 backward-Euler steps=3 stage fallbacks=1
R1_c0.0375 rc=0 FATAL=0 NOTADM=2 DTRETRY=0 wall 377 s builds=1 | hesdirk2: stage steps=130 backward-Euler steps=6 stage fallbacks=2
R1_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 537 s builds=1 | hesdirk2: stage steps=270 backward-Euler steps=1 stage fallbacks=0
R2_c0.15 rc=0 FATAL=0 NOTADM=1 DTRETRY=0 wall 205 s builds=35 | hesdirk2: stage steps=31 backward-Euler steps=3 stage fallbacks=1
R2_c0.075 rc=0 FATAL=0 NOTADM=1 DTRETRY=0 wall 308 s builds=69 | hesdirk2: stage steps=65 backward-Euler steps=3 stage fallbacks=1
R2_c0.0375 rc=0 FATAL=0 NOTADM=2 DTRETRY=0 wall 454 s builds=138 | hesdirk2: stage steps=130 backward-Euler steps=6 stage fallbacks=2
R2_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 677 s builds=271 | hesdirk2: stage steps=270 backward-Euler steps=1 stage fallbacks=0
R3_c0.15 rc=0 FATAL=0 NOTADM=1 DTRETRY=0 wall 227 s builds=40 | hesdirk2: stage steps=36 backward-Euler steps=3 stage fallbacks=1
R3_c0.075 rc=0 FATAL=0 NOTADM=2 DTRETRY=0 wall 342 s builds=80 | hesdirk2: stage steps=72 backward-Euler steps=6 stage fallbacks=2
R3_c0.0375 rc=0 FATAL=0 NOTADM=2 DTRETRY=0 wall 500 s builds=158 | hesdirk2: stage steps=150 backward-Euler steps=6 stage fallbacks=2
R3_c0.01875 rc=0 FATAL=0 NOTADM=0 DTRETRY=0 wall 811 s builds=312 | hesdirk2: stage steps=311 backward-Euler steps=1 stage fallbacks=0
```
E1 has no `closure = vet_col` line (eddington), so builds are empty there.

## Files
/orion/ptmp/jinma/raycheap11_1109: decks/ (F1 E1 N1 T1 R1 R2 R3 per level), rc11.sbatch, ana.sh, rc11_ana.out,
runs/*/run.log + wall.txt, athena_addb642d0 (kept). bins, rst, sources and the build dir deleted.
