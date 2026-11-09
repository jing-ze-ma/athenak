# TASK 2026-11-09 viper -> ORION: raycheap7 -- VET update order in the production-like regime (gas_feedback on)

Answer with `docs/handover/NOTE-2026-11-xx-orion-raycheap7.md`. Follow-up to raycheap6: there the gas was held with
gas_feedback = false, which zeroes the VET source relaxation (clq = 0, "held gas: T* = T_b", rad_m1_vetlat.cpp ~650)
and makes the rebuild's stage-2 inputs identical to predict's (RB == PR bitwise by construction). Here
gas_feedback = TRUE (the gas is still re-imposed every stage by gas_hold, so the problem stays defined), which turns
the dt-dependent source relaxation on (vet_source = relaxed_esc, the production default) and makes the closure
depend on the stage state.

## Binary and decks
MOD = `fork/raycheap2-1109` head `8532af080` (code = a1481b185, your raycheap6 build: reuse it).
Decks `tests_raycheap/gasfb/gfb_<arm>_c<cfl>.athinput` = raycheap6 gasmod decks + gas_feedback = true;
24 MeshBlocks = 24 ranks, t_end 400, cfl 0.15 / 0.075 / 0.0375 / 0.01875.

| arm | time2_vet_col | other |
|---|---|---|
| X9 | extrap2, x2_fields 9 (production) | |
| PR | predict | |
| RB | rebuild | must now DIFFER from PR (the stage-2 build reaches the operator?) |
| PRg | predict | vet_source = gas (no relaxation: is the relaxation the remaining dt term?) |
| X9t / PRt | as X9 / PR | vet_gd_tangential = true (the full tensor, the user's criterion) |

## Analysis and PASS
Per arm: `python3 tests_raycheap/conv.py <c0.15> <c0.075> <c0.0375> <c0.01875>`; report d, p, wall, the
`closure = vet_col: N tensor builds` counts and `time_scheme=hesdirk2:` counts; and max|dE|/max|E| of RB vs PR and
PRg vs PR at the finest level (`python3 tests_raycheap/cmpbin.py`). PASS (user's VET criterion): finest-pair p >= 1.8
for PRt (and/or X9t). If gas_feedback = true makes any arm stop (FATAL / NOT ADMISSIBLE storms), report and skip it.

## NOTE back
conv.py lines, counts, walls, the RB/PR and PRg/PR differences. No code changes.
